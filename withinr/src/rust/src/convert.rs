//! Shared R ↔ Rust coercion helpers and error plumbing for the bridge.
//!
//! R users pass **1-based** integer factor labels. This module validates that
//! all entries are positive and contain no `NA`, then subtracts 1 to produce
//! the caller labels used by the Rust solver.

use extendr_api::{error::Result, prelude::*};
use ndarray::{ArrayView2, ShapeBuilder};

use within::{Design, Effect, IntoDesign};

pub(crate) fn err(message: impl Into<String>) -> Error {
    Error::Other(message.into())
}

/// Raise `Err` as an R error, preserving its message.
///
/// extendr's default handling of `Result` return values unwraps them (a
/// panic), which surfaces in R as a generic "User function panicked" message.
/// Every `#[extendr]` entry point therefore wraps its fallible body in a
/// `?`-friendly closure and funnels the outcome through this helper instead.
pub(crate) fn or_throw<T>(result: Result<T>) -> T {
    match result {
        Ok(value) => value,
        Err(error) => throw_r_error(error.to_string()),
    }
}

pub(crate) fn usize_to_i32(value: usize, name: &str) -> Result<i32> {
    i32::try_from(value).map_err(|_| err(format!("{name} exceeds i32 range")))
}

pub(crate) fn parse_positive_f64(field: &Robj, name: &str) -> Result<f64> {
    let value = if let Some(value) = field.as_real() {
        value
    } else if let Some(value) = field.as_integer() {
        value as f64
    } else {
        return Err(err(format!("{name} must be a positive finite number")));
    };

    if !value.is_finite() || value <= 0.0 {
        return Err(err(format!("{name} must be a positive finite number")));
    }
    Ok(value)
}

pub(crate) fn parse_nonnegative_f64(field: &Robj, name: &str) -> Result<f64> {
    let value = if let Some(value) = field.as_real() {
        value
    } else if let Some(value) = field.as_integer() {
        value as f64
    } else {
        return Err(err(format!("{name} must be a non-negative finite number")));
    };

    if !value.is_finite() || value < 0.0 {
        return Err(err(format!("{name} must be a non-negative finite number")));
    }
    Ok(value)
}

pub(crate) fn parse_nonnegative_integer(field: &Robj, name: &str) -> Result<i64> {
    let value = if let Some(value) = field.as_integer() {
        value as f64
    } else if let Some(value) = field.as_real() {
        value
    } else {
        return Err(err(format!("{name} must be numeric")));
    };

    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 {
        return Err(err(format!("{name} must be a non-negative integer")));
    }
    Ok(value as i64)
}

pub(crate) fn parse_positive_u32(field: &Robj, name: &str) -> Result<u32> {
    let value = parse_nonnegative_integer(field, name)?;
    if value < 1 || value > u32::MAX as i64 {
        return Err(err(format!("{name} must be in 1..={}", u32::MAX)));
    }
    Ok(value as u32)
}

pub(crate) fn parse_optional_u32(field: &Robj, name: &str) -> Result<Option<u32>> {
    if field.is_null() {
        return Ok(None);
    }
    Ok(Some(parse_positive_u32(field, name)?))
}

pub(crate) fn get_field_or_null(obj: &Robj, field: &str) -> Robj {
    obj.dollar(field).unwrap_or_else(|_| nil_value())
}

/// Borrowed weights for the one-shot paths (`NULL` → `None`), avoiding a copy.
pub(crate) fn weights_slice(weights: &Robj) -> Result<Option<&[f64]>> {
    if weights.is_null() {
        return Ok(None);
    }
    weights
        .as_real_slice()
        .map(Some)
        .ok_or_else(|| err("weights must be a numeric vector or NULL"))
}

/// Convert an R integer matrix (1-based, column-major) to a 0-based `u32`
/// buffer preserving R's column-major layout.
pub(crate) fn cast_categories(data: &[i32]) -> Result<Vec<u32>> {
    let mut out = Vec::with_capacity(data.len());
    for (pos, &value) in data.iter().enumerate() {
        if value == i32::MIN {
            return Err(err(format!(
                "categories must not contain NA values (found at position {})",
                pos + 1
            )));
        }
        if value < 1 {
            return Err(err(format!(
                "categories must be >= 1 (1-based); found {value} at position {}",
                pos + 1
            )));
        }
        out.push((value - 1) as u32);
    }
    Ok(out)
}

pub(crate) fn categories_view<'a>(
    categories: &RMatrix<i32>,
    cats_u32: &'a [u32],
) -> Result<ArrayView2<'a, u32>> {
    ArrayView2::from_shape((categories.nrows(), categories.ncols()).f(), cats_u32)
        .map_err(|e| err(e.to_string()))
}

fn categories_design(categories: &RMatrix<i32>) -> Result<Design<'static>> {
    let cats_u32 = cast_categories(categories.data())?;
    let cats = categories_view(categories, &cats_u32)?;
    crate::threads::run(|| cats.into_design().map(Design::into_owned))?
        .map_err(|e| err(e.to_string()))
}

struct OwnedEffect {
    levels: Vec<u32>,
    intercept: bool,
    slopes: Vec<Vec<f64>>,
}

fn parse_effect(obj: Robj, term: usize) -> Result<OwnedEffect> {
    if !obj.inherits("within_effect") {
        return Err(err(format!(
            "design term {} must be created by Effect(...)",
            term + 1
        )));
    }

    let levels_obj = get_field_or_null(&obj, "levels");
    let levels = levels_obj
        .as_integer_slice()
        .ok_or_else(|| err(format!("design term {} levels must be integer", term + 1)))?;
    let levels = cast_categories(levels)?;

    let intercept = get_field_or_null(&obj, "intercept")
        .as_bool()
        .ok_or_else(|| {
            err(format!(
                "design term {} intercept must be TRUE or FALSE",
                term + 1
            ))
        })?;

    let slopes_obj = get_field_or_null(&obj, "slopes");
    let slopes_list = slopes_obj
        .as_list()
        .ok_or_else(|| err(format!("design term {} slopes must be a list", term + 1)))?;
    let slopes = slopes_list
        .values()
        .enumerate()
        .map(|(slope, value)| {
            value.as_real_vector().ok_or_else(|| {
                err(format!(
                    "design term {} slope {} must be numeric",
                    term + 1,
                    slope + 1
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(OwnedEffect {
        levels,
        intercept,
        slopes,
    })
}

/// Parse an R category matrix or list of `Effect()` terms into an owned native
/// design. Owning the columns keeps persistent solver handles independent of
/// the lifetime and mutability of their originating R objects.
pub(crate) fn parse_design(design: &Robj) -> Result<Design<'static>> {
    if design.inherits("within_design") {
        let ptr = ExternalPtr::<Design<'static>>::try_from(get_field_or_null(design, "ptr"))?;
        return Ok(ptr.try_addr()?.clone());
    }
    if design.is_matrix() {
        let categories = RMatrix::<i32>::try_from(design)
            .map_err(|_| err("design matrices must contain integer factor labels"))?;
        return categories_design(&categories);
    }

    let terms = design
        .as_list()
        .ok_or_else(|| err("design must be an integer matrix or a list of Effect() terms"))?;
    let owned = terms
        .values()
        .enumerate()
        .map(|(term, obj)| parse_effect(obj, term))
        .collect::<Result<Vec<_>>>()?;
    let effects = owned
        .iter()
        .map(|effect| {
            Effect::new(
                &effect.levels,
                effect.intercept,
                effect.slopes.iter().map(Vec::as_slice),
            )
        })
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|e| err(e.to_string()))?;

    crate::threads::run(|| Design::new(effects).map(Design::into_owned))?
        .map_err(|e| err(e.to_string()))
}
