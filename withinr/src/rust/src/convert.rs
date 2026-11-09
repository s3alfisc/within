//! Shared R ↔ Rust coercion helpers, error plumbing, and design extraction.
//!
//! Mirrors `crates/within-py/src/convert.rs` plus the design extraction from
//! `crates/within-py/src/api.rs`. Factor labels pass through unchanged: an R
//! label `7` is native label `7`, so coefficient addresses report the labels
//! the caller supplied. Term, column, and `x` indices are 1-based in R.

use extendr_api::{error::Result, prelude::*, AsTypedSlice};
use ndarray::{ArrayView2, ShapeBuilder};

use within::{Design, Effect};

pub(crate) fn err(message: impl Into<String>) -> Error {
    Error::Other(message.into())
}

/// Raise `Err` as an R error, preserving its message.
///
/// extendr's default handling of `Result` return values panics, which surfaces
/// in R as a generic "User function panicked" message. Every `#[extendr]`
/// entry point therefore funnels its fallible body through this helper.
pub(crate) fn or_throw<T>(result: Result<T>) -> T {
    match result {
        Ok(value) => value,
        Err(error) => throw_r_error(error.to_string()),
    }
}

pub(crate) fn usize_to_i32(value: usize, name: &str) -> Result<i32> {
    i32::try_from(value).map_err(|_| err(format!("{name} exceeds R's integer range")))
}

/// 0-based native position → 1-based R index.
pub(crate) fn one_based(value: usize, name: &str) -> Result<i32> {
    usize_to_i32(value, name)?
        .checked_add(1)
        .ok_or_else(|| err(format!("{name} exceeds R's integer range")))
}

pub(crate) fn get_field(obj: &Robj, field: &str) -> Robj {
    obj.dollar(field).unwrap_or_else(|_| nil_value())
}

// ---------------------------------------------------------------------------
// Scalars
// ---------------------------------------------------------------------------

/// A length-one, non-missing integer or double.
pub(crate) fn scalar_f64(obj: &Robj, name: &str) -> Result<f64> {
    if let Some(&[value]) = obj.as_integer_slice() {
        if value != i32::MIN {
            return Ok(f64::from(value));
        }
    } else if let Some(&[value]) = obj.as_real_slice() {
        if !value.is_nan() {
            return Ok(value);
        }
    }
    Err(err(format!("`{name}` must be a single non-missing number")))
}

/// A whole number in `[min, max]`.
fn scalar_whole(obj: &Robj, name: &str, min: f64, max: f64) -> Result<f64> {
    let value = scalar_f64(obj, name)?;
    if !value.is_finite() || value.fract() != 0.0 || value < min || value > max {
        return Err(err(format!(
            "`{name}` must be a whole number in [{min}, {max}]"
        )));
    }
    Ok(value)
}

pub(crate) fn scalar_usize(obj: &Robj, name: &str, min: usize) -> Result<usize> {
    // Keep within R's integer range so the value round-trips back to R exactly.
    scalar_whole(obj, name, min as f64, f64::from(i32::MAX)).map(|v| v as usize)
}

pub(crate) fn scalar_u32(obj: &Robj, name: &str, min: u32) -> Result<u32> {
    scalar_whole(obj, name, f64::from(min), f64::from(i32::MAX)).map(|v| v as u32)
}

/// Seeds are doubles in R; every integer up to 2^53 is exact.
pub(crate) fn scalar_u64(obj: &Robj, name: &str) -> Result<u64> {
    scalar_whole(obj, name, 0.0, 9_007_199_254_740_992.0).map(|v| v as u64)
}

pub(crate) fn optional<T>(obj: &Robj, parse: impl FnOnce(&Robj) -> Result<T>) -> Result<Option<T>> {
    if obj.is_null() {
        Ok(None)
    } else {
        parse(obj).map(Some)
    }
}

pub(crate) fn scalar_bool(obj: &Robj, name: &str) -> Result<bool> {
    obj.as_bool()
        .ok_or_else(|| err(format!("`{name}` must be TRUE or FALSE")))
}

// ---------------------------------------------------------------------------
// Vectors
// ---------------------------------------------------------------------------

/// Borrow a double vector (the R wrappers coerce before calling).
///
/// The slice lives as long as the R vector, which the calling R frame keeps
/// alive for the whole `.Call`, so it may outlive this (temporary) `Robj`
/// handle — the same contract as extendr's `as_integer_slice`.
pub(crate) fn real_slice<'a>(obj: &Robj, name: &str) -> Result<&'a [f64]> {
    if obj.rtype() != Rtype::Doubles {
        return Err(err(format!("`{name}` must be a double vector")));
    }
    <Robj as AsTypedSlice<'a, f64>>::as_typed_slice(obj)
        .ok_or_else(|| err(format!("`{name}` must be a double vector")))
}

/// Validated non-negative factor labels, copied to the native `u32` layout.
pub(crate) fn labels_to_u32(data: &[i32], name: &str) -> Result<Vec<u32>> {
    data.iter()
        .enumerate()
        .map(|(pos, &value)| {
            if value == i32::MIN {
                Err(err(format!(
                    "`{name}` must not contain NA (found at position {})",
                    pos + 1
                )))
            } else if value < 0 {
                Err(err(format!(
                    "`{name}` labels must be non-negative; found {value} at position {}",
                    pos + 1
                )))
            } else {
                Ok(value as u32)
            }
        })
        .collect()
}

pub(crate) fn optional_weights<'a>(weights: &Robj) -> Result<Option<&'a [f64]>> {
    optional(weights, |w| real_slice(w, "weights"))
}

/// `(nrow, ncol)` of a two-dimensional R array.
fn matrix_dims(obj: &Robj) -> Option<(usize, usize)> {
    match obj.get_attrib("dim")?.as_integer_slice()? {
        &[nrow, ncol] if nrow >= 0 && ncol >= 0 => Some((nrow as usize, ncol as usize)),
        _ => None,
    }
}

/// Borrow the columns of a column-major double matrix; R columns are always contiguous.
pub(crate) fn matrix_columns<'a>(matrix: &Robj, name: &str) -> Result<(usize, Vec<&'a [f64]>)> {
    let data = real_slice(matrix, name)?;
    let (nrow, ncol) =
        matrix_dims(matrix).ok_or_else(|| err(format!("`{name}` must be a matrix")))?;
    let columns = (0..ncol).map(|j| &data[j * nrow..(j + 1) * nrow]).collect();
    Ok((nrow, columns))
}

// ---------------------------------------------------------------------------
// Designs
// ---------------------------------------------------------------------------

/// Native columns behind one R `Effect()`: labels are converted, slopes borrowed.
pub(crate) struct EffectColumns<'a> {
    levels: Vec<u32>,
    intercept: bool,
    slopes: Vec<&'a [f64]>,
}

/// Native categories matrix, converted from R's 1-column-per-factor integers.
pub(crate) struct Categories {
    labels: Vec<u32>,
    nrow: usize,
    ncol: usize,
}

impl Categories {
    pub(crate) fn view(&self) -> Result<ArrayView2<'_, u32>> {
        ArrayView2::from_shape((self.nrow, self.ncol).f(), &self.labels)
            .map_err(|e| err(e.to_string()))
    }

    pub(crate) fn nrow(&self) -> usize {
        self.nrow
    }
}

/// A persistent, categorical, or effect-based design supplied from R.
pub(crate) enum DesignSource<'a> {
    /// A cheap clone of an existing persistent `Design()` handle.
    Persistent(Design<'static>),
    Categories(Categories),
    Effects(Vec<EffectColumns<'a>>),
}

impl<'a> DesignSource<'a> {
    pub(crate) fn effects<'c>(columns: &'c [EffectColumns<'a>]) -> Result<Vec<Effect<'c>>> {
        columns
            .iter()
            .map(|c| Effect::new(&c.levels, c.intercept, c.slopes.iter().copied()))
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|e| err(e.to_string()))
    }

    /// Observation count, read without building the design.
    pub(crate) fn n_obs(&self) -> Option<usize> {
        match self {
            Self::Persistent(design) => Some(design.n_obs()),
            Self::Categories(categories) => Some(categories.nrow()),
            Self::Effects(columns) => columns.first().map(|c| c.levels.len()),
        }
    }

    /// Build an owned native design (for `Design()` and `Solver()`).
    pub(crate) fn into_owned(self) -> Result<Design<'static>> {
        match self {
            Self::Persistent(design) => Ok(design),
            Self::Categories(categories) => {
                let view = categories.view()?;
                crate::threads::run(|| Design::from_categories(view).map(Design::into_owned))?
                    .map_err(|e| err(e.to_string()))
            }
            Self::Effects(columns) => {
                let effects = Self::effects(&columns)?;
                crate::threads::run(|| Design::new(effects).map(Design::into_owned))?
                    .map_err(|e| err(e.to_string()))
            }
        }
    }
}

fn parse_effect<'a>(obj: &Robj, term: usize) -> Result<EffectColumns<'a>> {
    if !obj.inherits("within_effect") {
        return Err(err(format!(
            "design term {} must be created by Effect()",
            term + 1
        )));
    }
    let levels = get_field(obj, "levels");
    let levels = levels
        .as_integer_slice()
        .ok_or_else(|| err(format!("design term {} levels must be integer", term + 1)))?;
    let levels = labels_to_u32(levels, &format!("Effect {} levels", term + 1))?;
    let intercept = scalar_bool(&get_field(obj, "intercept"), "intercept")?;
    let slopes = get_field(obj, "slopes");
    let slopes = slopes
        .as_list()
        .ok_or_else(|| err(format!("design term {} slopes must be a list", term + 1)))?
        .values()
        .enumerate()
        .map(|(slope, value)| {
            real_slice(&value, &format!("Effect {} slope {}", term + 1, slope + 1))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(EffectColumns {
        levels,
        intercept,
        slopes,
    })
}

/// Parse a `Design()` handle, an integer category matrix, or a list of `Effect()`s.
pub(crate) fn parse_design<'a>(design: &Robj) -> Result<DesignSource<'a>> {
    if design.inherits("within_design") {
        let ptr = ExternalPtr::<Design<'static>>::try_from(get_field(design, "ptr"))
            .map_err(|_| err("`design` holds an invalid Design() handle"))?;
        return Ok(DesignSource::Persistent(ptr.try_addr()?.clone()));
    }
    if design.is_matrix() {
        let data = design
            .as_integer_slice()
            .ok_or_else(|| err("a design matrix must hold integer factor labels"))?;
        let (nrow, ncol) =
            matrix_dims(design).ok_or_else(|| err("a design matrix must be two-dimensional"))?;
        return Ok(DesignSource::Categories(Categories {
            labels: labels_to_u32(data, "design")?,
            nrow,
            ncol,
        }));
    }
    let terms = design.as_list().ok_or_else(|| {
        err("`design` must be a Design(), an integer matrix, or a list of Effect()")
    })?;
    let columns = terms
        .values()
        .enumerate()
        .map(|(term, obj)| parse_effect(&obj, term))
        .collect::<Result<Vec<_>>>()?;
    Ok(DesignSource::Effects(columns))
}
