//! Shared R ↔ Rust coercion helpers, error plumbing, and design extraction.
//!
//! Mirrors `crates/within-py/src/convert.rs` plus the design extraction from
//! `crates/within-py/src/api.rs`. Factor labels pass through unchanged: an R
//! label `7` is native label `7`, so coefficient addresses report the labels
//! the caller supplied. Term, column, and `x` indices are 1-based in R, and so
//! are the positions named in error and warning messages.

use std::panic::{catch_unwind, resume_unwind, AssertUnwindSafe};

use extendr_api::{error::Result, prelude::*, AsTypedSlice};
use ndarray::{ArrayView2, ShapeBuilder};

use within::{BuildError, BuildWarning, ChannelPair, Design, Effect, SolveError, WithinError};

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Message prefixes the R wrapper (`native()` in `R/errors.R`) turns into the
/// condition classes `within_input_error` and `within_runtime_error`, the
/// counterparts of Python's `ValueError`/`TypeError` and `RuntimeError`.
const INPUT_TAG: &str = "[within:input] ";
const RUNTIME_TAG: &str = "[within:runtime] ";

/// Invalid input: a wrong type, value, or shape.
pub(crate) fn err(message: impl Into<String>) -> Error {
    Error::Other(message.into())
}

/// A failure while running on valid input, such as a diverging subdomain solve.
pub(crate) fn runtime_err(message: impl std::fmt::Display) -> Error {
    Error::Other(format!("{RUNTIME_TAG}{message}"))
}

/// Run an `#[extendr]` body and raise its error, or a panic, as an R error.
///
/// The error leaves by unwinding: extendr's wrapper catches the panic, drops
/// every argument it protected, and only then raises the R error. Raising from
/// inside the body instead (`throw_r_error`) would longjmp over those drops and
/// pin the arguments for the rest of the session. extendr passes the message
/// to `Rf_error` as a *format string*, so `%` is escaped.
pub(crate) fn guarded<T>(body: impl FnOnce() -> Result<T>) -> T {
    let message = match catch_unwind(AssertUnwindSafe(body)) {
        Ok(Ok(value)) => return value,
        Ok(Err(error)) => {
            let message = error.to_string();
            if message.starts_with(RUNTIME_TAG) {
                message
            } else {
                format!("{INPUT_TAG}{message}")
            }
        }
        Err(panic) => {
            let detail = panic
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| panic.downcast_ref::<&str>().map(|s| (*s).to_owned()))
                .unwrap_or_else(|| "unknown panic".to_owned());
            format!("{RUNTIME_TAG}internal error in within: {detail}")
        }
    };
    resume_unwind(Box::new(message.replace('%', "%%")))
}

fn describe_pair(pair: &ChannelPair) -> String {
    format!(
        "term {} column {} and term {} column {}",
        pair.rows.term + 1,
        pair.rows.column + 1,
        pair.cols.term + 1,
        pair.cols.column + 1
    )
}

/// The core's `Display`, with 0-based positions rewritten 1-based for R.
pub(crate) fn describe_build_error(error: &BuildError) -> String {
    match error {
        BuildError::ObservationCountMismatch {
            column,
            expected,
            got,
        } => format!(
            "design column {} has {got} observations, expected {expected}",
            column + 1
        ),
        BuildError::SlopeLengthMismatch {
            slope,
            expected,
            got,
        } => format!("slope {} has {got} values, expected {expected}", slope + 1),
        BuildError::InvalidLoading {
            slope,
            index,
            value,
        } => format!(
            "slope {} value at position {} must be finite, got {value}",
            slope + 1,
            index + 1
        ),
        BuildError::InvalidWeight { index, value } => format!(
            "weight at position {} must be finite and non-negative, got {value}",
            index + 1
        ),
        BuildError::UnscalableComponent { pair } => format!(
            "signed component between {} is not diagonally scalable to SDDM form",
            describe_pair(pair)
        ),
        BuildError::SingularDiagonal { index } => {
            format!("zero diagonal in preconditioner at position {}", index + 1)
        }
        other => other.to_string(),
    }
}

/// `BuildWarning`'s `Display`, with 1-based terms and columns.
pub(crate) fn describe_warning(warning: &BuildWarning) -> String {
    match warning {
        BuildWarning::UnscalableComponent {
            pair,
            iterations,
            violation,
        } => format!(
            "signed component between {}: dominance scaling uncertified after {iterations} \
             iterations (max relative violation {violation:.2e}); deficits clamped, \
             preconditioner quality may degrade",
            describe_pair(pair)
        ),
        other => other.to_string(),
    }
}

/// Classified like Python: build failures and invalid solve input are input
/// errors, every other solve failure is a runtime error.
pub(crate) fn build_err(error: BuildError) -> Error {
    err(describe_build_error(&error))
}

pub(crate) fn solve_err(error: SolveError) -> Error {
    match error {
        SolveError::InvalidInput { .. } => err(error.to_string()),
        other => runtime_err(other),
    }
}

pub(crate) fn within_err(error: WithinError) -> Error {
    match error {
        WithinError::Build(error) => build_err(error),
        WithinError::Solve(error) => solve_err(error),
        other => err(other.to_string()),
    }
}

// ---------------------------------------------------------------------------
// Native handles
// ---------------------------------------------------------------------------

/// The external pointer behind an R handle, rejecting the wrong type and the
/// null pointer a handle becomes after `saveRDS()`/`readRDS()`.
pub(crate) fn handle<T: 'static>(obj: &Robj, what: &str) -> Result<ExternalPtr<T>> {
    match ExternalPtr::<T>::try_from(obj.clone()) {
        Ok(ptr) => Ok(ptr),
        Err(Error::ExpectedExternalNonNullPtr(_)) => Err(err(format!(
            "this {what} is no longer valid: native handles do not survive saveRDS(), \
             serialize(), or a new R session; rebuild it"
        ))),
        Err(_) => Err(err(format!("expected a {what} handle"))),
    }
}

// ---------------------------------------------------------------------------
// Scalars
// ---------------------------------------------------------------------------

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

/// Capped at R's integer range, so the value round-trips back to R exactly.
pub(crate) fn scalar_usize(obj: &Robj, name: &str, min: usize) -> Result<usize> {
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
/// handle — the same contract as extendr's `as_integer_slice`. No R code runs
/// while native work holds it, so the garbage collector cannot free it.
pub(crate) fn real_slice<'a>(obj: &Robj, name: &str) -> Result<&'a [f64]> {
    if obj.rtype() != Rtype::Doubles {
        return Err(err(format!("`{name}` must be a double vector")));
    }
    <Robj as AsTypedSlice<'a, f64>>::as_typed_slice(obj)
        .ok_or_else(|| err(format!("`{name}` must be a double vector")))
}

/// Borrow an integer vector (the R wrappers coerce before calling).
pub(crate) fn integer_slice<'a>(obj: &Robj, name: &str) -> Result<&'a [i32]> {
    obj.as_integer_slice()
        .ok_or_else(|| err(format!("`{name}` must be an integer vector")))
}

/// Reject NA/NaN/Inf with the 1-based position, before the core sees it.
pub(crate) fn check_finite(values: &[f64], name: &str) -> Result<()> {
    match values.iter().position(|v| !v.is_finite()) {
        None => Ok(()),
        Some(pos) => Err(err(format!(
            "`{name}` must be finite; found {} at position {}",
            values[pos],
            pos + 1
        ))),
    }
}

/// Validated non-negative factor labels, copied to the native `u32` layout.
/// With `nrow`, offending entries are reported by row and column.
pub(crate) fn labels_to_u32(data: &[i32], name: &str, nrow: Option<usize>) -> Result<Vec<u32>> {
    let at = |pos: usize| match nrow {
        Some(nrow) if nrow > 0 => format!("row {}, column {}", pos % nrow + 1, pos / nrow + 1),
        _ => format!("position {}", pos + 1),
    };
    data.iter()
        .enumerate()
        .map(|(pos, &value)| {
            if value == i32::MIN {
                Err(err(format!(
                    "`{name}` must not contain NA (found at {})",
                    at(pos)
                )))
            } else if value < 0 {
                Err(err(format!(
                    "`{name}` labels must be non-negative; found {value} at {}",
                    at(pos)
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

/// Borrow the columns of a finite column-major double matrix; R columns are
/// always contiguous.
pub(crate) fn matrix_columns<'a>(matrix: &Robj, name: &str) -> Result<(usize, Vec<&'a [f64]>)> {
    let data = real_slice(matrix, name)?;
    let (nrow, ncol) =
        matrix_dims(matrix).ok_or_else(|| err(format!("`{name}` must be a matrix")))?;
    let columns: Vec<&'a [f64]> = (0..ncol).map(|j| &data[j * nrow..(j + 1) * nrow]).collect();
    for (j, column) in columns.iter().enumerate() {
        check_finite(column, &format!("{name}[, {}]", j + 1))?;
    }
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

impl EffectColumns<'_> {
    /// Validate through the native constructor, so the rules live in one place.
    pub(crate) fn as_effect(&self) -> Result<Effect<'_>> {
        Effect::new(&self.levels, self.intercept, self.slopes.iter().copied()).map_err(build_err)
    }
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
            .enumerate()
            .map(|(term, c)| {
                c.as_effect()
                    .map_err(|e| err(format!("design term {}: {e}", term + 1)))
            })
            .collect()
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
                    .map_err(build_err)
            }
            Self::Effects(columns) => {
                let effects = Self::effects(&columns)?;
                crate::threads::run(|| Design::new(effects).map(Design::into_owned))?
                    .map_err(build_err)
            }
        }
    }
}

fn parse_effect<'a>(obj: &Robj) -> Result<EffectColumns<'a>> {
    if !obj.inherits("within_effect") {
        return Err(err("a design term must be created by Effect()"));
    }
    let levels = integer_slice(&get_field(obj, "levels"), "levels")?;
    let slopes = get_field(obj, "slopes");
    let slopes = slopes
        .as_list()
        .ok_or_else(|| err("`slopes` must be a list"))?
        .values()
        .enumerate()
        .map(|(slope, value)| real_slice(&value, &format!("slope {}", slope + 1)))
        .collect::<Result<Vec<_>>>()?;
    Ok(EffectColumns {
        levels: labels_to_u32(levels, "levels", None)?,
        intercept: scalar_bool(&get_field(obj, "intercept"), "intercept")?,
        slopes,
    })
}

/// Parse a `Design()` handle, an integer category matrix, or a list of `Effect()`s.
pub(crate) fn parse_design<'a>(design: &Robj) -> Result<DesignSource<'a>> {
    if design.inherits("within_design") {
        let ptr = handle::<Design<'static>>(&get_field(design, "ptr"), "Design")?;
        return Ok(DesignSource::Persistent(ptr.try_addr()?.clone()));
    }
    if design.is_matrix() {
        let data = integer_slice(design, "design")?;
        let (nrow, ncol) =
            matrix_dims(design).ok_or_else(|| err("a design matrix must be two-dimensional"))?;
        return Ok(DesignSource::Categories(Categories {
            labels: labels_to_u32(data, "design", Some(nrow))?,
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
        .map(|(term, obj)| {
            parse_effect(&obj).map_err(|e| err(format!("design term {}: {e}", term + 1)))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(DesignSource::Effects(columns))
}

// Validate one `Effect()` through the native constructor (called by `Effect()`).
#[extendr]
fn effect_validate_impl(effect: Robj) -> bool {
    guarded(|| parse_effect(&effect)?.as_effect().map(|_| true))
}

extendr_module! {
    mod convert;
    fn effect_validate_impl;
}
