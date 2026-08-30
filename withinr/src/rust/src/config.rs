//! R config objects → native [`within::config`] conversions, plus the
//! built-preconditioner handle exposed to R.
//!
//! Mirrors `crates/within-py/src/config.rs`: this is where R-facing
//! configuration values (classed lists, shortcut strings) are parsed into
//! native `within` config types.

use extendr_api::{error::Result, prelude::*};

use within::config::{
    ApproxCholConfig, ApproxSchurConfig, LocalSolverConfig, LsmrOptions, PreconditionerConfig,
    ReductionStrategy, ScalingConfig, ScalingFailure, SchurMode,
};
use within::Preconditioner;

use crate::convert::{
    err, get_field_or_null, or_throw, parse_nonnegative_f64, parse_nonnegative_integer,
    parse_optional_u32, parse_positive_f64, parse_positive_u32, usize_to_i32,
};

/// Native interpretation of the R `preconditioner` argument.
///
/// A pre-built [`Preconditioner`] takes the reuse path; everything else is a
/// [`PreconditionerConfig`] (or `None` for the library default) to build from.
/// Mirrors `PrecondInput` on the Python side.
pub(crate) enum PreconditionerArg {
    Config(Option<PreconditionerConfig>),
    Built(Preconditioner),
}

// ---------------------------------------------------------------------------
// LSMR options
// ---------------------------------------------------------------------------

/// Resolve the R `options` argument into native [`LsmrOptions`].
///
/// `NULL` means library defaults; anything else must be created by
/// `LsmrOptions(...)`. Mirrors `resolve_lsmr_config` on the Python side;
/// value validation lives here, not in the R constructor.
pub(crate) fn parse_lsmr_options(options: &Robj) -> Result<LsmrOptions> {
    if options.is_null() {
        return Ok(LsmrOptions::default());
    }
    if !options.inherits("within_lsmr_options") {
        return Err(err("options must be created by LsmrOptions(...) or NULL"));
    }

    let tol = parse_positive_f64(&get_field_or_null(options, "tol"), "tol")?;

    let maxiter = parse_nonnegative_integer(&get_field_or_null(options, "maxiter"), "maxiter")?;
    if maxiter < 1 {
        return Err(err("maxiter must be >= 1"));
    }

    let local_size_obj = get_field_or_null(options, "local_size");
    let local_size = if local_size_obj.is_null() {
        None
    } else {
        let value = parse_nonnegative_integer(&local_size_obj, "local_size")?;
        if value < 1 {
            return Err(err("local_size must be NULL or >= 1"));
        }
        Some(value as usize)
    };

    Ok(LsmrOptions {
        tol,
        maxiter: maxiter as usize,
        local_size,
    })
}

// ---------------------------------------------------------------------------
// Preconditioner config conversion
// ---------------------------------------------------------------------------

fn parse_approx_chol_config(obj: &Robj) -> Result<ApproxCholConfig> {
    if obj.is_null() {
        return Ok(LocalSolverConfig::default().approx_chol);
    }
    if !obj.inherits("within_approx_chol_config") {
        return Err(err(
            "approx_chol must be an object created by ApproxCholConfig(...) or NULL",
        ));
    }

    let seed_obj = get_field_or_null(obj, "seed");
    let split_merge_obj = get_field_or_null(obj, "split_merge");

    let seed = parse_nonnegative_integer(&seed_obj, "approx_chol$seed")?;

    Ok(ApproxCholConfig {
        seed: seed as u64,
        split_merge: parse_optional_u32(&split_merge_obj, "approx_chol$split_merge")?,
    })
}

fn parse_approx_schur_config(obj: &Robj) -> Result<ApproxSchurConfig> {
    if obj.is_null() {
        return Ok(ApproxSchurConfig::default());
    }
    if !obj.inherits("within_approx_schur_config") {
        return Err(err(
            "schur$config must be an object created by ApproxSchurConfig(...) or NULL",
        ));
    }

    let seed_obj = get_field_or_null(obj, "seed");
    let split_obj = get_field_or_null(obj, "split");

    let seed = parse_nonnegative_integer(&seed_obj, "schur$config$seed")?;

    Ok(ApproxSchurConfig {
        seed: seed as u64,
        split: parse_positive_u32(&split_obj, "schur$config$split")?,
    })
}

fn parse_schur_mode(obj: &Robj) -> Result<SchurMode> {
    if obj.is_null() {
        return Ok(SchurMode::default());
    }
    if obj.inherits("within_schur_exact") {
        return Ok(SchurMode::Exact);
    }
    if obj.inherits("within_schur_approximate") {
        return parse_approx_schur_config(&get_field_or_null(obj, "config"))
            .map(SchurMode::Approximate);
    }
    Err(err(
        "schur must be Schur$Exact, Schur$Approximate(...), or NULL",
    ))
}

fn parse_scaling_config(obj: &Robj) -> Result<ScalingConfig> {
    if obj.is_null() {
        return Ok(ScalingConfig::default());
    }
    if !obj.inherits("within_scaling_config") {
        return Err(err(
            "scaling must be an object created by ScalingConfig(...) or NULL",
        ));
    }

    let defaults = ScalingConfig::default();
    let tolerance_obj = get_field_or_null(obj, "tolerance");
    let max_iterations_obj = get_field_or_null(obj, "max_iterations");
    let on_failure_obj = get_field_or_null(obj, "on_failure");

    let tolerance = if tolerance_obj.is_null() {
        defaults.tolerance
    } else {
        parse_positive_f64(&tolerance_obj, "scaling$tolerance")?
    };
    let max_iterations = if max_iterations_obj.is_null() {
        defaults.max_iterations
    } else {
        let value = parse_nonnegative_integer(&max_iterations_obj, "scaling$max_iterations")?;
        if value < 1 {
            return Err(err("scaling$max_iterations must be >= 1"));
        }
        usize::try_from(value).map_err(|_| err("scaling$max_iterations is too large"))?
    };
    let on_failure = if on_failure_obj.is_null() {
        defaults.on_failure
    } else {
        match on_failure_obj.as_str() {
            Some("warn") => ScalingFailure::Warn,
            Some("error") => ScalingFailure::Error,
            _ => {
                return Err(err(
                    "scaling$on_failure must be ScalingFailure$Warn or ScalingFailure$Error",
                ))
            }
        }
    };

    Ok(ScalingConfig {
        tolerance,
        max_iterations,
        on_failure,
    })
}

fn parse_local_solver_config(obj: &Robj) -> Result<LocalSolverConfig> {
    if obj.is_null() {
        return Ok(LocalSolverConfig::default());
    }
    if !obj.inherits("within_local_solver_config") {
        return Err(err(
            "local_solver must be an object created by LocalSolverConfig(...) or NULL",
        ));
    }

    let approx_chol_obj = get_field_or_null(obj, "approx_chol");
    let schur_obj = get_field_or_null(obj, "schur");
    let dense_threshold_obj = get_field_or_null(obj, "dense_threshold");
    let scaling_obj = get_field_or_null(obj, "scaling");
    let ridge_obj = get_field_or_null(obj, "ridge");

    let mut config = LocalSolverConfig::default();
    if !approx_chol_obj.is_null() {
        config.approx_chol = parse_approx_chol_config(&approx_chol_obj)?;
    }
    if !schur_obj.is_null() {
        config.schur = parse_schur_mode(&schur_obj)?;
    }
    if !dense_threshold_obj.is_null() {
        let value =
            parse_nonnegative_integer(&dense_threshold_obj, "local_solver$dense_threshold")?;
        config.dense_threshold =
            usize::try_from(value).map_err(|_| err("local_solver$dense_threshold is too large"))?;
    }
    if !scaling_obj.is_null() {
        config.scaling = parse_scaling_config(&scaling_obj)?;
    }
    if !ridge_obj.is_null() {
        config.ridge = parse_nonnegative_f64(&ridge_obj, "local_solver$ridge")?;
    }
    Ok(config)
}

fn parse_reduction_strategy(obj: &Robj) -> Result<ReductionStrategy> {
    if obj.is_null() {
        return Ok(ReductionStrategy::default());
    }
    let Some(name) = obj.as_str() else {
        return Err(err("reduction must be a character scalar"));
    };
    match name {
        "auto" => Ok(ReductionStrategy::Auto),
        "atomic_scatter" => Ok(ReductionStrategy::AtomicScatter),
        "parallel_reduction" => Ok(ReductionStrategy::ParallelReduction),
        other => Err(err(format!(
            "unknown reduction strategy '{other}'; use 'auto', 'atomic_scatter', or 'parallel_reduction'"
        ))),
    }
}

fn parse_preconditioner_string(name: &str) -> Result<PreconditionerArg> {
    match name {
        "additive" => Ok(PreconditionerArg::Config(Some(
            PreconditionerConfig::default(),
        ))),
        "off" => Ok(PreconditionerArg::Config(Some(PreconditionerConfig::Off))),
        "diagonal" => Ok(PreconditionerArg::Config(Some(
            PreconditionerConfig::Diagonal,
        ))),
        other => Err(err(format!(
            "unknown preconditioner '{other}'; use PreconditionerConfig$Additive, \
             PreconditionerConfig$Off, PreconditionerConfig$Diagonal, AdditiveSchwarz(...), \
             a Preconditioner, or NULL"
        ))),
    }
}

pub(crate) fn parse_preconditioner(preconditioner: Robj) -> Result<PreconditionerArg> {
    if preconditioner.is_null() {
        return Ok(PreconditionerArg::Config(None));
    }

    // A `within_preconditioner` R object is an environment holding the
    // external pointer in `$ptr`; a bare external pointer is accepted too.
    let ptr_obj = if preconditioner.inherits("within_preconditioner") {
        get_field_or_null(&preconditioner, "ptr")
    } else {
        preconditioner.clone()
    };
    if let Ok(ptr) = ExternalPtr::<Preconditioner>::try_from(ptr_obj) {
        return Ok(PreconditionerArg::Built(ptr.try_addr()?.clone()));
    }

    if let Some(name) = preconditioner.as_str() {
        return parse_preconditioner_string(name);
    }

    if preconditioner.inherits("within_additive_schwarz") {
        let local_solver =
            parse_local_solver_config(&get_field_or_null(&preconditioner, "local_solver"))?;
        let reduction = parse_reduction_strategy(&get_field_or_null(&preconditioner, "reduction"))?;
        return Ok(PreconditionerArg::Config(Some(
            PreconditionerConfig::Additive {
                local_solver,
                reduction,
            },
        )));
    }

    Err(err(
        "preconditioner must be NULL, 'additive', 'off', 'diagonal', AdditiveSchwarz(...), or a Preconditioner",
    ))
}

// ---------------------------------------------------------------------------
// Preconditioner handle API
// ---------------------------------------------------------------------------

// Apply a built preconditioner: y = M^{-1} x.
#[extendr]
fn preconditioner_apply_impl(preconditioner: ExternalPtr<Preconditioner>, x: &[f64]) -> Vec<f64> {
    or_throw((|| -> Result<Vec<f64>> {
        let inner = preconditioner.try_addr()?;
        if x.len() != inner.ncols() {
            return Err(err(format!(
                "x has length {} but preconditioner expects {}",
                x.len(),
                inner.ncols()
            )));
        }
        let mut y = vec![0.0; inner.nrows()];
        inner.apply(x, &mut y).map_err(|e| err(e.to_string()))?;
        Ok(y)
    })())
}

// Number of rows in a built preconditioner.
#[extendr]
fn preconditioner_nrows_impl(preconditioner: ExternalPtr<Preconditioner>) -> i32 {
    or_throw((|| -> Result<i32> {
        usize_to_i32(preconditioner.try_addr()?.nrows(), "nrows")
    })())
}

// Number of columns in a built preconditioner.
#[extendr]
fn preconditioner_ncols_impl(preconditioner: ExternalPtr<Preconditioner>) -> i32 {
    or_throw((|| -> Result<i32> {
        usize_to_i32(preconditioner.try_addr()?.ncols(), "ncols")
    })())
}

// Concrete preconditioner variant name (used by the R print method, matching
// the Python __repr__).
#[extendr]
fn preconditioner_variant_impl(preconditioner: ExternalPtr<Preconditioner>) -> String {
    or_throw((|| -> Result<String> {
        Ok(preconditioner.try_addr()?.variant_name().to_string())
    })())
}

fn reduction_name(reduction: ReductionStrategy) -> &'static str {
    match reduction {
        ReductionStrategy::Auto => "auto",
        ReductionStrategy::AtomicScatter => "atomic_scatter",
        ReductionStrategy::ParallelReduction => "parallel_reduction",
    }
}

fn scaling_failure_name(on_failure: ScalingFailure) -> &'static str {
    match on_failure {
        ScalingFailure::Warn => "warn",
        ScalingFailure::Error => "error",
    }
}

fn local_solver_to_list(config: &LocalSolverConfig) -> List {
    let schur = match &config.schur {
        SchurMode::Approximate(approx) => list!(
            variant = "approximate",
            seed = approx.seed as f64,
            split = approx.split as f64
        ),
        SchurMode::Exact => list!(variant = "exact"),
    };
    list!(
        approx_chol = list!(
            seed = config.approx_chol.seed as f64,
            split_merge = config.approx_chol.split_merge.map(f64::from)
        ),
        schur = schur,
        dense_threshold = config.dense_threshold as f64,
        scaling = list!(
            tolerance = config.scaling.tolerance,
            max_iterations = config.scaling.max_iterations as f64,
            on_failure = scaling_failure_name(config.scaling.on_failure)
        ),
        ridge = config.ridge
    )
}

// Complete normalized construction configuration for a built preconditioner.
// The wildcard keeps this bridge source-compatible when the non-exhaustive
// native enum gains a strategy: callers still receive the concrete variant.
#[extendr]
fn preconditioner_config_impl(preconditioner: ExternalPtr<Preconditioner>) -> List {
    or_throw((|| -> Result<List> {
        let inner = preconditioner.try_addr()?;
        let variant = inner.variant_name().to_string();
        Ok(match inner.config() {
            PreconditionerConfig::Off => list!(variant = variant),
            PreconditionerConfig::Diagonal => list!(variant = variant),
            PreconditionerConfig::Additive {
                local_solver,
                reduction,
            } => list!(
                variant = variant,
                local_solver = local_solver_to_list(local_solver),
                reduction = reduction_name(*reduction)
            ),
            _ => list!(variant = variant),
        })
    })())
}

// Original preconditioner build duration, preserved across serialization.
#[extendr]
fn preconditioner_build_duration_impl(preconditioner: ExternalPtr<Preconditioner>) -> f64 {
    or_throw((|| -> Result<f64> {
        Ok(preconditioner.try_addr()?.build_duration().as_secs_f64())
    })())
}

// Serialize a built preconditioner into raw bytes.
#[extendr]
fn preconditioner_serialize_impl(preconditioner: ExternalPtr<Preconditioner>) -> Raw {
    or_throw((|| -> Result<Raw> {
        let bytes =
            postcard::to_stdvec(preconditioner.try_addr()?).map_err(|e| err(e.to_string()))?;
        Ok(Raw::from_bytes(&bytes))
    })())
}

// Deserialize a built preconditioner from raw bytes.
#[extendr]
fn preconditioner_deserialize_impl(data: Raw) -> ExternalPtr<Preconditioner> {
    or_throw((|| -> Result<ExternalPtr<Preconditioner>> {
        let preconditioner: Preconditioner =
            postcard::from_bytes(data.as_slice()).map_err(|e| err(e.to_string()))?;
        Ok(ExternalPtr::new(preconditioner))
    })())
}

extendr_module! {
    mod config;
    fn preconditioner_apply_impl;
    fn preconditioner_nrows_impl;
    fn preconditioner_ncols_impl;
    fn preconditioner_variant_impl;
    fn preconditioner_config_impl;
    fn preconditioner_build_duration_impl;
    fn preconditioner_serialize_impl;
    fn preconditioner_deserialize_impl;
}
