//! R config objects ↔ native [`within::config`] conversions, plus the built
//! preconditioner handle.
//!
//! Mirrors `crates/within-py/src/config.rs`. Each R constructor
//! (`LsmrOptions()`, `LocalSolverConfig()`, `PreconditionerConfig$Adaptive()`,
//! ...) round-trips its arguments through [`config_normalize_impl`]: the
//! arguments are parsed into the native type — so validation and defaults live
//! only in the core — and converted back into a fully resolved, classed R
//! list. `Preconditioner$config` uses the same conversion, so two equal
//! configurations are `identical()` in R just as they compare equal in Python.

use extendr_api::{error::Result, prelude::*};

use within::config::{
    ApproxCholConfig, ApproxSchurConfig, LocalSolverConfig, LsmrOptions, PreconditionerConfig,
    ScalingConfig, ScalingFailure, SchurMode, Staleness,
};
use within::{Preconditioner, PreconditionerInput};

use crate::convert::{
    err, get_field, optional, or_throw, real_slice, scalar_f64, scalar_u32, scalar_u64,
    scalar_usize, usize_to_i32,
};

const LSMR_OPTIONS: &str = "within_lsmr_options";
const APPROX_CHOL: &str = "within_approx_chol_config";
const APPROX_SCHUR: &str = "within_approx_schur_config";
const SCHUR: &str = "within_schur";
const SCALING: &str = "within_scaling_config";
const LOCAL_SOLVER: &str = "within_local_solver_config";
const STALENESS: &str = "within_staleness";
const PRECONDITIONER_CONFIG: &str = "within_preconditioner_config";
const PRECONDITIONER: &str = "within_preconditioner";

/// Every config object also carries `within_config`, which hosts its `print` method.
fn classed(list: List, class: &str) -> Result<Robj> {
    let mut robj: Robj = list.into();
    robj.set_class([class, "within_config"])?;
    Ok(robj)
}

fn expect_class(obj: &Robj, class: &str, what: &str) -> Result<()> {
    if obj.inherits(class) {
        Ok(())
    } else {
        Err(err(format!(
            "`{what}` must be created by {}",
            constructor(class)
        )))
    }
}

fn constructor(class: &str) -> &'static str {
    match class {
        LSMR_OPTIONS => "LsmrOptions()",
        APPROX_CHOL => "ApproxCholConfig()",
        APPROX_SCHUR => "ApproxSchurConfig()",
        SCHUR => "Schur$approximate() or Schur$exact()",
        SCALING => "ScalingConfig()",
        LOCAL_SOLVER => "LocalSolverConfig()",
        STALENESS => "Staleness()",
        _ => "PreconditionerConfig$Off(), $Diagonal(), $Additive(), or $Adaptive()",
    }
}

fn int(value: usize, name: &str) -> Result<Robj> {
    usize_to_i32(value, name).map(Robj::from)
}

fn opt_robj(value: Option<Robj>) -> Robj {
    value.unwrap_or_else(nil_value)
}

// ---------------------------------------------------------------------------
// LSMR options
// ---------------------------------------------------------------------------

fn parse_lsmr(obj: &Robj) -> Result<LsmrOptions> {
    expect_class(obj, LSMR_OPTIONS, "options")?;
    let tol = scalar_f64(&get_field(obj, "tol"), "tol")?;
    if !tol.is_finite() || tol <= 0.0 {
        return Err(err("`tol` must be a positive finite number"));
    }
    Ok(LsmrOptions {
        tol,
        maxiter: scalar_usize(&get_field(obj, "maxiter"), "maxiter", 1)?,
        local_size: optional(&get_field(obj, "local_size"), |v| {
            scalar_usize(v, "local_size", 1)
        })?,
    })
}

fn lsmr_to_robj(options: &LsmrOptions) -> Result<Robj> {
    let local_size = options
        .local_size
        .map(|v| int(v, "local_size"))
        .transpose()?;
    classed(
        list!(
            tol = options.tol,
            maxiter = int(options.maxiter, "maxiter")?,
            local_size = opt_robj(local_size)
        ),
        LSMR_OPTIONS,
    )
}

/// `NULL` selects the library defaults.
pub(crate) fn resolve_lsmr(options: &Robj) -> Result<LsmrOptions> {
    if options.is_null() {
        Ok(LsmrOptions::default())
    } else {
        parse_lsmr(options)
    }
}

// ---------------------------------------------------------------------------
// Local solver configuration
// ---------------------------------------------------------------------------

fn parse_approx_chol(obj: &Robj) -> Result<ApproxCholConfig> {
    expect_class(obj, APPROX_CHOL, "approx_chol")?;
    Ok(ApproxCholConfig {
        seed: scalar_u64(&get_field(obj, "seed"), "seed")?,
        split_merge: optional(&get_field(obj, "split_merge"), |v| {
            scalar_u32(v, "split_merge", 1)
        })?,
    })
}

fn approx_chol_to_robj(config: &ApproxCholConfig) -> Result<Robj> {
    let split_merge = config
        .split_merge
        .map(|v| int(v as usize, "split_merge"))
        .transpose()?;
    classed(
        list!(
            seed = config.seed as f64,
            split_merge = opt_robj(split_merge)
        ),
        APPROX_CHOL,
    )
}

fn parse_approx_schur(obj: &Robj) -> Result<ApproxSchurConfig> {
    expect_class(obj, APPROX_SCHUR, "config")?;
    Ok(ApproxSchurConfig {
        split: scalar_u32(&get_field(obj, "split"), "split", 1)?,
    })
}

fn approx_schur_to_robj(config: &ApproxSchurConfig) -> Result<Robj> {
    classed(
        list!(split = int(config.split as usize, "split")?),
        APPROX_SCHUR,
    )
}

fn parse_schur(obj: &Robj) -> Result<SchurMode> {
    expect_class(obj, SCHUR, "schur")?;
    match get_field(obj, "mode").as_str() {
        Some("exact") => Ok(SchurMode::Exact),
        Some("approximate") => Ok(SchurMode::Approximate(
            optional(&get_field(obj, "config"), parse_approx_schur)?.unwrap_or_default(),
        )),
        _ => Err(err("`schur$mode` must be \"approximate\" or \"exact\"")),
    }
}

fn schur_to_robj(mode: &SchurMode) -> Result<Robj> {
    match mode {
        SchurMode::Exact => classed(list!(mode = "exact", config = nil_value()), SCHUR),
        SchurMode::Approximate(config) => classed(
            list!(mode = "approximate", config = approx_schur_to_robj(config)?),
            SCHUR,
        ),
    }
}

fn parse_scaling(obj: &Robj) -> Result<ScalingConfig> {
    expect_class(obj, SCALING, "scaling")?;
    let defaults = ScalingConfig::default();
    let on_failure = match optional(&get_field(obj, "on_failure"), |v| {
        v.as_str()
            .map(str::to_owned)
            .ok_or_else(|| err("`on_failure` must be \"warn\" or \"error\""))
    })?
    .as_deref()
    {
        None => defaults.on_failure,
        Some("warn") => ScalingFailure::Warn,
        Some("error") => ScalingFailure::Error,
        Some(other) => {
            return Err(err(format!(
                "`on_failure` must be \"warn\" or \"error\", got \"{other}\""
            )))
        }
    };
    Ok(ScalingConfig {
        tolerance: optional(&get_field(obj, "tolerance"), |v| scalar_f64(v, "tolerance"))?
            .unwrap_or(defaults.tolerance),
        max_iterations: optional(&get_field(obj, "max_iterations"), |v| {
            scalar_usize(v, "max_iterations", 0)
        })?
        .unwrap_or(defaults.max_iterations),
        on_failure,
    })
}

fn scaling_to_robj(config: &ScalingConfig) -> Result<Robj> {
    let on_failure = match config.on_failure {
        ScalingFailure::Warn => "warn",
        ScalingFailure::Error => "error",
    };
    classed(
        list!(
            tolerance = config.tolerance,
            max_iterations = int(config.max_iterations, "max_iterations")?,
            on_failure = on_failure
        ),
        SCALING,
    )
}

fn parse_local_solver(obj: &Robj) -> Result<LocalSolverConfig> {
    expect_class(obj, LOCAL_SOLVER, "local_solver")?;
    let defaults = LocalSolverConfig::default();
    Ok(LocalSolverConfig {
        approx_chol: optional(&get_field(obj, "approx_chol"), parse_approx_chol)?
            .unwrap_or(defaults.approx_chol),
        schur: optional(&get_field(obj, "schur"), parse_schur)?.unwrap_or(defaults.schur),
        dense_threshold: optional(&get_field(obj, "dense_threshold"), |v| {
            scalar_usize(v, "dense_threshold", 0)
        })?
        .unwrap_or(defaults.dense_threshold),
        scaling: optional(&get_field(obj, "scaling"), parse_scaling)?.unwrap_or(defaults.scaling),
    })
}

fn local_solver_to_robj(config: &LocalSolverConfig) -> Result<Robj> {
    classed(
        list!(
            approx_chol = approx_chol_to_robj(&config.approx_chol)?,
            schur = schur_to_robj(&config.schur)?,
            dense_threshold = int(config.dense_threshold, "dense_threshold")?,
            scaling = scaling_to_robj(&config.scaling)?
        ),
        LOCAL_SOLVER,
    )
}

fn parse_staleness(obj: &Robj) -> Result<Staleness> {
    expect_class(obj, STALENESS, "stall")?;
    let defaults = Staleness::default();
    let window = optional(&get_field(obj, "window"), |v| scalar_usize(v, "window", 0))?
        .unwrap_or(defaults.window());
    let threshold = optional(&get_field(obj, "threshold"), |v| scalar_f64(v, "threshold"))?
        .unwrap_or(defaults.threshold());
    Staleness::try_new(window, threshold).map_err(|e| err(e.to_string()))
}

fn staleness_to_robj(stall: &Staleness) -> Result<Robj> {
    classed(
        list!(
            window = int(stall.window(), "window")?,
            threshold = stall.threshold()
        ),
        STALENESS,
    )
}

// ---------------------------------------------------------------------------
// Preconditioner configuration
// ---------------------------------------------------------------------------

fn parse_preconditioner_config(obj: &Robj) -> Result<PreconditionerConfig> {
    expect_class(obj, PRECONDITIONER_CONFIG, "preconditioner")?;
    let local_solver = || -> Result<LocalSolverConfig> {
        Ok(optional(&get_field(obj, "local_solver"), parse_local_solver)?.unwrap_or_default())
    };
    match get_field(obj, "variant").as_str() {
        Some("Off") => Ok(PreconditionerConfig::Off),
        Some("Diagonal") => Ok(PreconditionerConfig::Diagonal),
        Some("Additive") => Ok(PreconditionerConfig::Additive {
            local_solver: local_solver()?,
        }),
        Some("Adaptive") => Ok(PreconditionerConfig::Adaptive {
            local_solver: local_solver()?,
            stall: optional(&get_field(obj, "stall"), parse_staleness)?.unwrap_or_default(),
        }),
        _ => Err(err(
            "`preconditioner$variant` must be \"Off\", \"Diagonal\", \"Additive\", or \"Adaptive\"",
        )),
    }
}

fn preconditioner_config_to_robj(config: &PreconditionerConfig) -> Result<Robj> {
    let list = match config {
        PreconditionerConfig::Off => list!(variant = "Off"),
        PreconditionerConfig::Diagonal => list!(variant = "Diagonal"),
        PreconditionerConfig::Additive { local_solver } => list!(
            variant = "Additive",
            local_solver = local_solver_to_robj(local_solver)?
        ),
        PreconditionerConfig::Adaptive {
            local_solver,
            stall,
        } => list!(
            variant = "Adaptive",
            local_solver = local_solver_to_robj(local_solver)?,
            stall = staleness_to_robj(stall)?
        ),
        _ => return Err(err("unsupported preconditioner configuration")),
    };
    classed(list, PRECONDITIONER_CONFIG)
}

/// Resolve the R `preconditioner` argument: `NULL` (library default), a
/// `PreconditionerConfig$...()` object, or a built `Preconditioner`.
pub(crate) fn resolve_preconditioner(obj: &Robj) -> Result<PreconditionerInput> {
    if obj.is_null() {
        return Ok(PreconditionerInput::Default);
    }
    if obj.inherits(PRECONDITIONER) {
        return Ok(PreconditionerInput::Prebuilt(preconditioner_handle(
            &get_field(obj, "ptr"),
        )?));
    }
    if obj.inherits(PRECONDITIONER_CONFIG) {
        return parse_preconditioner_config(obj).map(PreconditionerInput::Config);
    }
    Err(err(
        "`preconditioner` must be NULL, PreconditionerConfig$Off(), $Diagonal(), \
         $Additive(...), $Adaptive(...), or a built Preconditioner",
    ))
}

// Validate a config object and return it fully resolved, with library defaults filled in.
#[extendr]
fn config_normalize_impl(obj: Robj) -> Robj {
    or_throw((|| -> Result<Robj> {
        let class = obj
            .class()
            .and_then(|mut classes| classes.next().map(str::to_owned))
            .unwrap_or_default();
        match class.as_str() {
            LSMR_OPTIONS => lsmr_to_robj(&parse_lsmr(&obj)?),
            APPROX_CHOL => approx_chol_to_robj(&parse_approx_chol(&obj)?),
            APPROX_SCHUR => approx_schur_to_robj(&parse_approx_schur(&obj)?),
            SCHUR => schur_to_robj(&parse_schur(&obj)?),
            SCALING => scaling_to_robj(&parse_scaling(&obj)?),
            LOCAL_SOLVER => local_solver_to_robj(&parse_local_solver(&obj)?),
            STALENESS => staleness_to_robj(&parse_staleness(&obj)?),
            PRECONDITIONER_CONFIG => {
                preconditioner_config_to_robj(&parse_preconditioner_config(&obj)?)
            }
            other => Err(err(format!("unknown configuration class \"{other}\""))),
        }
    })())
}

// ---------------------------------------------------------------------------
// Built preconditioner handle
// ---------------------------------------------------------------------------

/// Clone the native preconditioner behind an R handle; the clone is O(1).
fn preconditioner_handle(ptr: &Robj) -> Result<Preconditioner> {
    let ptr = ExternalPtr::<Preconditioner>::try_from(ptr.clone())
        .map_err(|_| err("invalid Preconditioner handle"))?;
    Ok(ptr.try_addr()?.clone())
}

// Apply a built preconditioner: `y = M^{-1} x`.
#[extendr]
fn preconditioner_apply_impl(preconditioner: ExternalPtr<Preconditioner>, x: Robj) -> Vec<f64> {
    or_throw((|| -> Result<Vec<f64>> {
        let inner = preconditioner.try_addr()?;
        let x = real_slice(&x, "x")?;
        if x.len() != inner.ncols() {
            return Err(err(format!(
                "`x` has length {} but the preconditioner expects {}",
                x.len(),
                inner.ncols()
            )));
        }
        let mut y = vec![0.0; inner.nrows()];
        crate::threads::run(|| inner.apply(x, &mut y))?.map_err(|e| err(e.to_string()))?;
        Ok(y)
    })())
}

// Shape, variant, configuration, and build time of a built preconditioner.
#[extendr]
fn preconditioner_info_impl(preconditioner: ExternalPtr<Preconditioner>) -> List {
    or_throw((|| -> Result<List> {
        let inner = preconditioner.try_addr()?;
        Ok(list!(
            nrows = int(inner.nrows(), "nrows")?,
            ncols = int(inner.ncols(), "ncols")?,
            variant = inner.variant_name(),
            config = preconditioner_config_to_robj(&inner.config())?,
            build_duration_seconds = inner.build_duration().as_secs_f64()
        ))
    })())
}

// Serialize a built preconditioner to raw bytes (the R analogue of pickling).
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
        let preconditioner: Preconditioner = postcard::from_bytes(data.as_slice())
            .map_err(|e| err(format!("failed to deserialize preconditioner: {e}")))?;
        Ok(ExternalPtr::new(preconditioner))
    })())
}

extendr_module! {
    mod config;
    fn config_normalize_impl;
    fn preconditioner_apply_impl;
    fn preconditioner_info_impl;
    fn preconditioner_serialize_impl;
    fn preconditioner_deserialize_impl;
}
