//! Solve entry points and the persistent `Design` / `Solver` handles.
//!
//! Mirrors `crates/within-py/src/api.rs`: one-shot `solve` / `solve_batch`
//! borrow the R inputs for the duration of the call, while `Design()` and
//! `Solver()` own native copies so they stay valid however the R inputs change.

use std::sync::atomic::{AtomicUsize, Ordering};

use extendr_api::{error::Result, prelude::*};

use within::{Design, Solver};

use crate::config::{resolve_lsmr, resolve_preconditioner};
use crate::convert::{
    err, matrix_columns, optional_weights, or_throw, parse_design, real_slice, usize_to_i32,
    DesignSource,
};
use crate::results::{batch_result_to_list, result_to_list, warning_strings};

/// Run `$body` with `$design` bound to the native form of a [`DesignSource`].
macro_rules! with_design {
    ($source:expr, |$design:ident| $body:expr) => {
        match $source {
            DesignSource::Persistent($design) => $body,
            DesignSource::Categories(categories) => {
                let $design = categories.view()?;
                $body
            }
            DesignSource::Effects(columns) => {
                let $design = DesignSource::effects(&columns)?;
                $body
            }
        }
    };
}

fn check_rows(y_rows: usize, n_obs: Option<usize>) -> Result<()> {
    match n_obs {
        Some(n_obs) if n_obs != y_rows => Err(err(format!(
            "`Y` has {y_rows} rows but the design has {n_obs} observations"
        ))),
        _ => Ok(()),
    }
}

// ---------------------------------------------------------------------------
// One-shot solves
// ---------------------------------------------------------------------------

#[extendr]
fn solve_impl(design: Robj, y: Robj, weights: Robj, options: Robj, preconditioner: Robj) -> List {
    or_throw((|| -> Result<List> {
        let lsmr = resolve_lsmr(&options)?;
        let precond = resolve_preconditioner(&preconditioner)?;
        let y = real_slice(&y, "y")?;
        let weights = optional_weights(&weights)?;
        let source = parse_design(&design)?;
        let result = with_design!(source, |design| crate::threads::run(|| {
            within::solve(design, y, weights, &lsmr, precond)
        })?)
        .map_err(|e| err(e.to_string()))?;
        let warnings = warning_strings(&result.warnings);
        result_to_list(result, warnings)
    })())
}

#[extendr]
fn solve_batch_impl(
    design: Robj,
    y_matrix: Robj,
    weights: Robj,
    options: Robj,
    preconditioner: Robj,
) -> List {
    or_throw((|| -> Result<List> {
        let lsmr = resolve_lsmr(&options)?;
        let precond = resolve_preconditioner(&preconditioner)?;
        let (n_rows, columns) = matrix_columns(&y_matrix, "Y")?;
        let weights = optional_weights(&weights)?;
        let source = parse_design(&design)?;
        check_rows(n_rows, source.n_obs())?;
        let result = with_design!(source, |design| crate::threads::run(|| {
            within::solve_batch(design, &columns, weights, &lsmr, precond)
        })?)
        .map_err(|e| err(e.to_string()))?;
        let warnings = warning_strings(&result.warnings);
        batch_result_to_list(result, warnings)
    })())
}

// ---------------------------------------------------------------------------
// Persistent design
// ---------------------------------------------------------------------------

// Build an owned native design; returns `ptr`, `n_obs`, and `n_dofs`.
#[extendr]
fn design_new_impl(design: Robj) -> List {
    or_throw((|| -> Result<List> {
        let design = parse_design(&design)?.into_owned()?;
        Ok(list!(
            n_obs = usize_to_i32(design.n_obs(), "n_obs")?,
            n_dofs = usize_to_i32(design.n_dofs(), "n_dofs")?,
            ptr = ExternalPtr::<Design<'static>>::new(design)
        ))
    })())
}

// ---------------------------------------------------------------------------
// Persistent solver
// ---------------------------------------------------------------------------

/// A persistent solver plus how many of its build warnings R has already seen.
///
/// Under `Adaptive` the Schwarz rung is built mid-solve, so its warnings
/// surface from the solve that escalates — each one exactly once.
pub(crate) struct SolverHandle {
    solver: Solver<'static>,
    reported_warnings: AtomicUsize,
}

impl SolverHandle {
    fn take_new_warnings(&self) -> Vec<String> {
        let warnings = self.solver.warnings();
        let seen = self
            .reported_warnings
            .swap(warnings.len(), Ordering::Relaxed);
        warning_strings(warnings.get(seen..).unwrap_or_default())
    }
}

// Build a persistent solver; returns `ptr`, `n_obs`, `n_dofs`, and build `warnings`.
#[extendr]
fn solver_new_impl(design: Robj, weights: Robj, preconditioner: Robj) -> List {
    or_throw((|| -> Result<List> {
        let precond = resolve_preconditioner(&preconditioner)?;
        let weights = optional_weights(&weights)?;
        let design = parse_design(&design)?.into_owned()?;
        let solver = crate::threads::run(|| Solver::new(design, weights, precond))?
            .map_err(|e| err(e.to_string()))?;
        let handle = SolverHandle {
            solver,
            reported_warnings: AtomicUsize::new(0),
        };
        let warnings = handle.take_new_warnings();
        Ok(list!(
            n_obs = usize_to_i32(handle.solver.n_obs(), "n_obs")?,
            n_dofs = usize_to_i32(handle.solver.n_dofs(), "n_dofs")?,
            warnings = warnings,
            ptr = ExternalPtr::new(handle)
        ))
    })())
}

#[extendr]
fn solver_solve_impl(solver: ExternalPtr<SolverHandle>, y: Robj, options: Robj) -> List {
    or_throw((|| -> Result<List> {
        let lsmr = resolve_lsmr(&options)?;
        let handle = solver.try_addr()?;
        let y = real_slice(&y, "y")?;
        let result = crate::threads::run(|| handle.solver.solve(y, &lsmr))?
            .map_err(|e| err(e.to_string()))?;
        result_to_list(result, handle.take_new_warnings())
    })())
}

#[extendr]
fn solver_solve_batch_impl(
    solver: ExternalPtr<SolverHandle>,
    y_matrix: Robj,
    options: Robj,
) -> List {
    or_throw((|| -> Result<List> {
        let lsmr = resolve_lsmr(&options)?;
        let handle = solver.try_addr()?;
        let (n_rows, columns) = matrix_columns(&y_matrix, "Y")?;
        check_rows(n_rows, Some(handle.solver.n_obs()))?;
        let result = crate::threads::run(|| handle.solver.solve_batch(&columns, &lsmr))?
            .map_err(|e| err(e.to_string()))?;
        batch_result_to_list(result, handle.take_new_warnings())
    })())
}

// The preconditioner in use (`NULL` when unpreconditioned). Under `Adaptive`
// it is the Schwarz map once built, else the diagonal base carrying the ladder.
#[extendr]
fn solver_preconditioner_impl(solver: ExternalPtr<SolverHandle>) -> Robj {
    or_throw((|| -> Result<Robj> {
        Ok(solver
            .try_addr()?
            .solver
            .preconditioner()
            .map_or_else(nil_value, |p| ExternalPtr::new(p.clone()).into()))
    })())
}

// Whether an `Adaptive` solve has handed off to Schwarz; always `FALSE` otherwise.
#[extendr]
fn solver_has_escalated_impl(solver: ExternalPtr<SolverHandle>) -> bool {
    or_throw(
        solver
            .try_addr()
            .map(|handle| handle.solver.has_escalated()),
    )
}

extendr_module! {
    mod api;
    fn solve_impl;
    fn solve_batch_impl;
    fn design_new_impl;
    fn solver_new_impl;
    fn solver_solve_impl;
    fn solver_solve_batch_impl;
    fn solver_preconditioner_impl;
    fn solver_has_escalated_impl;
}
