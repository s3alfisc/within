//! Solve entry points and the persistent solver handle.
//!
//! Mirrors `crates/within-py/src/api.rs`: the one-shot `solve`/`solve_batch`
//! bindings plus the persistent `Solver` surface.

use extendr_api::{error::Result, prelude::*};

use within::{
    solve as solve_native, solve_batch as solve_batch_native, Preconditioner,
    Solver as NativeSolver,
};

use crate::config::{parse_lsmr_options, parse_preconditioner, PreconditionerArg};
use crate::convert::{err, or_throw, parse_design, usize_to_i32, weights_slice};
use crate::results::{batch_result_to_list, result_to_list};

/// Owned-store solver held behind the R external pointer.
pub(crate) type SolverHandle = NativeSolver<'static>;

/// Borrow the columns of an R numeric matrix as slices.
///
/// R matrices are column-major, so every RHS column is already contiguous —
/// no copy needed (the Python glue needs Cow columns for strided input; R
/// input never is).
fn matrix_columns(matrix: &RMatrix<f64>) -> Vec<&[f64]> {
    let data = matrix.data();
    let nrow = matrix.nrows();
    (0..matrix.ncols())
        .map(|j| &data[j * nrow..(j + 1) * nrow])
        .collect()
}

// Persistent designs share owned observation storage across weighted solvers.
#[extendr]
fn design_new_impl(design: Robj) -> ExternalPtr<within::Design<'static>> {
    or_throw(parse_design(&design).map(ExternalPtr::new))
}

#[extendr]
fn solver_has_escalated_impl(solver: ExternalPtr<SolverHandle>) -> bool {
    or_throw(solver.try_addr().map(|handle| handle.has_escalated()))
}

// ---------------------------------------------------------------------------
// One-shot solve API
// ---------------------------------------------------------------------------

// Solve fixed-effects normal equations for a single response vector.
#[extendr]
fn solve_impl(design: Robj, y: &[f64], options: Robj, weights: Robj, preconditioner: Robj) -> List {
    or_throw((|| -> Result<List> {
        let design = parse_design(&design)?;
        let lsmr = parse_lsmr_options(&options)?;
        let weights = weights_slice(&weights)?;

        match parse_preconditioner(preconditioner)? {
            PreconditionerArg::Config(config) => {
                crate::threads::run(|| solve_native(design, y, weights, &lsmr, config.as_ref()))?
                    .map_err(|e| err(e.to_string()))
                    .and_then(result_to_list)
            }
            PreconditionerArg::Built(built) => {
                crate::threads::run(|| solve_native(design, y, weights, &lsmr, built))?
                    .map_err(|e| err(e.to_string()))
                    .and_then(result_to_list)
            }
        }
    })())
}

// Solve fixed-effects normal equations for multiple response vectors.
#[extendr]
fn solve_batch_impl(
    design: Robj,
    y_matrix: RMatrix<f64>,
    options: Robj,
    weights: Robj,
    preconditioner: Robj,
) -> List {
    or_throw((|| -> Result<List> {
        let design = parse_design(&design)?;
        if y_matrix.nrows() != design.n_obs() {
            return Err(err(format!(
                "Y has {} rows but design has {} observations",
                y_matrix.nrows(),
                design.n_obs()
            )));
        }
        let lsmr = parse_lsmr_options(&options)?;
        let weights = weights_slice(&weights)?;
        let column_refs = matrix_columns(&y_matrix);

        match parse_preconditioner(preconditioner)? {
            PreconditionerArg::Config(config) => crate::threads::run(|| {
                solve_batch_native(design, &column_refs, weights, &lsmr, config.as_ref())
            })?
            .map_err(|e| err(e.to_string()))
            .and_then(batch_result_to_list),
            PreconditionerArg::Built(built) => crate::threads::run(|| {
                solve_batch_native(design, &column_refs, weights, &lsmr, built)
            })?
            .map_err(|e| err(e.to_string()))
            .and_then(batch_result_to_list),
        }
    })())
}

// ---------------------------------------------------------------------------
// Persistent solver API
// ---------------------------------------------------------------------------

// Build a persistent solver that can be reused across multiple solves.
#[extendr]
fn solver_new_impl(design: Robj, weights: Robj, preconditioner: Robj) -> ExternalPtr<SolverHandle> {
    or_throw((|| -> Result<ExternalPtr<SolverHandle>> {
        let weights = weights_slice(&weights)?;
        let design = parse_design(&design)?;

        let solver = match parse_preconditioner(preconditioner)? {
            PreconditionerArg::Config(config) => {
                crate::threads::run(|| NativeSolver::new(design, weights, config.as_ref()))?
            }
            PreconditionerArg::Built(built) => {
                crate::threads::run(|| NativeSolver::new(design, weights, built))?
            }
        }
        .map_err(|e| err(e.to_string()))?;

        Ok(ExternalPtr::new(solver))
    })())
}

// Solve one response vector with a persistent solver.
#[extendr]
fn solver_solve_impl(solver: ExternalPtr<SolverHandle>, y: &[f64], options: Robj) -> List {
    or_throw((|| -> Result<List> {
        let lsmr = parse_lsmr_options(&options)?;
        let handle = solver.try_addr()?;
        crate::threads::run(|| handle.solve(y, &lsmr))?
            .map_err(|e| err(e.to_string()))
            .and_then(result_to_list)
    })())
}

// Solve multiple response vectors with a persistent solver.
#[extendr]
fn solver_solve_batch_impl(
    solver: ExternalPtr<SolverHandle>,
    y_matrix: RMatrix<f64>,
    options: Robj,
) -> List {
    or_throw((|| -> Result<List> {
        let lsmr = parse_lsmr_options(&options)?;
        let handle = solver.try_addr()?;

        if y_matrix.nrows() != handle.n_obs() {
            return Err(err(format!(
                "Y has {} rows but solver has {} observations",
                y_matrix.nrows(),
                handle.n_obs()
            )));
        }

        let column_refs = matrix_columns(&y_matrix);

        crate::threads::run(|| handle.solve_batch(&column_refs, &lsmr))?
            .map_err(|e| err(e.to_string()))
            .and_then(batch_result_to_list)
    })())
}

// Return the built preconditioner from a persistent solver, or NULL.
#[extendr]
fn solver_preconditioner_impl(
    solver: ExternalPtr<SolverHandle>,
) -> Option<ExternalPtr<Preconditioner>> {
    or_throw((|| -> Result<Option<ExternalPtr<Preconditioner>>> {
        Ok(solver
            .try_addr()?
            .preconditioner()
            .map(|preconditioner| ExternalPtr::new(preconditioner.clone())))
    })())
}

// Number of DOFs (coefficients) in the persistent solver.
#[extendr]
fn solver_n_dofs_impl(solver: ExternalPtr<SolverHandle>) -> i32 {
    or_throw((|| -> Result<i32> {
        usize_to_i32(solver.try_addr()?.n_dofs(), "n_dofs")
    })())
}

// Number of observations in the persistent solver.
#[extendr]
fn solver_n_obs_impl(solver: ExternalPtr<SolverHandle>) -> i32 {
    or_throw((|| -> Result<i32> {
        usize_to_i32(solver.try_addr()?.n_obs(), "n_obs")
    })())
}

extendr_module! {
    mod api;
    fn design_new_impl;
    fn solver_has_escalated_impl;
    fn solve_impl;
    fn solve_batch_impl;
    fn solver_new_impl;
    fn solver_solve_impl;
    fn solver_solve_batch_impl;
    fn solver_preconditioner_impl;
    fn solver_n_dofs_impl;
    fn solver_n_obs_impl;
}
