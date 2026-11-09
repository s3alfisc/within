//! Native result → R list conversions and the coefficient-layout handle.
//!
//! Mirrors `crates/within-py/src/results.rs`. Term, column, and flat `x`
//! indices are 1-based in R; levels are the caller's own factor labels.

use extendr_api::{error::Result, prelude::*};

use within::{
    BatchSolveResult, BuildWarning, Channel, CoefficientAddress, CoefficientLayout, SolveResult,
};

use crate::convert::{err, one_based, or_throw, usize_to_i32};

fn label_to_i32(level: u32) -> Result<i32> {
    i32::try_from(level).map_err(|_| err("level label exceeds R's integer range"))
}

/// Columns `term`, `level`, `column` (the R wrapper turns them into a data frame).
fn addresses_to_list(addresses: &[CoefficientAddress]) -> Result<List> {
    let mut term = Vec::with_capacity(addresses.len());
    let mut level = Vec::with_capacity(addresses.len());
    let mut column = Vec::with_capacity(addresses.len());
    for address in addresses {
        term.push(one_based(address.channel.term, "term")?);
        level.push(label_to_i32(address.level)?);
        column.push(one_based(address.channel.column, "column")?);
    }
    Ok(list!(term = term, level = level, column = column))
}

pub(crate) fn warning_strings(warnings: &[BuildWarning]) -> Vec<String> {
    warnings.iter().map(ToString::to_string).collect()
}

/// `warnings` are re-emitted as R warnings by the R wrapper, then dropped.
pub(crate) fn result_to_list(result: SolveResult, warnings: Vec<String>) -> Result<List> {
    Ok(list!(
        x = result.x,
        unidentified = addresses_to_list(&result.unidentified)?,
        layout = ExternalPtr::new(result.layout),
        demeaned = result.demeaned,
        converged = result.converged,
        iterations = usize_to_i32(result.iterations, "iterations")?,
        residual = result.residual,
        time_total = result.time_total,
        time_setup = result.time_setup,
        time_solve = result.time_solve,
        warnings = warnings
    ))
}

fn column_major(nrow: usize, ncol: usize, data: &[f64]) -> RMatrix<f64> {
    let mut matrix = RMatrix::new(nrow, ncol);
    matrix.data_mut().copy_from_slice(data);
    matrix
}

pub(crate) fn batch_result_to_list(
    result: BatchSolveResult,
    warnings: Vec<String>,
) -> Result<List> {
    let n_rhs = result.converged.len();
    let iterations = result
        .iterations
        .iter()
        .map(|&value| usize_to_i32(value, "iterations"))
        .collect::<Result<Vec<_>>>()?;
    Ok(list!(
        x = column_major(result.n_dofs, n_rhs, &result.x),
        unidentified = addresses_to_list(&result.unidentified)?,
        layout = ExternalPtr::new(result.layout),
        demeaned = column_major(result.n_obs, n_rhs, &result.demeaned),
        converged = result.converged,
        iterations = iterations,
        residual = result.residual,
        time_solve = result.time_solve,
        time_setup = result.time_setup,
        time_total = result.time_total,
        warnings = warnings
    ))
}

// ---------------------------------------------------------------------------
// Coefficient layout
// ---------------------------------------------------------------------------

// `n_dofs`, `n_terms`, and per-term `n_levels` / `n_columns`, read once.
#[extendr]
fn layout_info_impl(layout: ExternalPtr<CoefficientLayout>) -> List {
    or_throw((|| -> Result<List> {
        let layout = layout.try_addr()?;
        let n_terms = layout.n_terms();
        let per_term = |f: fn(&CoefficientLayout, usize) -> Option<usize>, name| {
            (0..n_terms)
                .map(|t| usize_to_i32(f(layout, t).unwrap_or(0), name))
                .collect::<Result<Vec<_>>>()
        };
        Ok(list!(
            n_dofs = usize_to_i32(layout.n_dofs(), "n_dofs")?,
            n_terms = usize_to_i32(n_terms, "n_terms")?,
            n_levels = per_term(CoefficientLayout::n_levels, "n_levels")?,
            n_columns = per_term(CoefficientLayout::n_columns, "n_columns")?
        ))
    })())
}

fn one_based_arg(value: i32, name: &str) -> Result<usize> {
    if value == i32::MIN || value < 1 {
        return Err(err(format!("`{name}` must be a positive integer")));
    }
    Ok(value as usize - 1)
}

// Flat 1-based `x` index of each `(term, level, column)`; the inputs have equal length.
#[extendr]
fn layout_index_impl(
    layout: ExternalPtr<CoefficientLayout>,
    term: &[i32],
    level: &[i32],
    column: &[i32],
) -> Vec<i32> {
    or_throw((|| -> Result<Vec<i32>> {
        let layout = layout.try_addr()?;
        term.iter()
            .zip(level)
            .zip(column)
            .map(|((&t, &l), &c)| {
                if l == i32::MIN || l < 0 {
                    return Err(err("`level` must be a non-negative label"));
                }
                let at = CoefficientAddress {
                    channel: Channel {
                        term: one_based_arg(t, "term")?,
                        column: one_based_arg(c, "column")?,
                    },
                    level: l as u32,
                };
                let index = layout.index(at).ok_or_else(|| {
                    err(format!(
                        "coefficient address (term={t}, level={l}, column={c}) out of range"
                    ))
                })?;
                one_based(index, "index")
            })
            .collect()
    })())
}

// `(term, level, column)` of each 1-based flat `x` index.
#[extendr]
fn layout_address_impl(layout: ExternalPtr<CoefficientLayout>, index: &[i32]) -> List {
    or_throw((|| -> Result<List> {
        let layout = layout.try_addr()?;
        let addresses = index
            .iter()
            .map(|&i| {
                let at = one_based_arg(i, "index")?;
                layout.address(at).ok_or_else(|| {
                    err(format!(
                        "x index {i} out of range (n_dofs={})",
                        layout.n_dofs()
                    ))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        addresses_to_list(&addresses)
    })())
}

extendr_module! {
    mod results;
    fn layout_info_impl;
    fn layout_index_impl;
    fn layout_address_impl;
}
