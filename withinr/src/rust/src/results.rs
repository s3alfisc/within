//! Native result → R list conversions.

use extendr_api::{error::Result, prelude::*};

use within::{BatchSolveResult, CoefficientAddress, CoefficientLayout, SolveResult};

use crate::convert::usize_to_i32;

fn one_based_i32<T>(value: T, name: &str) -> Result<i32>
where
    i32: TryFrom<T>,
{
    i32::try_from(value)
        .map_err(|_| Error::Other(format!("{name} exceeds i32 range")))?
        .checked_add(1)
        .ok_or_else(|| Error::Other(format!("{name} exceeds i32 range")))
}

fn addresses_to_list(addresses: &[CoefficientAddress]) -> Result<List> {
    let term = addresses
        .iter()
        .map(|address| one_based_i32(address.channel.term, "term"))
        .collect::<Result<Vec<_>>>()?;
    let level = addresses
        .iter()
        .map(|address| one_based_i32(address.level, "level"))
        .collect::<Result<Vec<_>>>()?;
    let column = addresses
        .iter()
        .map(|address| one_based_i32(address.channel.column, "column"))
        .collect::<Result<Vec<_>>>()?;
    Ok(list!(term = term, level = level, column = column))
}

fn one_based_arg(value: i32, name: &str) -> Result<usize> {
    let zero_based = value
        .checked_sub(1)
        .ok_or_else(|| Error::Other(format!("{name} must be >= 1")))?;
    usize::try_from(zero_based).map_err(|_| Error::Other(format!("{name} must be >= 1")))
}

#[extendr]
fn layout_n_dofs_impl(layout: ExternalPtr<CoefficientLayout>) -> i32 {
    crate::convert::or_throw((|| -> Result<i32> {
        usize_to_i32(layout.try_addr()?.n_dofs(), "n_dofs")
    })())
}

#[extendr]
fn layout_n_terms_impl(layout: ExternalPtr<CoefficientLayout>) -> i32 {
    crate::convert::or_throw((|| -> Result<i32> {
        usize_to_i32(layout.try_addr()?.n_terms(), "n_terms")
    })())
}

#[extendr]
fn layout_n_levels_impl(layout: ExternalPtr<CoefficientLayout>, term: i32) -> i32 {
    crate::convert::or_throw((|| -> Result<i32> {
        let term = one_based_arg(term, "term")?;
        let value = layout
            .try_addr()?
            .n_levels(term)
            .ok_or_else(|| Error::Other(format!("term {} is out of range", term + 1)))?;
        usize_to_i32(value, "n_levels")
    })())
}

#[extendr]
fn layout_n_columns_impl(layout: ExternalPtr<CoefficientLayout>, term: i32) -> i32 {
    crate::convert::or_throw((|| -> Result<i32> {
        let term = one_based_arg(term, "term")?;
        let value = layout
            .try_addr()?
            .n_columns(term)
            .ok_or_else(|| Error::Other(format!("term {} is out of range", term + 1)))?;
        usize_to_i32(value, "n_columns")
    })())
}

#[extendr]
fn layout_index_impl(
    layout: ExternalPtr<CoefficientLayout>,
    term: i32,
    level: i32,
    column: i32,
) -> i32 {
    crate::convert::or_throw((|| -> Result<i32> {
        let term = one_based_arg(term, "term")?;
        let column = one_based_arg(column, "column")?;
        let level = level
            .checked_sub(1)
            .ok_or_else(|| Error::Other("level must be >= 1".to_string()))?
            .try_into()
            .map_err(|_| Error::Other("level is out of range".to_string()))?;
        let address = CoefficientAddress {
            channel: within::Channel { term, column },
            level,
        };
        let index = layout.try_addr()?.index(address).ok_or_else(|| {
            Error::Other("coefficient address is out of range for this layout".to_string())
        })?;
        one_based_i32(index, "index")
    })())
}

#[extendr]
fn layout_address_impl(layout: ExternalPtr<CoefficientLayout>, index: i32) -> List {
    crate::convert::or_throw((|| -> Result<List> {
        let index = one_based_arg(index, "index")?;
        let address = layout
            .try_addr()?
            .address(index)
            .ok_or_else(|| Error::Other(format!("index {} is out of range", index + 1)))?;
        Ok(list!(
            term = one_based_i32(address.channel.term, "term")?,
            level = one_based_i32(address.level, "level")?,
            column = one_based_i32(address.channel.column, "column")?
        ))
    })())
}

pub(crate) fn result_to_list(result: SolveResult) -> Result<List> {
    let unidentified = addresses_to_list(&result.unidentified)?;
    let warnings = result
        .warnings
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let layout = ExternalPtr::new(result.layout);
    Ok(list!(
        x = result.x,
        unidentified = unidentified,
        warnings = warnings,
        layout = layout,
        demeaned = result.demeaned,
        converged = result.converged,
        iterations = usize_to_i32(result.iterations, "iterations")?,
        residual = result.residual,
        time_total = result.time_total,
        time_setup = result.time_setup,
        time_solve = result.time_solve
    ))
}

pub(crate) fn batch_result_to_list(result: BatchSolveResult) -> Result<List> {
    let n_rhs = result.converged.len();

    let mut x = RMatrix::new(result.n_dofs, n_rhs);
    x.data_mut().copy_from_slice(&result.x);

    let mut demeaned = RMatrix::new(result.n_obs, n_rhs);
    demeaned.data_mut().copy_from_slice(&result.demeaned);

    let iterations = result
        .iterations
        .iter()
        .map(|&value| usize_to_i32(value, "iterations"))
        .collect::<Result<Vec<_>>>()?;
    let unidentified = addresses_to_list(&result.unidentified)?;
    let warnings = result
        .warnings
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let layout = ExternalPtr::new(result.layout);

    Ok(list!(
        x = x,
        unidentified = unidentified,
        warnings = warnings,
        layout = layout,
        demeaned = demeaned,
        converged = result.converged,
        iterations = iterations,
        residual = result.residual,
        time_solve = result.time_solve,
        time_setup = result.time_setup,
        time_total = result.time_total
    ))
}

extendr_module! {
    mod results;
    fn layout_n_dofs_impl;
    fn layout_n_terms_impl;
    fn layout_n_levels_impl;
    fn layout_n_columns_impl;
    fn layout_index_impl;
    fn layout_address_impl;
}
