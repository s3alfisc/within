//! Cross-term collinearity screen (#281): a slope covariate that is (nearly) a per-level
//! combination of another term's columns adds a null direction whitening cannot see.

use rayon::prelude::*;

use super::{Design, Loading};
use crate::channel::Channel;
use crate::BuildWarning;

/// Residual share below which a covariate counts as reproduced by the other term.
const COLLINEARITY_TOL: f64 = 1e-3;

/// Pivots below this share of the largest diagonal count as a dependent column.
const PIVOT_FLOOR: f64 = 1e-12;

pub(crate) fn detect_collinear_slopes(
    design: &Design<'_>,
    weights: Option<&[f64]>,
) -> Vec<BuildWarning> {
    // Grouped by the reproducing term: its per-level Gram is shared by every candidate.
    (0..design.n_factors())
        .into_par_iter()
        .flat_map_iter(|term| {
            let candidates: Vec<(Channel, u32)> = (0..design.n_factors())
                .filter(|&t| t != term)
                .flat_map(|t| design.channels(t))
                .filter_map(|slope| match design.loading(slope) {
                    Loading::Covariate(column) => Some((slope, column)),
                    Loading::Constant => None,
                })
                .collect();
            let residuals = relative_residuals(design, weights, &candidates, term);
            candidates
                .into_iter()
                .zip(residuals)
                .filter_map(move |((slope, _), relative_residual)| {
                    (relative_residual <= COLLINEARITY_TOL).then_some(
                        BuildWarning::CollinearSlopeCovariate {
                            slope,
                            term,
                            relative_residual,
                        },
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Share of each candidate covariate's (weighted) sum of squares outside `term`'s column span.
fn relative_residuals(
    design: &Design<'_>,
    weights: Option<&[f64]>,
    candidates: &[(Channel, u32)],
    term: usize,
) -> Vec<f64> {
    let s = candidates.len();
    if s == 0 {
        return Vec::new();
    }
    let targets: Vec<&[f64]> = candidates
        .iter()
        .map(|&(_, column)| design.frame.loading_column(column as usize))
        .collect();
    let levels = design.frame.level_column(term);
    let columns = &design.terms[term].columns;
    let p = columns.len();

    let covariate_slots: Vec<(usize, &[f64])> = columns
        .iter()
        .enumerate()
        .filter_map(|(slot, c)| {
            c.covariate()
                .map(|&k| (slot, design.frame.loading_column(k as usize)))
        })
        .collect();

    let mut per_level_gram = vec![0.0f64; design.terms[term].n_levels * p * p];
    let mut per_level_cross = vec![0.0f64; design.terms[term].n_levels * s * p];
    let mut ss_total = vec![0.0f64; s];
    let mut x = vec![1.0f64; p];
    for (obs, &level) in levels.iter().enumerate() {
        let w = weights.map_or(1.0, |w| w[obs]);
        for &(slot, loadings) in &covariate_slots {
            x[slot] = loadings[obs];
        }
        let gram = &mut per_level_gram[level as usize * p * p..][..p * p];
        for q in 0..p {
            let wx = w * x[q];
            for r in 0..p {
                gram[q * p + r] += wx * x[r];
            }
        }
        let cross = &mut per_level_cross[level as usize * s * p..][..s * p];
        for (j, target) in targets.iter().enumerate() {
            let value = target[obs];
            ss_total[j] += w * value * value;
            for q in 0..p {
                cross[j * p + q] += w * x[q] * value;
            }
        }
    }

    let mut ss_projected = vec![0.0f64; s];
    let mut active = Vec::with_capacity(p);
    for (gram, crosses) in per_level_gram
        .chunks_exact_mut(p * p)
        .zip(per_level_cross.chunks_exact_mut(s * p))
    {
        projected_ss(gram, crosses, p, &mut ss_projected, &mut active);
    }
    ss_total
        .iter()
        .zip(&ss_projected)
        .map(|(&total, &projected)| {
            if total <= 0.0 {
                1.0
            } else {
                ((total - projected) / total).max(0.0)
            }
        })
        .collect()
}

/// Adds `cᵀ G⁺ c` per cross vector to `ss`, by pivoted outer-product Cholesky of one level's Gram.
fn projected_ss(
    gram: &mut [f64],
    crosses: &mut [f64],
    p: usize,
    ss: &mut [f64],
    active: &mut Vec<usize>,
) {
    let floor = PIVOT_FLOOR * (0..p).map(|q| gram[q * p + q]).fold(0.0f64, f64::max);
    active.clear();
    active.extend(0..p);
    while let Some((position, &pivot)) = active
        .iter()
        .enumerate()
        .max_by(|a, b| gram[a.1 * p + a.1].total_cmp(&gram[b.1 * p + b.1]))
    {
        let diagonal = gram[pivot * p + pivot];
        if diagonal <= floor {
            break;
        }
        active.swap_remove(position);
        for (contribution, cross) in ss.iter_mut().zip(crosses.chunks_exact(p)) {
            *contribution += cross[pivot] * cross[pivot] / diagonal;
        }
        for &row in active.iter() {
            let multiplier = gram[row * p + pivot] / diagonal;
            for cross in crosses.chunks_exact_mut(p) {
                cross[row] -= multiplier * cross[pivot];
            }
            for &col in active.iter() {
                gram[row * p + col] -= multiplier * gram[pivot * p + col];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Effect;

    fn warn_pairs(design: &Design<'_>, weights: Option<&[f64]>) -> Vec<(Channel, usize)> {
        detect_collinear_slopes(design, weights)
            .into_iter()
            .map(|w| match w {
                BuildWarning::CollinearSlopeCovariate { slope, term, .. } => (slope, term),
                other => panic!("unexpected warning {other:?}"),
            })
            .collect()
    }

    fn two_factor_levels(n: usize) -> (Vec<u32>, Vec<u32>) {
        let a = (0..n).map(|i| (i % 40) as u32).collect();
        let b = (0..n).map(|i| ((i / 40) % 25) as u32).collect();
        (a, b)
    }

    use crate::test_rng::pseudo_noise;

    #[test]
    fn shared_covariate_across_two_terms_warns_both_ways() {
        let n = 4000;
        let (a, b) = two_factor_levels(n);
        let z = pseudo_noise(n, 3);
        let design = Design::new(vec![
            Effect::new(&a, true, [&z[..]]).unwrap(),
            Effect::new(&b, true, [&z[..]]).unwrap(),
        ])
        .unwrap();
        let pairs = warn_pairs(&design, None);
        assert!(
            pairs.contains(&(Channel { term: 0, column: 1 }, 1)),
            "{pairs:?}"
        );
        assert!(
            pairs.contains(&(Channel { term: 1, column: 1 }, 0)),
            "{pairs:?}"
        );
    }

    #[test]
    fn per_level_affine_transform_of_a_shared_covariate_warns() {
        let n = 4000;
        let (a, b) = two_factor_levels(n);
        let z = pseudo_noise(n, 5);
        let scaled: Vec<f64> = z
            .iter()
            .zip(&b)
            .map(|(&v, &level)| v * (1.0 + level as f64) - 0.25 * level as f64)
            .collect();
        let design = Design::new(vec![
            Effect::new(&a, true, [&z[..]]).unwrap(),
            Effect::new(&b, true, [&scaled[..]]).unwrap(),
        ])
        .unwrap();
        let pairs = warn_pairs(&design, None);
        assert!(
            pairs.contains(&(Channel { term: 0, column: 1 }, 1)),
            "{pairs:?}"
        );
    }

    #[test]
    fn factor_measurable_covariate_warns_against_the_intercept_term() {
        let n = 4000;
        let (a, b) = two_factor_levels(n);
        let of_b: Vec<f64> = b.iter().map(|&level| level as f64 + 1.0).collect();
        let design = Design::new(vec![
            Effect::new(&a, true, [&of_b[..]]).unwrap(),
            Effect::new(&b, true, []).unwrap(),
        ])
        .unwrap();
        let pairs = warn_pairs(&design, None);
        assert_eq!(pairs, vec![(Channel { term: 0, column: 1 }, 1)]);
    }

    #[test]
    fn independent_covariates_stay_silent() {
        let n = 4000;
        let (a, b) = two_factor_levels(n);
        let z1 = pseudo_noise(n, 7);
        let z2 = pseudo_noise(n, 11);
        let design = Design::new(vec![
            Effect::new(&a, true, [&z1[..]]).unwrap(),
            Effect::new(&b, true, [&z2[..]]).unwrap(),
        ])
        .unwrap();
        assert_eq!(warn_pairs(&design, None), Vec::new());
    }

    #[test]
    fn zero_weight_rows_are_excluded_from_the_screen() {
        // The covariates disagree only on rows whose weight is zero (row removal).
        let n = 4000;
        let (a, b) = two_factor_levels(n);
        let z = pseudo_noise(n, 13);
        let mut spoiled = z.clone();
        let mut weights = vec![1.0; n];
        for i in (0..n).step_by(7) {
            spoiled[i] += 10.0;
            weights[i] = 0.0;
        }
        let design = Design::new(vec![
            Effect::new(&a, true, [&z[..]]).unwrap(),
            Effect::new(&b, true, [&spoiled[..]]).unwrap(),
        ])
        .unwrap();
        // The screen reads frame order, so caller weights go through the locality permutation.
        let weights = design.permute_obs_in(&weights).into_owned();
        assert!(!warn_pairs(&design, Some(&weights)).is_empty());
        assert_eq!(warn_pairs(&design, None), Vec::new());
    }
}
