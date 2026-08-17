//! Behavior of the `Adaptive` (diagonal→Schwarz) preconditioner strategy (#260).

use within::config::{LocalSolverConfig, ReductionStrategy};
use within::{LsmrOptions, PreconditionerConfig, Solver, Staleness};

#[path = "common/orchestrate_helpers.rs"]
mod common;

/// Escalates after any single non-vanishing contraction, so a handoff is deterministic.
fn eager_stall() -> Staleness {
    Staleness::try_new(1, 0.0).expect("valid staleness")
}

fn adaptive(stall: Staleness) -> PreconditionerConfig {
    PreconditionerConfig::Adaptive {
        local_solver: LocalSolverConfig::default(),
        reduction: ReductionStrategy::Auto,
        stall,
    }
}

fn additive() -> PreconditionerConfig {
    PreconditionerConfig::Additive {
        local_solver: LocalSolverConfig::default(),
        reduction: ReductionStrategy::Auto,
    }
}

fn params() -> LsmrOptions {
    LsmrOptions {
        tol: 1e-8,
        maxiter: 1000,
        ..Default::default()
    }
}

#[test]
fn escalation_is_visible_and_answer_matches_additive() {
    let y = common::make_deterministic_y(&common::make_test_design());

    let adaptive_solver =
        Solver::new(common::make_test_design(), None, adaptive(eager_stall())).expect("solver");
    let adaptive_result = adaptive_solver
        .solve(&y, &params())
        .expect("adaptive solve");
    let additive_solver =
        Solver::new(common::make_test_design(), None, additive()).expect("solver");
    let additive_result = additive_solver
        .solve(&y, &params())
        .expect("additive solve");

    // Visible in the solver's lifecycle...
    assert!(adaptive_solver.has_escalated(), "eager stall must hand off");
    assert!(!additive_solver.has_escalated());

    // ...invisible in the answer: the fitted values `Dx` agree with a cold Schwarz solve
    // of the same system. Raw coefficients agree only up to `null(D)` (the reference-level
    // indeterminacy of a rank-deficient FE model): a warm-started resume shares the unique
    // fitted values but its min-norm correction picks a different null-space representative.
    common::assert_converged_with_small_residual(&adaptive_result, 1e-6);
    common::assert_normal_equations_satisfied(
        &common::test_categories(),
        None,
        &y,
        &adaptive_result,
        1e-6,
    );
    common::assert_solutions_close(&adaptive_result.demeaned, &additive_result.demeaned, 1e-6);
}

#[test]
fn weighted_answer_matches_additive() {
    let cats = vec![vec![0, 1, 0, 1, 2], vec![0, 0, 1, 1, 0]];
    let weights = vec![1.0, 2.0, 1.5, 0.5, 3.0];
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let solve_with = |precond: PreconditionerConfig| {
        Solver::new(
            common::make_design(cats.clone()).expect("design"),
            Some(weights.clone()),
            precond,
        )
        .expect("solver")
        .solve(&y, &params())
        .expect("solve")
    };

    let adaptive_solver = Solver::new(
        common::make_design(cats.clone()).expect("design"),
        Some(weights.clone()),
        adaptive(eager_stall()),
    )
    .expect("solver");
    let adaptive_result = adaptive_solver.solve(&y, &params()).expect("solve");
    let additive_result = solve_with(additive());

    assert!(adaptive_solver.has_escalated());
    common::assert_normal_equations_satisfied(
        &cats,
        Some(weights.as_slice()),
        &y,
        &adaptive_result,
        1e-6,
    );
    common::assert_solutions_close(&adaptive_result.demeaned, &additive_result.demeaned, 1e-6);
}

#[test]
fn without_a_target_it_behaves_like_diagonal() {
    // A single factor has no cross-factor pair, so there is nothing to escalate to.
    let design = || common::make_design(vec![vec![0, 1, 2, 0, 1, 2]]).expect("design");
    let y = common::make_deterministic_y(&design());

    let solver = Solver::new(design(), None, adaptive(eager_stall())).expect("solver");
    // Two solves: the first probes and caches "no target"; the second takes the fast path.
    let first = solver.solve(&y, &params()).expect("first solve");
    let second = solver.solve(&y, &params()).expect("second solve");
    let diagonal = Solver::new(design(), None, &PreconditionerConfig::Diagonal)
        .expect("solver")
        .solve(&y, &params())
        .expect("diagonal solve");

    assert!(
        !solver.has_escalated(),
        "no factor-pair target means no escalation"
    );
    common::assert_solutions_close(&first.demeaned, &diagonal.demeaned, 1e-6);
    common::assert_solutions_close(&second.demeaned, &diagonal.demeaned, 1e-6);
}

#[test]
fn total_iteration_budget_is_honored_across_rungs() {
    let y = common::make_deterministic_y(&common::make_test_design());
    let budget = LsmrOptions {
        tol: 1e-12,
        maxiter: 5,
        ..Default::default()
    };
    let result = Solver::new(common::make_test_design(), None, adaptive(eager_stall()))
        .expect("solver")
        .solve(&y, &budget)
        .expect("solve");

    assert!(
        result.iterations <= budget.maxiter,
        "ladder used {} iterations, exceeding the budget of {}",
        result.iterations,
        budget.maxiter
    );
}

#[test]
fn batch_escalation_builds_schwarz_once() {
    let y0 = common::make_deterministic_y(&common::make_test_design());
    let y1: Vec<f64> = y0.iter().map(|v| v * 0.5 + 0.3).collect();

    let solver =
        Solver::new(common::make_test_design(), None, adaptive(eager_stall())).expect("solver");
    let result = solver
        .solve_batch(&[&y0[..], &y1[..]], &params())
        .expect("batch solve");

    assert!(result.converged.iter().all(|&c| c));
    assert!(solver.has_escalated());
    // The single deferred build has landed, so the active preconditioner is Schwarz.
    assert_eq!(
        solver.preconditioner().expect("built").variant_name(),
        "Additive"
    );
}
