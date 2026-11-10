assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

assert_equal <- function(a, b, tol = 0, msg = "Values are not equal.") {
  ok <- isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
  if (!ok) stop(msg, call. = FALSE)
}

assert_error <- function(expr, contains = NULL, msg = "Expected an error, but none occurred.") {
  got <- tryCatch({
    force(expr)
    NULL
  }, error = function(e) e)
  if (is.null(got)) stop(msg, call. = FALSE)
  if (!is.null(contains) && !grepl(contains, conditionMessage(got), fixed = TRUE)) {
    stop(sprintf("Expected error containing '%s', got: %s", contains, conditionMessage(got)),
         call. = FALSE)
  }
}

printed <- function(x) paste(utils::capture.output(print(x)), collapse = "\n")

withinr_run_tests <- function(verbose = TRUE) {
  if (!requireNamespace("withinr", quietly = TRUE)) {
    stop("Package 'withinr' must be installed or loaded first.", call. = FALSE)
  }
  P <- withinr::PreconditionerConfig

  cats_2x2 <- matrix(c(1L, 1L, 2L, 2L,
                       1L, 2L, 1L, 2L), ncol = 2)
  y_simple <- c(1.0, 2.0, 3.0, 4.0)

  # -------------------------------------------------------------------------
  # Configuration objects: validated natively, fully resolved, comparable
  # -------------------------------------------------------------------------
  opts <- withinr::LsmrOptions(tol = 1e-10, maxiter = 2000, local_size = 2)
  assert_equal(opts$maxiter, 2000L, msg = "LsmrOptions maxiter not normalized to integer.")
  assert_true(is.integer(opts$local_size), "LsmrOptions local_size not integer.")
  assert_true(is.null(withinr::LsmrOptions()$local_size), "LsmrOptions default local_size not NULL.")
  # Like Python's constructors, these accept any value of the right type; the
  # core judges them when they are used.
  assert_equal(withinr::LsmrOptions(maxiter = 0, local_size = 0)$maxiter, 0L,
               msg = "LsmrOptions rejected maxiter = 0.")
  assert_error(withinr::within_solve(cats_2x2, y_simple, options = withinr::LsmrOptions(tol = -1)),
               "tolerance")
  assert_error(withinr::LsmrOptions(maxiter = 1.5), "maxiter")

  local <- withinr::LocalSolverConfig()
  assert_equal(local$approx_chol$split_merge, 2L, msg = "LocalSolverConfig default approx_chol not resolved.")
  assert_equal(local$schur$mode, "approximate", msg = "Default Schur mode is not approximate.")
  assert_equal(local$scaling$on_failure, "warn", msg = "Default scaling policy is not warn.")
  assert_equal(withinr::Schur$exact()$mode, "exact", msg = "Schur$exact() mode mismatch.")
  assert_equal(withinr::Schur$approximate(withinr::ApproxSchurConfig(split = 3))$config$split, 3L,
               msg = "Schur$approximate() lost its config.")
  assert_error(withinr::ApproxSchurConfig(split = 0), "split")
  assert_error(withinr::ScalingConfig(on_failure = "explode"), "on_failure")
  assert_error(withinr::LocalSolverConfig(schur = "exact"), "Schur$approximate()")

  stall <- withinr::Staleness()
  assert_equal(c(stall$window, stall$threshold), c(4, 0.7), msg = "Staleness defaults mismatch.")
  assert_error(withinr::Staleness(window = 0), "positive")
  assert_error(withinr::Staleness(threshold = 1), "threshold")

  assert_true(identical(P$Adaptive(), P$Adaptive(withinr::LocalSolverConfig(), withinr::Staleness())),
              "Equal Adaptive configs are not identical().")
  assert_true(!identical(P$Adaptive(), P$Additive()), "Adaptive and Additive compare equal.")
  assert_equal(P$Off()$variant, "Off", msg = "Off variant mismatch.")
  assert_true(grepl("PreconditionerConfig$Adaptive(", printed(P$Adaptive()), fixed = TRUE),
              "PreconditionerConfig print mismatch.")
  assert_true(grepl("Schur$exact()", printed(withinr::LocalSolverConfig(schur = withinr::Schur$exact())),
                    fixed = TRUE), "Schur print mismatch.")

  # -------------------------------------------------------------------------
  # One-shot solve: result shape, layout, correctness
  # -------------------------------------------------------------------------
  r <- withinr::within_solve(cats_2x2, y_simple)
  assert_true(inherits(r, "within_solve_result"), "within_solve() result class mismatch.")
  assert_equal(names(r), c("x", "unidentified", "layout", "demeaned", "converged", "iterations",
                           "residual", "time_total", "time_setup", "time_solve"),
               msg = "within_solve() returned unexpected fields.")
  assert_true(is.data.frame(r$unidentified), "unidentified is not a data frame.")
  assert_equal(names(r$unidentified), c("term", "level", "column"), msg = "unidentified columns.")
  assert_true(isTRUE(r$converged), "within_solve() did not converge on the smoke case.")
  assert_true(is.double(r$x) && is.double(r$demeaned), "x / demeaned are not double.")
  assert_true(r$iterations >= 0L && r$residual >= 0 && r$time_total >= 0, "Bad solve metadata.")
  assert_true(identical(coef(r), r$x) && identical(residuals(r), r$demeaned), "coef()/residuals().")
  assert_true(grepl("<within_solve_result>", printed(r), fixed = TRUE), "Result print mismatch.")

  layout <- r$layout
  assert_true(inherits(layout, "within_coefficient_layout"), "layout class mismatch.")
  assert_equal(layout$n_dofs(), length(r$x), msg = "layout$n_dofs() mismatch.")
  assert_equal(layout$n_terms(), 2L, msg = "layout$n_terms() mismatch.")
  assert_equal(c(layout$n_levels(1), layout$n_columns(1)), c(2L, 1L), msg = "layout per-term sizes.")
  addresses <- layout$address(seq_len(layout$n_dofs()))
  assert_equal(layout$index(addresses$term, addresses$level, addresses$column), seq_len(layout$n_dofs()),
               msg = "layout index/address round trip failed.")
  assert_error(layout$n_levels(3), "out of range")
  assert_error(layout$index(1, 9, 1), "out of range")
  assert_error(layout$address(0), "positive")
  assert_true(grepl("CoefficientLayout(n_terms=2", printed(layout), fixed = TRUE), "layout print.")

  d <- r$demeaned
  assert_equal(c(d[1] + d[2], d[3] + d[4], d[1] + d[3], d[2] + d[4]), rep(0, 4), tol = 1e-6,
               msg = "Group means not centered.")
  y_hat <- r$x[layout$index(1, cats_2x2[, 1], 1)] + r$x[layout$index(2, cats_2x2[, 2], 1)]
  assert_equal(y_hat + r$demeaned, y_simple, tol = 1e-6, msg = "y != D x + demeaned.")

  # Preconditioner choices, options as LsmrOptions or a plain list, weights.
  for (precond in list(P$Off(), P$Diagonal(), P$Additive(), P$Adaptive(), "off", "diagonal",
                       "additive", "adaptive")) {
    assert_true(withinr::within_solve(cats_2x2, y_simple, preconditioner = precond)$converged,
                "A preconditioner choice failed to converge.")
  }
  schwarz <- P$Additive(withinr::LocalSolverConfig(
    approx_chol = withinr::ApproxCholConfig(seed = 7, split_merge = 2),
    schur = withinr::Schur$approximate(withinr::ApproxSchurConfig(split = 1)),
    dense_threshold = 0,
    scaling = withinr::ScalingConfig(tolerance = 1e-8, max_iterations = 100, on_failure = "warn")
  ))
  assert_true(withinr::within_solve(cats_2x2, y_simple, preconditioner = schwarz)$converged,
              "Tuned Additive config failed.")
  exact <- P$Additive(withinr::LocalSolverConfig(schur = withinr::Schur$exact()))
  assert_true(withinr::within_solve(cats_2x2, y_simple, preconditioner = exact)$converged,
              "Exact Schur config failed.")
  assert_true(withinr::within_solve(cats_2x2, y_simple, options = opts)$converged, "LsmrOptions failed.")
  assert_true(withinr::within_solve(cats_2x2, y_simple, options = list(tol = 1e-12))$converged,
              "List options failed.")
  assert_true(withinr::within_solve(cats_2x2, y_simple, weights = c(1, 2, 1, 2))$converged,
              "Weighted solve failed.")
  assert_error(withinr::within_solve(cats_2x2, y_simple, preconditioner = "bogus"), "must be one of")
  assert_error(withinr::within_solve(cats_2x2, y_simple, preconditioner = 3), "preconditioner")
  assert_error(withinr::within_solve(cats_2x2, y_simple, options = 3), "options")

  # -------------------------------------------------------------------------
  # Batch solve
  # -------------------------------------------------------------------------
  Y <- cbind(a = y_simple, b = rev(y_simple))
  b <- withinr::within_solve_batch(cats_2x2, Y)
  assert_true(inherits(b, "within_batch_solve_result"), "Batch result class mismatch.")
  assert_true(is.matrix(b$x) && is.matrix(b$demeaned), "Batch x / demeaned are not matrices.")
  assert_equal(dim(b$demeaned), c(4L, 2L), msg = "Batch demeaned shape mismatch.")
  assert_equal(colnames(b$demeaned), c("a", "b"), msg = "Batch colnames not carried over.")
  assert_equal(names(b$converged), c("a", "b"), msg = "Batch per-RHS names not carried over.")
  assert_true(all(b$converged) && b$time_setup >= 0, "Batch metadata mismatch.")
  s1 <- withinr::within_solve(cats_2x2, Y[, 1])
  s2 <- withinr::within_solve(cats_2x2, Y[, 2])
  assert_equal(b$demeaned, cbind(s1$demeaned, s2$demeaned), tol = 1e-6, msg = "Batch != singles.")
  assert_equal(withinr::within_solve_batch(cats_2x2, as.data.frame(Y))$demeaned, b$demeaned,
               tol = 1e-12, msg = "Data-frame Y mismatch.")
  assert_equal(dim(withinr::within_solve_batch(cats_2x2, Y[, 0, drop = FALSE])$x), c(4L, 0L),
               msg = "Empty batch not well shaped.")
  assert_error(withinr::within_solve_batch(cats_2x2, y_simple), "numeric matrix")
  assert_error(withinr::within_solve_batch(cats_2x2, Y[-1, ]), "rows")

  # -------------------------------------------------------------------------
  # R-native design inputs
  # -------------------------------------------------------------------------
  df <- data.frame(firm = factor(c("x", "x", "y", "y")), year = c("a", "b", "a", "b"))
  assert_equal(withinr::within_solve(df, y_simple)$demeaned, s1$demeaned, tol = 1e-8,
               msg = "Data-frame design mismatch.")
  cats_dbl <- matrix(as.double(cats_2x2), ncol = 2)
  assert_equal(withinr::within_solve(cats_dbl, y_simple)$demeaned, s1$demeaned, tol = 1e-8,
               msg = "Numeric matrix design mismatch.")
  assert_equal(withinr::within_solve(cats_2x2 - 1L, y_simple)$demeaned, s1$demeaned, tol = 1e-8,
               msg = "0-based labels mismatch.")
  one_way <- withinr::within_solve(c(1, 1, 2, 2), y_simple)
  assert_equal(one_way$demeaned, c(-0.5, 0.5, -0.5, 0.5), tol = 1e-8, msg = "Vector design mismatch.")

  bad <- cats_2x2
  bad[1, 1] <- NA_integer_
  assert_error(withinr::within_solve(bad, y_simple), "NA")
  assert_error(withinr::within_solve(-cats_2x2, y_simple), "non-negative")
  cats_fractional <- cats_dbl
  cats_fractional[1, 1] <- 1.5
  assert_error(withinr::within_solve(cats_fractional, y_simple), "whole numbers")
  assert_error(withinr::within_solve(list(1, 2), y_simple), "Effect() terms")
  assert_error(withinr::within_solve(cats_2x2, y_simple[-1]), "length")
  assert_error(withinr::within_solve(cats_2x2, "a"), "numeric vector")

  # -------------------------------------------------------------------------
  # Errors: classed conditions, 1-based messages, no leaks, no format strings
  # -------------------------------------------------------------------------
  caught <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) e)
  e <- caught(withinr::within_solve(cats_2x2, c(1, NA, 3, 4)))
  assert_true(inherits(e, "within_input_error") && inherits(e, "within_error"),
              "Native input errors are not classed.")
  assert_true(grepl("position 2", conditionMessage(e), fixed = TRUE), "y position not 1-based.")
  assert_true(inherits(caught(withinr::Effect(1:3, intercept = FALSE)), "within_input_error"),
              "R-side input errors are not classed.")
  assert_error(withinr::within_solve(cats_2x2, y_simple, weights = c(1, -1, 1, 1)), "weight at position 2")
  assert_error(withinr::within_solve_batch(cats_2x2, cbind(y_simple, c(1, 2, Inf, 4))),
               "`Y[, 2]` must be finite")
  na_matrix <- cats_2x2
  na_matrix[2, 2] <- NA
  assert_error(withinr::within_solve(na_matrix, y_simple), "row 2, column 2")
  assert_error(withinr::within_solve(matrix(c(1, 1, 2, 2, 1, 2.5, 1, 2), 4), y_simple),
               "row 2, column 2")
  assert_error(withinr::within_solve(data.frame(), y_simple), "no factor columns")

  # `%` reaches R verbatim (the message is no format string, which used to segfault).
  assert_error(withinr::ScalingConfig(on_failure = "%d %d %s"), "\"%d %d %s\"")
  assert_error(withinr::ScalingConfig(on_failure = "100%"), "\"100%\"")

  # A failing call releases its arguments.
  invisible(gc())
  before <- sum(gc()[, 2])
  big <- numeric(2e6)
  for (i in 1:10) caught(withinr::within_solve(cats_2x2, big))
  rm(big)
  invisible(gc())
  assert_true(sum(gc()[, 2]) - before < 40, "Failing native calls leak their arguments.")

  # -------------------------------------------------------------------------
  # Build warnings: classed, 1-based, each emitted exactly once
  # -------------------------------------------------------------------------
  f <- c(0, 0, 0, 1, 1, 1)
  g <- c(0, 1, 2, 0, 1, 2)
  z_signed <- c(-2, 1, 1, -1, -1, 2)
  signed <- list(withinr::Effect(f, slopes = z_signed), withinr::Effect(g))
  uncertified <- P$Additive(withinr::LocalSolverConfig(
    scaling = withinr::ScalingConfig(max_iterations = 0)
  ))
  collect_warnings <- function(expr) {
    seen <- list()
    value <- withCallingHandlers(expr, warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = seen)
  }
  one_shot <- collect_warnings(withinr::within_solve(signed, c(1, 3, 2, 5, 4, 6),
                                                     preconditioner = uncertified))
  assert_equal(length(one_shot$warnings), 1L, msg = "One-shot build warning not emitted once.")
  w <- one_shot$warnings[[1]]
  assert_true(inherits(w, "within_build_warning"), "Build warning class mismatch.")
  assert_true(grepl("term 1 column 2 and term 2 column 1", conditionMessage(w), fixed = TRUE),
              "Build warning positions are not 1-based.")
  built <- collect_warnings(withinr::Solver(signed, preconditioner = uncertified))
  assert_equal(length(built$warnings), 1L, msg = "Solver build warning not emitted once.")
  later <- collect_warnings(built$value$solve(c(1, 3, 2, 5, 4, 6)))
  assert_equal(length(later$warnings), 0L, msg = "Solver re-emitted a build warning.")
  assert_true(later$value$converged, "Solve after an uncertified scaling failed.")

  # -------------------------------------------------------------------------
  # Handles and layout edge cases
  # -------------------------------------------------------------------------
  stale <- unserialize(serialize(withinr::Solver(cats_2x2), NULL))
  assert_error(stale$solve(y_simple), "no longer valid")
  assert_error(withinr::within_solve(unserialize(serialize(withinr::Design(cats_2x2), NULL)), y_simple),
               "rebuild it")
  assert_equal(r$layout$index(1, integer(), 1), integer(), msg = "Zero-length index lookup.")
  empty_batch <- withinr::Solver(cats_2x2)$solve_batch(matrix(0, 4, 0))
  assert_true(!grepl("Inf", printed(empty_batch), fixed = TRUE), "Empty batch print is broken.")

  # Character matrices code each column separately, like data frames.
  chars <- cbind(c("a", "a", "b", "b"), c("x", "y", "x", "y"))
  char_fit <- withinr::within_solve(chars, y_simple)
  assert_equal(char_fit$layout$address(1:4)$level, c(1L, 2L, 1L, 2L),
               msg = "Character matrix columns share level codes.")
  assert_equal(char_fit$demeaned, s1$demeaned, tol = 1e-8, msg = "Character matrix design mismatch.")

  # -------------------------------------------------------------------------
  # Persistent Solver and Preconditioner reuse
  # -------------------------------------------------------------------------
  solver <- withinr::Solver(cats_2x2, preconditioner = P$Additive())
  assert_true(inherits(solver, "within_solver"), "Solver() class mismatch.")
  assert_equal(c(solver$n_obs, solver$n_dofs), c(4L, 4L), msg = "Solver sizes mismatch.")
  assert_true(identical(solver$has_escalated, FALSE), "Additive solver reports escalation.")
  assert_error(solver$n_obs <- 1L, "locked")
  p <- solver$preconditioner
  assert_true(inherits(p, "within_preconditioner"), "solver$preconditioner class mismatch.")
  assert_true(grepl("Preconditioner(Additive, n=4)", printed(p), fixed = TRUE), "Preconditioner print.")
  assert_equal(c(p$nrows, p$ncols), c(4L, 4L), msg = "Preconditioner shape mismatch.")
  assert_equal(length(p$apply(rep(1, p$ncols))), p$nrows, msg = "Preconditioner apply length.")
  assert_error(p$apply(1), "length")
  assert_true(p$build_duration_seconds >= 0, "Negative build duration.")
  assert_true(identical(p$config, P$Additive()), "Preconditioner config does not round-trip.")
  assert_true(is.null(withinr::Solver(cats_2x2, preconditioner = "off")$preconditioner),
              "Off solver has a preconditioner.")
  assert_equal(withinr::Solver(cats_2x2, preconditioner = "diagonal")$preconditioner$variant, "Diagonal",
               msg = "Diagonal preconditioner variant mismatch.")

  assert_equal(solver$solve(y_simple)$demeaned, s1$demeaned, tol = 1e-6, msg = "Persistent solve.")
  assert_equal(solver$solve_batch(Y)$demeaned, b$demeaned, tol = 1e-6, msg = "Persistent batch.")
  assert_error(solver$solve_batch(Y[-1, ]), "rows")

  p2 <- withinr::Preconditioner(p$serialize())
  assert_true(identical(p2$config, p$config), "Deserialized config mismatch.")
  assert_equal(p2$build_duration_seconds, p$build_duration_seconds, msg = "Build duration not preserved.")
  assert_equal(withinr::Solver(cats_2x2, preconditioner = p2)$solve(y_simple)$demeaned, s1$demeaned,
               tol = 1e-6, msg = "Preconditioner reuse mismatch.")
  assert_error(withinr::Preconditioner(as.raw(c(1, 2, 3))), "deserialize")
  assert_error(withinr::Preconditioner("bytes"), "raw vector")

  # -------------------------------------------------------------------------
  # Larger fixture: preconditioners agree, caller order kept, base R oracle
  # -------------------------------------------------------------------------
  set.seed(20260704)
  n <- 180L
  cats <- cbind(sample.int(29L, n, TRUE), sample.int(17L, n, TRUE), sample.int(11L, n, TRUE))
  y <- rnorm(n)
  Y_big <- cbind(y, y_alt = 0.25 * y + rnorm(n))
  weights <- runif(n, min = 0.5, max = 2.0)
  tight <- withinr::LsmrOptions(tol = 1e-10, maxiter = 4000L, local_size = 4L)

  fits <- lapply(list(P$Adaptive(), P$Additive(), P$Diagonal(), P$Off()), function(precond) {
    withinr::within_solve(cats, y, weights, tight, precond)
  })
  for (fit in fits) {
    assert_true(fit$converged, "A preconditioner failed on the larger design.")
    assert_equal(fit$demeaned, fits[[1]]$demeaned, tol = 1e-6, msg = "Preconditioners disagree.")
  }

  D <- do.call(cbind, lapply(seq_len(ncol(cats)), function(j) model.matrix(~ factor(cats[, j]) - 1)))
  oracle <- lm.wfit(D, Y_big, w = weights)$residuals
  batch_fresh <- withinr::within_solve_batch(cats, Y_big, weights, tight)
  assert_equal(batch_fresh$demeaned, oracle, tol = 1e-7, msg = "Weighted demeaning != lm.wfit.")

  pre_big <- withinr::Solver(cats, weights = weights, preconditioner = "additive")$preconditioner
  reuse <- withinr::within_solve_batch(cats, Y_big, weights, tight, preconditioner = pre_big)
  assert_equal(reuse$demeaned, oracle, tol = 1e-7, msg = "Batch preconditioner reuse mismatch.")

  ord <- order(cats[, 1], cats[, 2], cats[, 3])
  sorted <- withinr::within_solve_batch(cats[ord, ], Y_big[ord, ], weights[ord], tight)
  assert_equal(batch_fresh$demeaned[ord, ], sorted$demeaned, tol = 1e-6,
               msg = "Unsorted input did not return results in caller order.")

  # Designs own their inputs, keep gappy caller labels, and are reusable.
  gappy <- cats * 100L
  design <- withinr::Design(gappy)
  assert_true(inherits(design, "within_design"), "Design() class mismatch.")
  assert_equal(c(design$n_obs, design$n_dofs), c(n, ncol(D)), msg = "Design sizes.")
  assert_true(identical(withinr::Design(design), design), "Design(Design) is not a pass-through.")
  gappy[, ] <- 1L
  invisible(gc())
  reused <- withinr::within_solve_batch(design, Y_big, weights, tight)
  assert_equal(reused$demeaned, oracle, tol = 1e-7, msg = "Persistent design mismatch.")
  assert_equal(reused$layout$address(1)$level, 100L, msg = "Gappy labels not preserved.")
  assert_equal(reused$layout$index(1, 100, 1), 1L, msg = "Gappy label index mismatch.")
  assert_error(reused$layout$index(1, 1, 1), "out of range")
  assert_equal(withinr::Solver(design)$solve(y, tight)$demeaned, lm.fit(D, y)$residuals, tol = 1e-7,
               msg = "Solver on a persistent design mismatch.")

  # Adaptive: force a stall to exercise the deferred native Schwarz build.
  adaptive <- withinr::Solver(design, weights = weights,
                              preconditioner = P$Adaptive(stall = withinr::Staleness(1, 0)))
  assert_true(!adaptive$has_escalated, "Adaptive built Schwarz eagerly.")
  assert_equal(adaptive$preconditioner$config$variant, "Adaptive", msg = "Ladder base config mismatch.")
  adaptive_result <- adaptive$solve(y, tight)
  assert_true(adaptive_result$converged, "Adaptive failed to converge.")
  assert_true(adaptive$has_escalated, "Forced adaptive escalation did not occur.")
  assert_equal(adaptive_result$demeaned, oracle[, 1], tol = 1e-7, msg = "Adaptive result mismatch.")
  assert_equal(adaptive$preconditioner$variant, "Additive", msg = "Escalated map is not Additive.")

  # -------------------------------------------------------------------------
  # Varying slopes, checked against independently built dummy designs
  # -------------------------------------------------------------------------
  z <- rnorm(n)
  z2 <- rnorm(n)
  D1 <- model.matrix(~ factor(cats[, 1]) - 1)
  D2 <- model.matrix(~ factor(cats[, 2]) - 1)

  slope_design <- list(withinr::Effect(cats[, 1], slopes = z), withinr::Effect(cats[, 2]))
  slope_oracle <- lm.wfit(cbind(D1, D1 * z, D2), Y_big, weights)$residuals
  slope_result <- withinr::within_solve_batch(slope_design, Y_big, weights, tight)
  assert_true(all(slope_result$converged), "Varying slopes failed to converge.")
  assert_equal(slope_result$demeaned, slope_oracle, tol = 1e-7, msg = "Varying slopes != lm.wfit.")
  assert_equal(slope_result$layout$n_columns(1), 2L, msg = "Slope term column count mismatch.")
  slope_solver <- withinr::Solver(slope_design, weights = weights)
  assert_equal(slope_solver$solve_batch(Y_big, tight)$demeaned, slope_oracle, tol = 1e-7,
               msg = "Persistent varying-slope solve mismatch.")

  # Two slopes from a matrix, and a slopes-only term (no intercept).
  multi <- list(withinr::Effect(cats[, 1], slopes = cbind(z, z2)),
                withinr::Effect(cats[, 2], intercept = FALSE, slopes = list(z2)),
                withinr::Effect(cats[, 3]))
  D3 <- model.matrix(~ factor(cats[, 3]) - 1)
  multi_oracle <- lm.fit(cbind(D1, D1 * z, D1 * z2, D2 * z2, D3), y)$residuals
  multi_fit <- withinr::within_solve(multi, y, options = tight)
  assert_true(multi_fit$converged, "Multi-slope design failed to converge.")
  assert_equal(multi_fit$demeaned, multi_oracle, tol = 1e-7, msg = "Multi-slope design != lm.fit.")
  assert_equal(sapply(1:3, multi_fit$layout$n_columns), c(3L, 1L, 1L), msg = "Multi-slope layout.")
  assert_equal(withinr::within_solve(withinr::Design(multi), y, options = tight)$demeaned, multi_oracle,
               tol = 1e-7, msg = "Design(list of Effect) mismatch.")

  # A Design owns its slopes: the R vectors may change or be collected.
  z_copy <- z + 0
  owned <- withinr::Design(list(withinr::Effect(cats[, 1], slopes = z_copy), withinr::Effect(cats[, 2])))
  z_copy[] <- 0
  rm(z_copy)
  invisible(gc())
  assert_equal(withinr::within_solve_batch(owned, Y_big, weights, tight)$demeaned, slope_oracle,
               tol = 1e-7, msg = "Design does not own its slope columns.")

  # A level whose slope covariate is identically zero is unidentified; its slot is 0.
  lv <- c(1, 1, 2, 2, 3, 3)
  zero_slope <- withinr::within_solve(list(withinr::Effect(lv, FALSE, c(0, 0, 1, 2, 3, 1))),
                                      c(1, 2, 3, 4, 5, 7))
  assert_equal(zero_slope$unidentified, data.frame(term = 1L, level = 1L, column = 1L),
               msg = "Unidentified direction not reported.")
  assert_equal(zero_slope$x[zero_slope$layout$index(1, 1, 1)], 0, msg = "Unidentified slot not 0.")

  assert_error(withinr::Effect(lv, slopes = 1:3), "slope 1 has 3 values, expected 6")
  assert_error(withinr::Effect(lv, intercept = FALSE), "intercept or at least one slope")
  assert_error(withinr::Effect(lv, slopes = c(1, 2, NA, 4, 5, 6)), "slope 1 value at position 3")
  assert_error(withinr::Effect(lv, intercept = NA), "TRUE or FALSE")

  # Native handles do not survive serialization.
  restored <- unserialize(serialize(design, NULL))
  assert_error(withinr::Solver(restored))

  if (verbose) message("withinr tests: OK")
  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  withinr_run_tests(verbose = TRUE)
}
