validate_levels <- function(levels, name = "levels") {
  if (!is.numeric(levels)) {
    stop(sprintf("`%s` must be numeric", name), call. = FALSE)
  }
  if (anyNA(levels) || any(!is.finite(levels)) || any(levels != trunc(levels))) {
    stop(sprintf("`%s` must contain finite integers with no NAs", name), call. = FALSE)
  }
  if (any(levels < 1) || any(levels > .Machine$integer.max)) {
    stop(sprintf("`%s` must contain 1-based labels in integer range", name), call. = FALSE)
  }
  as.integer(levels)
}

validate_design <- function(design) {
  if (inherits(design, "within_design")) return(design)
  if (is.matrix(design)) {
    return(matrix(
      validate_levels(design, "design"),
      nrow = nrow(design),
      ncol = ncol(design)
    ))
  }
  if (is.list(design) && all(vapply(design, inherits, logical(1), "within_effect"))) {
    return(design)
  }
  stop("`design` must be a category matrix or a list of Effect() terms", call. = FALSE)
}

#' Fixed-effect term
#'
#' A term contains 1-based factor labels, an optional intercept column, and
#' zero or more varying-slope covariates.
#'
#' @param levels Numeric vector of positive integer factor labels.
#' @param intercept Whether to include the factor intercept.
#' @param slopes `NULL`, a numeric vector, a numeric matrix (one slope per
#'   column), or a list of numeric vectors.
#' @return A term accepted in the `design` list passed to [within_solve()],
#'   [within_solve_batch()], or [Solver()].
#' @export
Effect <- function(levels, intercept = TRUE, slopes = NULL) {
  levels <- validate_levels(levels)
  if (length(intercept) != 1L || is.na(intercept)) {
    stop("`intercept` must be TRUE or FALSE", call. = FALSE)
  }
  intercept <- as.logical(intercept)
  if (is.null(slopes)) {
    slopes <- list()
  } else if (is.matrix(slopes)) {
    slopes <- lapply(seq_len(ncol(slopes)), function(column) {
      as.double(slopes[, column])
    })
  } else if (is.numeric(slopes) && !is.list(slopes)) {
    slopes <- list(as.double(slopes))
  } else if (is.list(slopes)) {
    slopes <- lapply(slopes, function(slope) {
      if (!is.numeric(slope)) {
        stop("each slope must be numeric", call. = FALSE)
      }
      as.double(slope)
    })
  } else {
    stop("`slopes` must be NULL, numeric, a matrix, or a list", call. = FALSE)
  }
  if (any(lengths(slopes) != length(levels))) {
    stop("each slope must have the same length as `levels`", call. = FALSE)
  }
  structure(
    list(levels = levels, intercept = intercept, slopes = slopes),
    class = "within_effect"
  )
}

validate_weights <- function(weights) {
  if (is.null(weights)) {
    NULL
  } else {
    as.double(weights)
  }
}

new_coefficient_layout <- function(ptr) {
  layout <- new.env(parent = emptyenv())
  layout$ptr <- ptr
  layout$n_dofs <- layout_n_dofs_impl(ptr)
  layout$n_terms <- layout_n_terms_impl(ptr)
  layout$n_levels <- function(term) {
    layout_n_levels_impl(ptr, as.integer(term))
  }
  layout$n_columns <- function(term) {
    layout_n_columns_impl(ptr, as.integer(term))
  }
  layout$index <- function(term, level, column) {
    layout_index_impl(
      ptr,
      as.integer(term),
      as.integer(level),
      as.integer(column)
    )
  }
  layout$address <- function(index) {
    layout_address_impl(ptr, as.integer(index))
  }
  class(layout) <- "within_coefficient_layout"
  layout
}

decorate_result <- function(result) {
  result$layout <- new_coefficient_layout(result$layout)
  result
}

#' Modified LSMR solver options
#'
#' Uses Modified Golub-Kahan bidiagonalization to solve the least-squares
#' problem directly. The preconditioner approximates `A^T A` and is applied as
#' one `M^{-1}` solve per iteration.
#'
#' Values are validated by the Rust bridge when the options are used.
#'
#' @param tol Positive finite convergence tolerance. Default `1e-8`.
#' @param maxiter Maximum number of LSMR iterations. Default `1000L`.
#' @param local_size Optional window size for modified Gram-Schmidt
#'   reorthogonalization, or `NULL` to use the short recurrence.
#' @return A solver options object accepted by [within_solve()],
#'   [within_solve_batch()], and persistent solver methods.
#' @export
LsmrOptions <- function(tol = 1e-8, maxiter = 1000L, local_size = NULL) {
  structure(
    list(tol = tol, maxiter = maxiter, local_size = local_size),
    class = "within_lsmr_options"
  )
}

#' Preconditioner shortcut values
#'
#' Use `PreconditionerConfig$Additive`, `PreconditionerConfig$Off`, or
#' `PreconditionerConfig$Diagonal`, or `PreconditionerConfig$Adaptive` as the
#' `preconditioner` argument. `Adaptive` is the library default.
#' `Additive` builds the default additive Schwarz preconditioner, `Off`
#' disables preconditioning, and `Diagonal` uses diagonal/Jacobi scaling.
#'
#' @export
PreconditionerConfig <- structure(
  list(Adaptive = "adaptive", Additive = "additive", Off = "off", Diagonal = "diagonal"),
  class = "within_preconditioner_config"
)

#' Reduction strategy shortcut values
#'
#' Use these values inside [AdditiveSchwarz()].
#'
#' @export
ReductionStrategy <- structure(
  list(
    Auto = "auto",
    AtomicScatter = "atomic_scatter",
    ParallelReduction = "parallel_reduction"
  ),
  class = "within_reduction_strategy"
)

#' Scaling failure policy shortcut values
#'
#' Use these values in [ScalingConfig()].
#'
#' @export
ScalingFailure <- structure(
  list(Warn = "warn", Error = "error"),
  class = "within_scaling_failure"
)

#' Approximate Cholesky configuration
#'
#' @param seed Non-negative integer random seed.
#' @param split_merge Optional positive integer split/merge parameter, or `NULL`.
#' @return A local-solver configuration object.
#' @export
ApproxCholConfig <- function(seed = 0, split_merge = NULL) {
  structure(
    list(seed = seed, split_merge = split_merge),
    class = "within_approx_chol_config"
  )
}

#' Approximate Schur complement configuration
#'
#' @param seed Non-negative integer random seed.
#' @param split Positive integer edge split factor.
#' @return A local-solver configuration object.
#' @export
ApproxSchurConfig <- function(seed = 0, split = 1L) {
  structure(
    list(seed = seed, split = split),
    class = "within_approx_schur_config"
  )
}

#' Schur reduction modes
#'
#' Use `Schur$Approximate()` for sampled Schur reduction or `Schur$Exact` for
#' the exact complement. `NULL` in [LocalSolverConfig()] selects the library
#' default (approximate).
#'
#' @export
Schur <- structure(
  list(
    Approximate = function(config = NULL) {
      structure(
        list(config = config),
        class = "within_schur_approximate"
      )
    },
    Exact = structure(list(), class = "within_schur_exact")
  ),
  class = "within_schur_modes"
)

#' Signed-component scaling configuration
#'
#' @param tolerance Non-negative certification tolerance, or `NULL` for the
#'   library default.
#' @param max_iterations Positive reduced-CG iteration budget, or `NULL`.
#' @param on_failure `ScalingFailure$Warn`, `ScalingFailure$Error`, or `NULL`.
#' @return A scaling configuration accepted by [LocalSolverConfig()].
#' @export
ScalingConfig <- function(tolerance = NULL,
                          max_iterations = NULL,
                          on_failure = NULL) {
  structure(
    list(
      tolerance = tolerance,
      max_iterations = max_iterations,
      on_failure = on_failure
    ),
    class = "within_scaling_config"
  )
}

#' Local solver configuration for additive Schwarz subdomains
#'
#' @param approx_chol `NULL` for the library default, or an
#'   [ApproxCholConfig()] object.
#' @param schur `NULL` for the library default, `Schur$Exact`, or an object
#'   created by `Schur$Approximate()`.
#' @param dense_threshold Optional non-negative integer dense Schur threshold.
#' @param scaling `NULL` for the library default, or a [ScalingConfig()] object.
#' @param ridge Optional non-negative local spectral floor.
#' @return A local-solver configuration object accepted by [AdditiveSchwarz()].
#' @export
LocalSolverConfig <- function(approx_chol = NULL,
                              schur = NULL,
                              dense_threshold = NULL,
                              scaling = NULL,
                              ridge = NULL) {
  structure(
    list(
      approx_chol = approx_chol,
      schur = schur,
      dense_threshold = dense_threshold,
      scaling = scaling,
      ridge = ridge
    ),
    class = "within_local_solver_config"
  )
}

#' Additive Schwarz preconditioner configuration
#'
#' @param local_solver `NULL` for the library default, or a
#'   [LocalSolverConfig()] object.
#' @param reduction One of `ReductionStrategy$Auto`,
#'   `ReductionStrategy$AtomicScatter`, or `ReductionStrategy$ParallelReduction`.
#' @return A preconditioner configuration object.
#' @export
AdditiveSchwarz <- function(local_solver = NULL,
                            reduction = ReductionStrategy$Auto) {
  structure(
    list(local_solver = local_solver, reduction = reduction),
    class = "within_additive_schwarz"
  )
}

#' Solve fixed-effects normal equations
#'
#' Computes fixed-effect coefficients by solving the normal equations
#' \eqn{D^T W D x = D^T W y} where \eqn{D} is the dummy-variable design
#' matrix implied by `design` and \eqn{W} is the diagonal weight matrix.
#'
#' @param design Either an integer category matrix of shape
#'   `(n_obs, n_factors)` with 1-based labels, a list of [Effect()] terms,
#'   or a reusable [Design()].
#' @param y Numeric vector of length `n_obs`.
#' @param options `NULL` for default [LsmrOptions()] or an options object.
#' @param weights Numeric vector of length `n_obs` or `NULL`.
#' @param preconditioner Controls preconditioning. The following input forms are
#'   accepted: `NULL` uses adaptive diagonal-to-Schwarz preconditioning,
#'   `PreconditionerConfig$Off` disables preconditioning,
#'   `PreconditionerConfig$Diagonal` uses diagonal/Jacobi scaling,
#'   [AdditiveSchwarz()] builds Schwarz immediately, [Adaptive()] tunes the
#'   default adaptive strategy, and a
#'   built [Preconditioner()] object reuses an existing factorization.
#' @return A named list with fields `x`, `unidentified`, `warnings`, `layout`,
#'   `demeaned`, `converged`, `iterations`, `residual`, `time_total`,
#'   `time_setup`, and `time_solve`.
#' @details Factor labels can have gaps. Coefficient addresses use the original
#'   1-based labels. Check `converged` before using `demeaned`; `warnings` contains
#'   diagnostic messages from native design and preconditioner construction.
#'   Native computations use one worker thread by default (adaptive builds may
#'   use a second). Set the environment
#'   variable `WITHINR_NUM_THREADS` before the first call to change this limit.
#' @examples
#' groups <- cbind(c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L))
#' fit <- within_solve(groups, c(1, 3, 4, 8))
#' fit$demeaned
#' @export
within_solve <- function(design,
                         y,
                         options = NULL,
                         weights = NULL,
                         preconditioner = NULL) {
  design <- validate_design(design)
  y <- as.double(y)
  weights <- validate_weights(weights)

  decorate_result(solve_impl(design, y, options, weights, preconditioner))
}

#' Solve fixed-effects normal equations for multiple response vectors
#'
#' Builds the operator and preconditioner once, then solves for each column
#' of `Y` in parallel.
#'
#' @inheritParams within_solve
#' @param Y Numeric matrix of shape `(n_obs, k)`.
#' @return A named list with matrix fields `x` and `demeaned`, plus
#'   `unidentified`, `warnings`, `layout`, `converged`, `iterations`,
#'   `residual`, `time_solve`, `time_setup`, and `time_total`.
#' @examples
#' groups <- cbind(c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L))
#' within_solve_batch(groups, cbind(c(1, 3, 4, 8), c(2, 5, 1, 7)))$demeaned
#' @export
within_solve_batch <- function(design,
                               Y,
                               options = NULL,
                               weights = NULL,
                               preconditioner = NULL) {
  design <- validate_design(design)
  if (!is.matrix(Y)) {
    stop("`Y` must be a matrix", call. = FALSE)
  }
  Y <- matrix(as.double(Y), nrow = nrow(Y), ncol = ncol(Y))
  weights <- validate_weights(weights)

  decorate_result(solve_batch_impl(design, Y, options, weights, preconditioner))
}

#' Persistent fixed-effects solver
#'
#' Retains the design and weight-dependent state across solves. Adaptive
#' preconditioning builds Schwarz only if progress stalls, then reuses it.
#'
#' @inheritParams within_solve
#' @return A `within_solver` object with `$solve()`, `$solve_batch()`, and
#'   `$preconditioner()` and `$has_escalated()` methods, plus `$n_dofs` and `$n_obs`
#'   fields. Native handles must be rebuilt after saving/loading an R session.
#'   A reused preconditioner must correspond to the same design and weights.
#' @export
Solver <- function(design, weights = NULL, preconditioner = NULL) {
  design <- validate_design(design)
  weights <- validate_weights(weights)
  ptr <- solver_new_impl(design, weights, preconditioner)

  solver <- new.env(parent = emptyenv())
  solver$ptr <- ptr
  solver$n_dofs <- solver_n_dofs_impl(ptr)
  solver$n_obs <- solver_n_obs_impl(ptr)
  solver$solve <- function(y, options = NULL) {
    decorate_result(solver_solve_impl(ptr, as.double(y), options))
  }
  solver$solve_batch <- function(Y, options = NULL) {
    if (!is.matrix(Y)) {
      stop("`Y` must be a matrix", call. = FALSE)
    }
    Y <- matrix(as.double(Y), nrow = nrow(Y), ncol = ncol(Y))
    decorate_result(solver_solve_batch_impl(ptr, Y, options))
  }
  solver$has_escalated <- function() solver_has_escalated_impl(ptr)
  solver$preconditioner <- function() {
    ptr <- solver_preconditioner_impl(solver$ptr)
    if (is.null(ptr)) {
      NULL
    } else {
      new_preconditioner(ptr)
    }
  }
  class(solver) <- "within_solver"
  solver
}

new_preconditioner <- function(ptr) {
  preconditioner <- new.env(parent = emptyenv())
  preconditioner$ptr <- ptr
  preconditioner$nrows <- preconditioner_nrows_impl(ptr)
  preconditioner$ncols <- preconditioner_ncols_impl(ptr)
  preconditioner$build_duration_seconds <- preconditioner_build_duration_impl(ptr)
  preconditioner$config <- preconditioner_config_impl(ptr)
  preconditioner$apply <- function(x) {
    preconditioner_apply_impl(ptr, as.double(x))
  }
  preconditioner$serialize <- function() {
    preconditioner_serialize_impl(ptr)
  }
  class(preconditioner) <- "within_preconditioner"
  preconditioner
}

#' Built preconditioner handle
#'
#' Deserializes a preconditioner from raw bytes produced by
#' `solver$preconditioner()$serialize()`. Bytes are specific to the native
#' library's wire format and are not a permanent cross-version storage format.
#'
#' @param data Raw vector containing serialized preconditioner bytes.
#' @return A `within_preconditioner` object with `$apply()` and `$serialize()`
#'   methods and `$nrows`, `$ncols`, `$config`, and
#'   `$build_duration_seconds` fields.
#' @export
Preconditioner <- function(data) {
  if (!is.raw(data)) {
    stop("`data` must be a raw vector", call. = FALSE)
  }
  new_preconditioner(preconditioner_deserialize_impl(data))
}

#' @export
print.within_solver <- function(x, ...) {
  cat(sprintf("<within_solver: n_obs=%d, n_dofs=%d>\n", x$n_obs, x$n_dofs))
  invisible(x)
}

# Matches the Python __repr__: Preconditioner(<variant>, n=<nrows>).
#' @export
print.within_preconditioner <- function(x, ...) {
  cat(sprintf(
    "Preconditioner(%s, n=%d)\n",
    preconditioner_variant_impl(x$ptr),
    x$nrows
  ))
  invisible(x)
}

#' @export
print.within_coefficient_layout <- function(x, ...) {
  cat(sprintf(
    "CoefficientLayout(n_terms=%d, n_dofs=%d)\n",
    x$n_terms,
    x$n_dofs
  ))
  invisible(x)
}

#' Reusable fixed-effects design
#'
#' Owns preprocessed observation storage independently of the R inputs. Share
#' a design across solvers with different weights without repeating preprocessing.
#' Rebuild this native handle after saving/loading an R session.
#' @inheritParams within_solve
#' @return A native design accepted wherever a category matrix is accepted.
#' @export
Design <- function(design) {
  structure(list(ptr = design_new_impl(validate_design(design))), class = "within_design")
}

#' Adaptive diagonal-to-Schwarz preconditioning
#'
#' Starts with diagonal scaling and builds additive Schwarz only when progress
#' stalls. This is the default strategy. A reused built preconditioner is fixed.
#' @inheritParams AdditiveSchwarz
#' @param window Positive number of consecutive stalled iterations.
#' @param threshold Contraction threshold in `[0, 1)`.
#' @return A preconditioner configuration object.
#' @export
Adaptive <- function(local_solver = NULL, reduction = ReductionStrategy$Auto,
                     window = 4L, threshold = 0.7) {
  structure(list(local_solver = local_solver, reduction = reduction,
                 window = window, threshold = threshold), class = "within_adaptive")
}
