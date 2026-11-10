# Configuration objects. Each constructor round-trips its arguments through the
# native library (`config_normalize_impl`), which applies the core's defaults
# and validation, so every object holds its fully resolved, effective
# configuration. Two equal configurations are `identical()`, as they compare
# equal in Python. The R side only coerces; it adds no rules of its own.

normalize_config <- function(fields, class) {
  native(config_normalize_impl(structure(fields, class = c(class, "within_config"))))
}

#' Modified LSMR solver options
#'
#' Uses Modified Golub-Kahan bidiagonalization to solve the rectangular
#' least-squares problem directly. The preconditioner \eqn{M \approx A^T A} is
#' applied as one \eqn{M^{-1}} solve per iteration.
#'
#' @param tol Convergence tolerance.
#' @param maxiter Maximum number of LSMR iterations.
#' @param local_size Number of past `v` vectors to reorthogonalize against via
#'   windowed modified Gram-Schmidt, or `NULL` (default) for the plain short
#'   recurrence. `5` to `20` is cheap insurance for ill-conditioned problems.
#' @details Count arguments here and in the other configuration objects are
#'   whole numbers no larger than `.Machine$integer.max`, so they round-trip
#'   through R integers.
#' @return A `within_lsmr_options` object, accepted as `options` by
#'   [within_solve()], [within_solve_batch()], and `Solver()$solve()`.
#' @examples
#' LsmrOptions(tol = 1e-10, maxiter = 5000)
#' @export
LsmrOptions <- function(tol = 1e-8, maxiter = 1000L, local_size = NULL) {
  normalize_config(
    list(tol = tol, maxiter = maxiter, local_size = local_size),
    "within_lsmr_options"
  )
}

#' Approximate Cholesky configuration
#'
#' @param seed Seed for the factorization's sampler, a whole number in
#'   `[0, 2^53]` (the doubles that represent integers exactly)
#'   (also used by the approximate Schur complement).
#' @param split_merge Split/merge parameter, or `NULL`.
#' @return A `within_approx_chol_config` object for [LocalSolverConfig()].
#' @export
ApproxCholConfig <- function(seed = 0, split_merge = NULL) {
  normalize_config(
    list(seed = seed, split_merge = split_merge),
    "within_approx_chol_config"
  )
}

#' Approximate Schur complement configuration
#'
#' Configures the clique-tree sampler; it draws from [ApproxCholConfig()]'s seed.
#'
#' @param split Positive edge split factor.
#' @return A `within_approx_schur_config` object for `Schur$approximate()`.
#' @export
ApproxSchurConfig <- function(split = 1L) {
  normalize_config(list(split = split), "within_approx_schur_config")
}

#' Schur-complement reduction modes
#'
#' `Schur$approximate(config = NULL)` samples the complement (the library
#' default) and `Schur$exact()` forms it exactly (higher fidelity, slower per
#' subdomain). Omitting `schur` in [LocalSolverConfig()] selects the default.
#'
#' @format A list of two constructor functions, `approximate` and `exact`.
#' @examples
#' Schur$approximate(ApproxSchurConfig(split = 2))
#' Schur$exact()
#' @export
Schur <- list(
  approximate = function(config = NULL) {
    normalize_config(list(mode = "approximate", config = config), "within_schur")
  },
  exact = function() {
    normalize_config(list(mode = "exact", config = NULL), "within_schur")
  }
)

#' Signed-component scaling certification
#'
#' @param tolerance Certification tolerance in `[0, 1)`, or `NULL` for the
#'   library default.
#' @param max_iterations Reduced-CG iteration budget, or `NULL`.
#' @param on_failure `"warn"` (clamp residual deficits, which only affects
#'   preconditioner quality, and emit a warning) or `"error"` (fail the
#'   build), or `NULL`.
#' @return A `within_scaling_config` object for [LocalSolverConfig()].
#' @export
ScalingConfig <- function(tolerance = NULL, max_iterations = NULL, on_failure = NULL) {
  normalize_config(
    list(tolerance = tolerance, max_iterations = max_iterations, on_failure = on_failure),
    "within_scaling_config"
  )
}

#' Local solver configuration for Schwarz subdomains
#'
#' Schur reduction plus approximate Cholesky. Arguments left `NULL` take the
#' library defaults; the returned object reports the effective values.
#'
#' @param approx_chol An [ApproxCholConfig()], or `NULL`.
#' @param schur `Schur$approximate(...)` or `Schur$exact()`, or `NULL`
#'   (approximate).
#' @param dense_threshold Non-negative dense Schur threshold, or `NULL`.
#' @param scaling A [ScalingConfig()], or `NULL`.
#' @return A `within_local_solver_config` object for
#'   `PreconditionerConfig$Additive()` and `PreconditionerConfig$Adaptive()`.
#' @examples
#' LocalSolverConfig(
#'   approx_chol = ApproxCholConfig(split_merge = 2),
#'   schur = Schur$approximate(ApproxSchurConfig(split = 2))
#' )
#' @export
LocalSolverConfig <- function(approx_chol = NULL, schur = NULL,
                              dense_threshold = NULL, scaling = NULL) {
  normalize_config(
    list(
      approx_chol = approx_chol, schur = schur,
      dense_threshold = dense_threshold, scaling = scaling
    ),
    "within_local_solver_config"
  )
}

#' Adaptive escalation rule
#'
#' Escalates after `window` consecutive contraction ratios exceed `threshold`.
#'
#' @param window Positive number of consecutive stalled iterations, or `NULL`
#'   (4).
#' @param threshold Contraction threshold in `[0, 1)`, or `NULL` (0.7).
#' @return A `within_staleness` object for `PreconditionerConfig$Adaptive()`.
#' @export
Staleness <- function(window = NULL, threshold = NULL) {
  normalize_config(list(window = window, threshold = threshold), "within_staleness")
}

#' Preconditioner configurations
#'
#' A tagged union mirroring Python's `within.PreconditionerConfig`:
#'
#' * `PreconditionerConfig$Off()`: no preconditioner.
#' * `PreconditionerConfig$Diagonal()`: diagonal (Jacobi) scaling by
#'   \eqn{diag(D^T W D)^{-1}}.
#' * `PreconditionerConfig$Additive(local_solver = NULL)`: additive Schwarz,
#'   built up front.
#' * `PreconditionerConfig$Adaptive(local_solver = NULL, stall = NULL)`: the
#'   library default. Starts on the diagonal and builds additive Schwarz only
#'   when the contraction stalls (see [Staleness()]), so a build error surfaces
#'   from the solve, not the constructor.
#'
#' Wherever a `preconditioner` is accepted, the strings `"off"`, `"diagonal"`,
#' `"additive"`, and `"adaptive"` are shorthands for these defaults. Each object
#' has a `variant` field (`"Off"`, `"Diagonal"`, `"Additive"`, `"Adaptive"`) and,
#' where applicable, `local_solver` and `stall` fields.
#'
#' @format A list of four constructor functions.
#' @examples
#' PreconditionerConfig$Adaptive(stall = Staleness(window = 8))
#' PreconditionerConfig$Additive(LocalSolverConfig(schur = Schur$exact()))
#' identical(PreconditionerConfig$Adaptive(), PreconditionerConfig$Adaptive())
#' @export
PreconditionerConfig <- list(
  Off = function() {
    normalize_config(list(variant = "Off"), "within_preconditioner_config")
  },
  Diagonal = function() {
    normalize_config(list(variant = "Diagonal"), "within_preconditioner_config")
  },
  Additive = function(local_solver = NULL) {
    normalize_config(
      list(variant = "Additive", local_solver = local_solver),
      "within_preconditioner_config"
    )
  },
  Adaptive = function(local_solver = NULL, stall = NULL) {
    normalize_config(
      list(variant = "Adaptive", local_solver = local_solver, stall = stall),
      "within_preconditioner_config"
    )
  }
)

as_preconditioner <- function(preconditioner) {
  if (is.character(preconditioner) && length(preconditioner) == 1L && !is.na(preconditioner)) {
    variants <- c(off = "Off", diagonal = "Diagonal", additive = "Additive",
                  adaptive = "Adaptive")
    variant <- variants[tolower(preconditioner)]
    if (is.na(variant)) {
      abort(sprintf(
        "`preconditioner` must be one of \"off\", \"diagonal\", \"additive\", or \"adaptive\", not \"%s\"",
        preconditioner
      ))
    }
    return(PreconditionerConfig[[variant]]())
  }
  preconditioner
}

as_options <- function(options) {
  if (is.list(options) && !inherits(options, "within_config")) {
    return(do.call(LsmrOptions, options))
  }
  options
}

config_call <- function(x) {
  switch(class(x)[[1L]],
    within_lsmr_options = "LsmrOptions",
    within_approx_chol_config = "ApproxCholConfig",
    within_approx_schur_config = "ApproxSchurConfig",
    within_scaling_config = "ScalingConfig",
    within_local_solver_config = "LocalSolverConfig",
    within_staleness = "Staleness",
    within_schur = paste0("Schur$", x$mode),
    within_preconditioner_config = paste0("PreconditionerConfig$", x$variant)
  )
}

#' @export
format.within_config <- function(x, ...) {
  fields <- unclass(x)
  fields$variant <- NULL
  fields$mode <- NULL
  if (inherits(x, "within_schur")) fields <- Filter(Negate(is.null), fields)
  args <- vapply(names(fields), function(name) {
    value <- fields[[name]]
    shown <- if (is.null(value)) {
      "NULL"
    } else if (inherits(value, "within_config")) {
      format(value)
    } else if (is.character(value)) {
      encodeString(value, quote = "\"")
    } else {
      format(value)
    }
    paste(name, "=", shown)
  }, character(1))
  paste0(config_call(x), "(", paste(args, collapse = ", "), ")")
}

#' @export
print.within_config <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}
