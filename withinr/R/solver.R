# Persistent solvers and built preconditioners. Both are locked environments
# holding a native handle, mirroring the Python classes' methods and
# read-only properties.

#' Persistent solver with a cached preconditioner
#'
#' Build once for a design and weights, then call `$solve()` or
#' `$solve_batch()` repeatedly with different responses. The preconditioner is
#' built once: at construction or, under `PreconditionerConfig$Adaptive()`, at
#' the first escalation.
#'
#' @inheritParams within_solve
#' @return A `within_solver` object with
#'   * `$solve(y, options = NULL)` and `$solve_batch(Y, options = NULL)`, which
#'     return the same results as [within_solve()] and [within_solve_batch()];
#'   * `$preconditioner`: the [Preconditioner()] in use (`NULL` when
#'     unpreconditioned), for reuse or serialization. Under `Adaptive` it is
#'     the Schwarz map once built and, before that, the diagonal base carrying
#'     the ladder, which escalates when reused;
#'   * `$has_escalated`: whether an `Adaptive` solve has handed off to Schwarz;
#'   * `$n_dofs` and `$n_obs`.
#'
#'   A reused preconditioner must come from the same design and weights.
#'   Native handles do not survive `saveRDS()`; rebuild them in a new session.
#' @examples
#' groups <- cbind(c(1, 1, 2, 2, 3, 3), c(1, 2, 1, 2, 1, 2))
#' solver <- Solver(groups, weights = c(1, 2, 1, 2, 1, 2))
#' solver$solve(c(1, 3, 4, 8, 2, 5))$converged
#' solver$solve_batch(cbind(c(1, 3, 4, 8, 2, 5), c(2, 5, 1, 7, 3, 3)))$iterations
#' solver$has_escalated
#' @export
Solver <- function(design, weights = NULL, preconditioner = NULL) {
  built <- native(solver_new_impl(
    as_design(design), as_weights(weights), as_preconditioner(preconditioner)
  ))
  emit_build_warnings(built$warnings)
  ptr <- built$ptr

  solver <- new.env(parent = emptyenv())
  solver$ptr <- ptr
  solver$n_dofs <- built$n_dofs
  solver$n_obs <- built$n_obs
  solver$solve <- function(y, options = NULL) {
    solve_result(native(solver_solve_impl(ptr, as_response(y), as_options(options))))
  }
  solver$solve_batch <- function(Y, options = NULL) {
    Y <- as_response_matrix(Y)
    batch_result(native(solver_solve_batch_impl(ptr, Y, as_options(options))), Y)
  }
  makeActiveBinding("preconditioner", function() {
    handle <- native(solver_preconditioner_impl(ptr))
    if (is.null(handle)) NULL else new_preconditioner(handle)
  }, solver)
  makeActiveBinding("has_escalated", function() native(solver_has_escalated_impl(ptr)), solver)
  lockEnvironment(solver, bindings = TRUE)
  structure(solver, class = "within_solver")
}

#' @export
print.within_solver <- function(x, ...) {
  cat(sprintf("<within_solver: n_obs=%d, n_dofs=%d>\n", x$n_obs, x$n_dofs))
  invisible(x)
}

new_preconditioner <- function(ptr) {
  info <- native(preconditioner_info_impl(ptr))
  preconditioner <- new.env(parent = emptyenv())
  preconditioner$ptr <- ptr
  preconditioner$nrows <- info$nrows
  preconditioner$ncols <- info$ncols
  preconditioner$variant <- info$variant
  preconditioner$config <- info$config
  preconditioner$build_duration_seconds <- info$build_duration_seconds
  preconditioner$apply <- function(x) native(preconditioner_apply_impl(ptr, as_response(x, "x")))
  preconditioner$serialize <- function() native(preconditioner_serialize_impl(ptr))
  lockEnvironment(preconditioner, bindings = TRUE)
  structure(preconditioner, class = "within_preconditioner")
}

#' Built fixed-effects preconditioner
#'
#' A preconditioner is built once per design (see [Solver()]) and can be
#' reused across solvers and sessions. `Preconditioner(data)` restores one from
#' the raw bytes written by `$serialize()`, the R counterpart of pickling in
#' Python; the bytes are tied to the native wire-format version.
#'
#' @param data A raw vector produced by `$serialize()`.
#' @return A `within_preconditioner` object with
#'   * `$apply(x)`: \eqn{y = M^{-1} x};
#'   * `$serialize()`: the raw bytes;
#'   * `$nrows`, `$ncols`, `$variant` (`"Diagonal"`, `"Additive"`, or
#'     `"Adaptive"`), `$config` (the normalized [PreconditionerConfig] it was
#'     built from), and `$build_duration_seconds`, preserved across
#'     serialization.
#'
#'   Pass it as the `preconditioner` argument to skip the factorization.
#' @examples
#' groups <- cbind(c(1, 1, 2, 2, 3, 3), c(1, 2, 1, 2, 1, 2))
#' p <- Solver(groups, preconditioner = "additive")$preconditioner
#' bytes <- p$serialize()
#' restored <- Preconditioner(bytes)
#' within_solve(groups, c(1, 3, 4, 8, 2, 5), preconditioner = restored)$converged
#' @export
Preconditioner <- function(data) {
  if (!is.raw(data)) {
    abort("`data` must be a raw vector produced by `$serialize()`")
  }
  new_preconditioner(native(preconditioner_deserialize_impl(data)))
}

#' @export
print.within_preconditioner <- function(x, ...) {
  cat(sprintf("Preconditioner(%s, n=%d)\n", x$variant, x$nrows))
  invisible(x)
}
