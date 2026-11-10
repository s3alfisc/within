# One-shot solves and their results.

as_response <- function(y, name = "y") {
  if (!is.numeric(y) || !is.null(dim(y))) {
    abort(sprintf("`%s` must be a numeric vector", name))
  }
  if (!is.double(y)) y <- as.double(y)
  y
}

as_response_matrix <- function(Y) {
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  if (!is.matrix(Y) || !is.numeric(Y)) {
    abort("`Y` must be a numeric matrix with one response per column")
  }
  if (!is.double(Y)) storage.mode(Y) <- "double"
  Y
}

as_weights <- function(weights) {
  if (is.null(weights)) NULL else as_response(weights, "weights")
}

emit_build_warnings <- function(messages) {
  for (message in messages) {
    warning(warningCondition(message, class = "within_build_warning", call = NULL))
  }
}

new_result <- function(raw, class, response_names = NULL) {
  emit_build_warnings(raw$warnings)
  raw$warnings <- NULL
  raw$unidentified <- as.data.frame(raw$unidentified)
  raw$layout <- new_coefficient_layout(raw$layout)
  if (!is.null(response_names)) {
    colnames(raw$x) <- response_names
    colnames(raw$demeaned) <- response_names
    for (field in c("converged", "iterations", "residual", "time_solve")) {
      names(raw[[field]]) <- response_names
    }
  }
  structure(raw, class = class)
}

solve_result <- function(raw) new_result(raw, "within_solve_result")

batch_result <- function(raw, Y) {
  new_result(raw, "within_batch_solve_result", colnames(Y))
}

#' Solve fixed-effects normal equations
#'
#' Computes fixed-effect coefficients by solving the normal equations
#' \eqn{D^T W D x = D^T W y}, where \eqn{D} is the design matrix implied by
#' `design` (dummies for each factor's intercept, and per-level slope columns
#' for varying slopes) and \eqn{W} is the diagonal weight matrix.
#'
#' @param design A matrix or data frame of factor labels with one column per
#'   factor, a list of [Effect()] terms (for varying slopes), or a persistent
#'   [Design()]. Factor and character columns use their level codes; numeric
#'   labels are non-negative whole numbers and need not be contiguous.
#' @param y Numeric response vector, one value per observation.
#' @param weights Observation weights, or `NULL` for unit weights.
#' @param options [LsmrOptions()], a list of its arguments, or `NULL` for the
#'   defaults (`tol = 1e-8`, `maxiter = 1000`).
#' @param preconditioner `NULL` (the default) is
#'   `PreconditionerConfig$Adaptive()`, which starts on the diagonal and
#'   escalates to additive Schwarz on a stalled contraction. Pass any
#'   [PreconditionerConfig] variant (or its shorthand string, such as
#'   `"diagonal"`) to choose another, or a built [Preconditioner()] to reuse an
#'   existing factorization.
#' @return A `within_solve_result`: a list with
#'   * `x`: coefficients, term-major (column `c` of the level at compact
#'     position `p` sits at `term_offset + (c - 1) * n_levels + p`). Use
#'     `layout` to translate. Unidentified slots hold `0`, never `NaN`.
#'   * `unidentified`: data frame (`term`, `level`, `column`) of per-level
#'     directions the data cannot identify.
#'   * `layout`: the coefficient layout, see below.
#'   * `demeaned`: the response minus the fitted fixed effects, `y - D x`.
#'   * `converged`, `iterations`, `residual` (relative normal-equation
#'     residual estimated by LSMR).
#'   * `time_total`, `time_setup`, `time_solve`: wall-clock seconds.
#'
#'   `coef()` returns `x` and `residuals()` returns `demeaned`.
#' @section Coefficient layout:
#' `layout` mirrors Python's `CoefficientLayout`. Terms and columns are
#' 1-based (columns are the intercept, if present, then the slopes in order);
#' levels are your own labels (factor codes for factor input). Methods:
#' `n_dofs()`, `n_terms()`, `n_levels(term)`, `n_columns(term)`,
#' `index(term, level, column)` (vectorized; the position in `x`), and
#' `address(index)` (a data frame of `term`, `level`, `column`).
#' @section Warnings and threads:
#' Non-fatal preconditioner build diagnostics are signalled as warnings of
#' class `within_build_warning`. Native work runs on a private pool of one
#' thread (an adaptive build may add a second); set the environment variable
#' `WITHINR_NUM_THREADS` before the first solve to use more.
#' @examples
#' groups <- cbind(firm = c(1, 1, 2, 2), year = c(1, 2, 1, 2))
#' fit <- within_solve(groups, c(1, 3, 4, 8))
#' fit$demeaned
#' fit$layout$index(term = 1, level = 2, column = 1)
#' @export
within_solve <- function(design, y, weights = NULL, options = NULL, preconditioner = NULL) {
  solve_result(native(solve_impl(
    as_design(design), as_response(y), as_weights(weights),
    as_options(options), as_preconditioner(preconditioner)
  )))
}

#' Solve fixed-effects normal equations for several responses
#'
#' Equivalent to [within_solve()] on each column of `Y`, but builds the design
#' and preconditioner once and solves the columns in parallel on the native
#' worker pool (one thread unless `WITHINR_NUM_THREADS` says otherwise; see
#' [within_solve()]).
#'
#' @inheritParams within_solve
#' @param Y Numeric matrix (or data frame) with one response per column.
#' @return A `within_batch_solve_result`: like [within_solve()]'s result, but
#'   `x` and `demeaned` are matrices with one column per response, and
#'   `converged`, `iterations`, `residual`, and `time_solve` are vectors.
#'   `time_setup` is the shared setup time. Column names of `Y` carry over.
#' @examples
#' groups <- cbind(c(1, 1, 2, 2), c(1, 2, 1, 2))
#' within_solve_batch(groups, cbind(y = c(1, 3, 4, 8), x = c(2, 5, 1, 7)))$demeaned
#' @export
within_solve_batch <- function(design, Y, weights = NULL, options = NULL,
                               preconditioner = NULL) {
  Y <- as_response_matrix(Y)
  batch_result(native(solve_batch_impl(
    as_design(design), Y, as_weights(weights),
    as_options(options), as_preconditioner(preconditioner)
  )), Y)
}

#' @export
print.within_solve_result <- function(x, ...) {
  cat("<within_solve_result>\n")
  cat(sprintf(
    "  converged: %s after %d iterations (residual %.3g)\n",
    x$converged, x$iterations, x$residual
  ))
  cat(sprintf(
    "  %d observations, %d coefficients in %d term(s), %d unidentified\n",
    length(x$demeaned), x$layout$n_dofs(), x$layout$n_terms(), nrow(x$unidentified)
  ))
  cat(sprintf(
    "  time: %.4gs (setup %.4gs, solve %.4gs)\n",
    x$time_total, x$time_setup, x$time_solve
  ))
  invisible(x)
}

#' @export
print.within_batch_solve_result <- function(x, ...) {
  cat("<within_batch_solve_result>\n")
  iterations <- if (length(x$iterations)) {
    paste0("; iterations ", paste(unique(range(x$iterations)), collapse = "-"))
  } else {
    ""
  }
  cat(sprintf(
    "  %d responses, %d converged%s\n",
    length(x$converged), sum(x$converged), iterations
  ))
  cat(sprintf(
    "  %d observations, %d coefficients in %d term(s), %d unidentified\n",
    nrow(x$demeaned), x$layout$n_dofs(), x$layout$n_terms(), nrow(x$unidentified)
  ))
  cat(sprintf("  time: %.4gs (setup %.4gs)\n", x$time_total, x$time_setup))
  invisible(x)
}

#' @export
coef.within_solve_result <- function(object, ...) object$x

#' @export
coef.within_batch_solve_result <- function(object, ...) object$x

#' @export
residuals.within_solve_result <- function(object, ...) object$demeaned

#' @export
residuals.within_batch_solve_result <- function(object, ...) object$demeaned

# ---------------------------------------------------------------------------
# Coefficient layout
# ---------------------------------------------------------------------------

new_coefficient_layout <- function(ptr) {
  info <- native(layout_info_impl(ptr))
  check_term <- function(term) {
    if (length(term) != 1L || is.na(term) || term < 1 || term > info$n_terms ||
        term != trunc(term)) {
      abort(sprintf("term %s out of range (n_terms=%d)", format(term), info$n_terms))
    }
    as.integer(term)
  }
  as_int <- function(v, name) {
    if (!is.numeric(v) || anyNA(v) || any(v != trunc(v)) || any(abs(v) > .Machine$integer.max)) {
      abort(sprintf("`%s` must hold whole numbers", name))
    }
    as.integer(v)
  }
  layout <- new.env(parent = emptyenv())
  layout$ptr <- ptr
  layout$n_dofs <- function() info$n_dofs
  layout$n_terms <- function() info$n_terms
  layout$n_levels <- function(term) info$n_levels[[check_term(term)]]
  layout$n_columns <- function(term) info$n_columns[[check_term(term)]]
  layout$index <- function(term, level, column) {
    lengths <- c(length(term), length(level), length(column))
    if (any(lengths == 0L)) return(integer())
    n <- max(lengths)
    recycle <- function(v, name) rep_len(as_int(v, name), n)
    native(layout_index_impl(
      ptr, recycle(term, "term"), recycle(level, "level"), recycle(column, "column")
    ))
  }
  layout$address <- function(index) {
    as.data.frame(native(layout_address_impl(ptr, as_int(index, "index"))))
  }
  lockEnvironment(layout, bindings = TRUE)
  structure(layout, class = "within_coefficient_layout")
}

#' @export
print.within_coefficient_layout <- function(x, ...) {
  cat(sprintf("CoefficientLayout(n_terms=%d, n_dofs=%d)\n", x$n_terms(), x$n_dofs()))
  invisible(x)
}
