# Designs: `Effect()` terms, persistent `Design()` handles, and the coercion
# of R-native inputs (factors, data frames, numeric matrices) to the integer
# labels the native library reads.

# Factor labels: factors and character vectors use their level codes (1..K);
# numeric labels are non-negative whole numbers, passed through unchanged.
as_labels <- function(x, name) {
  if (is.factor(x)) {
    return(as.integer(x))
  }
  if (is.character(x)) {
    return(as.integer(factor(x)))
  }
  if (is.integer(x)) {
    return(x) # NA and negative labels are rejected natively.
  }
  if (is.double(x)) {
    if (anyNA(x) || any(!is.finite(x)) || any(x != trunc(x)) ||
        any(x < 0) || any(x > .Machine$integer.max)) {
      stop(sprintf(
        "`%s` labels must be non-negative whole numbers without NA", name
      ), call. = FALSE)
    }
    return(as.integer(x))
  }
  stop(sprintf(
    "`%s` must hold factor labels (integer, numeric, factor, or character)", name
  ), call. = FALSE)
}

as_slopes <- function(slopes) {
  if (is.null(slopes)) {
    return(list())
  }
  if (is.data.frame(slopes)) {
    slopes <- as.list(slopes)
  } else if (is.matrix(slopes)) {
    slopes <- stats::setNames(
      lapply(seq_len(ncol(slopes)), function(j) slopes[, j]),
      colnames(slopes)
    )
  } else if (is.numeric(slopes)) {
    slopes <- list(slopes)
  } else if (!is.list(slopes)) {
    stop("`slopes` must be NULL, a numeric vector, a matrix, a data frame, or a list",
         call. = FALSE)
  }
  lapply(seq_along(slopes), function(j) {
    slope <- slopes[[j]]
    if (!is.numeric(slope)) stop(sprintf("slope %d must be numeric", j), call. = FALSE)
    slope <- as.double(slope)
    if (any(!is.finite(slope))) {
      stop(sprintf("slope %d must be finite (no NA, NaN, or Inf)", j), call. = FALSE)
    }
    slope
  }) |> stats::setNames(names(slopes))
}

#' One fixed-effect term, with optional varying slopes
#'
#' A term pairs a factor with the columns its levels carry: an optional
#' per-level intercept and zero or more per-level slopes on covariates. For
#' example, `Effect(firm, intercept = TRUE, slopes = tenure)` absorbs a firm
#' fixed effect and a firm-specific slope on `tenure`.
#'
#' @param levels Factor labels, one per observation: a factor or character
#'   vector (its level codes are used) or non-negative whole numbers (used as
#'   the labels). Labels need not be contiguous.
#' @param intercept Whether to include the per-level intercept.
#' @param slopes Slope covariates: `NULL`, a numeric vector, a numeric matrix
#'   or data frame (one slope per column), or a list of numeric vectors. Each
#'   slope has one finite value per observation.
#' @return A `within_effect` object. A list of these is a design, accepted by
#'   [within_solve()], [within_solve_batch()], [Design()], and [Solver()].
#' @examples
#' firm <- c(1, 1, 2, 2, 3, 3)
#' year <- c(1, 2, 1, 2, 1, 2)
#' tenure <- c(1, 2, 1, 3, 2, 5)
#' design <- list(Effect(firm, slopes = tenure), Effect(year))
#' within_solve(design, c(1, 2, 3, 5, 4, 8))$demeaned
#' @export
Effect <- function(levels, intercept = TRUE, slopes = NULL) {
  levels <- as_labels(levels, "levels")
  if (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept)) {
    stop("`intercept` must be TRUE or FALSE", call. = FALSE)
  }
  slopes <- as_slopes(slopes)
  if (any(lengths(slopes) != length(levels))) {
    stop("each slope must have one value per observation (the length of `levels`)",
         call. = FALSE)
  }
  if (!intercept && length(slopes) == 0L) {
    stop("an Effect needs an intercept or at least one slope", call. = FALSE)
  }
  structure(
    list(levels = levels, intercept = intercept, slopes = slopes),
    class = "within_effect"
  )
}

#' @export
print.within_effect <- function(x, ...) {
  cat(sprintf(
    "<within_effect: %d obs, %d levels, intercept=%s, %d slope(s)>\n",
    length(x$levels), length(unique(x$levels)), x$intercept, length(x$slopes)
  ))
  invisible(x)
}

# Coerce any accepted design form to what the native library reads: a
# `Design()` handle, an integer matrix, or a list of `Effect()`s.
as_design <- function(design) {
  if (inherits(design, "within_design")) {
    return(design)
  }
  if (inherits(design, "within_effect")) {
    return(list(design))
  }
  if (is.data.frame(design)) {
    columns <- lapply(seq_along(design), function(j) {
      as_labels(design[[j]], sprintf("design column %d", j))
    })
    return(matrix(unlist(columns, use.names = FALSE), ncol = length(columns)))
  }
  if (is.matrix(design)) {
    if (is.integer(design)) return(design)
    labels <- as_labels(as.vector(design), "design")
    return(matrix(labels, nrow = nrow(design), ncol = ncol(design)))
  }
  if (is.list(design)) {
    if (!length(design) || !all(vapply(design, inherits, logical(1), "within_effect"))) {
      stop("a list `design` must contain only Effect() terms", call. = FALSE)
    }
    return(design)
  }
  if (is.atomic(design) && is.null(dim(design))) {
    return(matrix(as_labels(design, "design"), ncol = 1L))
  }
  stop(
    "`design` must be a Design(), a matrix or data frame of factor labels, ",
    "or a list of Effect() terms",
    call. = FALSE
  )
}

#' Persistent fixed-effects design
#'
#' Builds and owns the native observation storage once, so several solvers
#' (for example with different weights) share the preprocessing. The handle
#' is independent of later changes to the R inputs.
#'
#' @param design A matrix or data frame of factor labels with one column per
#'   factor (factor and character columns use their level codes), a list of
#'   [Effect()] terms, or an existing `Design()`.
#' @return A `within_design` object with fields `n_obs` and `n_dofs`, accepted
#'   wherever a design is. Native handles do not survive `saveRDS()`; rebuild
#'   them in a new session.
#' @examples
#' design <- Design(data.frame(firm = c("a", "a", "b", "b"), year = c(1, 2, 1, 2)))
#' design$n_dofs
#' @export
Design <- function(design) {
  design <- as_design(design)
  if (inherits(design, "within_design")) {
    return(design)
  }
  structure(design_new_impl(design), class = "within_design")
}

#' @export
print.within_design <- function(x, ...) {
  cat(sprintf("<within_design: n_obs=%d, n_dofs=%d>\n", x$n_obs, x$n_dofs))
  invisible(x)
}
