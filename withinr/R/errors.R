# Errors are conditions of class `within_error`, refined as
# `within_input_error` (invalid input: Python's ValueError/TypeError) or
# `within_runtime_error` (a failure on valid input: Python's RuntimeError).

abort <- function(..., kind = "input") {
  stop(errorCondition(
    paste0(...),
    class = c(paste0("within_", kind, "_error"), "within_error"),
    call = NULL
  ))
}

# Evaluate a native call, turning the bridge's tagged error messages (see
# `convert::guarded` in src/rust/src/convert.rs) into classed conditions.
native <- function(expr) {
  tryCatch(expr, error = function(e) {
    if (inherits(e, "within_error")) stop(e)
    message <- conditionMessage(e)
    for (kind in c("input", "runtime")) {
      tag <- sprintf("[within:%s] ", kind)
      if (startsWith(message, tag)) {
        abort(substring(message, nchar(tag) + 1L), kind = kind)
      }
    }
    stop(e)
  })
}
