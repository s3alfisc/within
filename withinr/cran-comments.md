## Local validation

Tested on macOS arm64, R 4.5.3, Rust/cargo 1.93.1:

`CARGO_NET_OFFLINE=true R CMD check --as-cran --no-manual withinr_0.1.0.tar.gz`

0 errors, 0 warnings, 2 notes:

* Unable to verify current time in the network-restricted test environment.
* `xcrun_db` left in the test temporary directory by the macOS SDK tooling.

Installation, loading/unloading, documentation, examples, and all tests passed.
The numerical tests compare weighted and varying-slope residuals with base R
least squares. Rust formatting and clippy checks also passed.

The repository CI is configured for Linux, macOS, and Windows, plus Rust 1.85
on Linux. Those remote jobs have not been run as part of this local validation.
Before submission, run them and the full manual check on a machine with LaTeX,
then update these comments with those results.

## Bundled Rust sources

The source tarball is approximately 4.1 MB because it includes locked Rust
dependencies in `src/rust/vendor.tar.xz`. Compilation is offline; no code is
downloaded during installation. The within core is included from the local
workspace, so installation does not require an unpublished crates.io release.
Precompiled WASI-only fixtures and objects are omitted from the source bundle.
Authorship and licenses are recorded in `inst/COPYRIGHTS`, with original notices
preserved in the vendored sources.

The Rust/R FFI, routine registration, C entrypoint, and R wrappers are generated
by extendr/rextendr. Compilation uses at most two jobs. The default native worker
pool has one thread; a deferred adaptive build may create a second worker.
