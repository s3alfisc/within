# withinr

R bindings for the [within](https://github.com/py-econometrics/within) Rust
fixed-effects solver, using [extendr](https://extendr.rs/) and its R development
helper [rextendr](https://extendr.rs/contributing/rextendr-internals/).

## Generate, build, and check

Run these commands from the repository root. Prerequisites are R >= 4.2,
Rust/cargo >= 1.85, the usual R package compilation tools (Rtools on Windows),
and `tar` with xz support. Full CRAN checks also require LaTeX and, on Unix,
`checkbashisms` (Debian's `devscripts` or Homebrew's `checkbashisms`).
Use a native Rust toolchain; on Windows use the
`x86_64-pc-windows-gnu` target supported by Rtools.

```sh
Rscript -e 'install.packages(c("rextendr", "usethis", "roxygen2", "pkgload"), repos="https://cloud.r-project.org")'
Rscript scripts/prepare-r-package.R --vendor
R CMD build --no-build-vignettes withinr
R CMD check --as-cran withinr_0.1.0.tar.gz
R CMD INSTALL withinr_0.1.0.tar.gz
```

The generator requires **rextendr 0.5.0** to reproduce the build scaffolding.
If a later version is installed, install 0.5.0 from the CRAN archive using
`remotes::install_version("rextendr", version = "0.5.0")`. Roxygen generates
`NAMESPACE` and `man/`; its version is recorded in `DESCRIPTION`.
Development tools are not runtime dependencies of the installed package.

`--vendor` downloads locked registry dependencies at preparation time. The
resulting tarball installs with Cargo offline, without this repository, Git,
Python, rextendr, or roxygen2. To explicitly verify that Cargo cannot download:

```sh
CARGO_NET_OFFLINE=true R CMD check --as-cran --no-manual withinr_0.1.0.tar.gz
```

`--no-manual` is useful on machines without LaTeX; run the full check on a
machine with LaTeX before submission. The CI workflow runs native checks on
Linux, macOS, and Windows, including the declared Rust minimum on Linux.
CRAN acceptance still requires maintainer submission and review.

## What is generated

The maintained sources are `R/solve.R` (the R API) and `src/rust/src/` (typed
adapters to the native solver). Like Python's PyO3 annotations, `#[extendr]`
and `extendr_module!` generate the FFI trampolines and registration metadata.
There are no handwritten `.Call` wrappers or registration tables.

`scripts/prepare-r-package.R` does the following:

1. Copies the current workspace's Rust core into `src/rust/within-core/`,
   retaining its workspace metadata and omitting development dependencies.
   This uses the exact local sources, including unreleased changes.
2. Calls `rextendr::use_extendr()` in a temporary package and copies its
   generated entrypoint, build templates, configure/cleanup scripts, and Rust
   `document` binary. It adds locked dependency resolution and a two-job build
   limit to the generated build templates.
3. With `--vendor`, calls `rextendr::vendor_crates()`, removes precompiled
   WASI-only fixtures/objects and updates their per-file checksums, then bundles
   source in `vendor.tar.xz`. It records dependency authors and licenses in
   `inst/COPYRIGHTS`; original notices remain in the vendored sources.
4. Runs roxygen to bootstrap the namespace, then compiles the Rust adapter.
   The helper-generated `document` binary writes `R/extendr-wrappers.R` from
   extendr metadata, and roxygen finishes the namespace and reference pages.

The generated core snapshot, vendor archive, and compiled objects are ignored
by Git. Regenerate the snapshot after core edits and use `--vendor` after any
Rust dependency change. If dependencies change, first prepare the snapshot,
update the bridge lockfile, then vendor and generate:

```sh
Rscript scripts/prepare-r-package.R --no-document
cargo check --manifest-path withinr/src/rust/Cargo.toml
Rscript scripts/prepare-r-package.R --vendor
```

To demonstrate regeneration from scratch, delete only the generated files and
run the same command (CI also does this):

```sh
rm -f withinr/R/extendr-wrappers.R withinr/NAMESPACE
rm -f withinr/src/entrypoint.c withinr/src/rust/document.rs
rm -f withinr/src/Makevars.in withinr/src/Makevars.win.in withinr/src/withinr-win.def
rm -f withinr/configure withinr/configure.win withinr/cleanup withinr/cleanup.win
rm -rf withinr/tools withinr/man
Rscript scripts/prepare-r-package.R --vendor
```

## R API

Category matrices contain positive, 1-based integer labels. Gaps are allowed;
coefficient addresses retain the original labels. For character or factor
identifiers, explicitly encode each column with `as.integer(factor(ids))`.

```r
library(withinr)
groups <- cbind(c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L))
y <- c(1, 3, 4, 8)
x <- c(2, 5, 1, 7)

fit <- within_solve(groups, y)
fit$demeaned
batch <- within_solve_batch(groups, cbind(y, x), weights = c(1, 2, 1, 2))
stopifnot(all(batch$converged))
```

Results contain coefficients `x`, ordinary residuals `demeaned = y - D x`,
`converged`, `iterations`, `residual`, timings, `warnings`, unidentified
coefficient addresses, and a lazy coefficient `layout`. Batch coefficient and
residual matrices have one column per response. Always check convergence;
nonconvergence is reported in the result rather than silently treated as success.
As with other fixed-effects solvers, coefficients may not be uniquely identified.

Use `within_solve_batch(groups, cbind(y, X))$demeaned` to absorb fixed effects
before a regression; this is the integration surface for a future fixest adapter.
This package does not modify or integrate into fixest itself.

```r
design <- Design(groups)
solver <- Solver(design, weights = c(1, 2, 1, 2))
solver$solve(y, LsmrOptions(tol = 1e-10))
solver$solve_batch(cbind(y, x))
solver$has_escalated()

slopes <- list(Effect(groups[, 1], slopes = x), Effect(groups[, 2]))
within_solve(slopes, y)
```

`Design` owns its inputs and shares preprocessing across solvers with different
weights. `Solver` retains weight-dependent state for repeated solves.

`NULL` preconditioning uses the native **adaptive diagonal-to-Schwarz** default:
Schwarz is built only when progress stalls. Choose an explicit strategy with
`PreconditionerConfig$Off`, `$Diagonal`, `$Additive`, or `$Adaptive`.
`AdditiveSchwarz(local_solver, reduction)` and
`Adaptive(local_solver, reduction, window = 4L, threshold = 0.7)` allow tuning.
See `?LocalSolverConfig` for approximate Cholesky, Schur, scaling, and ridge options.

A solver's `$preconditioner()` returns its current built map (`NULL` for Off).
Before adaptive escalation this is a diagonal map. Reusing a built map fixes that
map; it does not recreate the adaptive strategy. Reuse only with the same design
and weights. For persistence use `p$serialize()` and `Preconditioner(bytes)`;
bytes are tied to the native wire-format version. Native design/solver/layout
handles must be rebuilt after saving and loading an R session.

Native computation uses a private one-thread pool; adaptive builds can create
a second pool, keeping combined default use within CRAN's two-core limit.
To request more threads,
set `Sys.setenv(WITHINR_NUM_THREADS = "4")` before the first native call.
Only Rust computation runs on workers; R allocation and error handling stay on
the R thread. This does not alter another package's Rayon pool.

## Local development checks

After preparing the package once:

```sh
Rscript -e 'pkgload::load_all("withinr"); source("withinr/tests/run_tests.R")'
cargo fmt --manifest-path withinr/src/rust/Cargo.toml -- --check
cargo clippy --locked --manifest-path withinr/src/rust/Cargo.toml -- -D warnings
```

The tests also run under `R CMD check`. They compare weighted and varying-slope
residuals with base R least squares, and exercise batches, persistent handles,
serialization, gappy labels, adaptive escalation, and invalid inputs.
