# within

`within` provides high-performance solvers for projecting out high-dimensional fixed effects from regression problems.

By the Frisch-Waugh-Lovell theorem, estimating a regression of the form *y = Xβ + Dα + ε* reduces to a sequence of least-squares projections, one for y and one for each column of X, followed by a cheap regression fit on the resulting residuals. The projection step of solving the normal equations *D'Dx = D'z* is the computational bottleneck, which is the problem `within` is designed to solve.

`within`'s solvers are tailored to the structure of fixed effects problems, which can be represented as a graph (as first noted by Correia, 2016). Concretely, `within` uses modified LSMR with a domain decomposition (Schwarz) preconditioner, backed by approximate Cholesky local solvers (Gao et al, 2025).

## Scope

`within` is a low-level fixed-effects kernel. Callers pass pre-factorized categorical codes: contiguous 0-based `uint32` level codes in F-order (column-major) arrays. Formula-level convenience — DataFrames, string/object categoricals, `pandas.factorize`, and formula parsing — is intentionally out of scope and belongs to a frontend layer built on top. The pyfixest-style workflow is served by such a frontend calling `within` underneath.

## Installation

You can install Python bindings from PyPi by running 

```bash
pip install within_py
```

Wheels ship for CPython 3.9+ (one `abi3` wheel per platform), free-threaded CPython 3.14+ (`cp314t`; PyO3 dropped experimental 3.13t support), and PyPy 3.11 (PyO3's minimum supported PyPy version; earlier PyPy releases are rejected at build time). CPython 3.15's `abi3t` (PEP 803, [PyO3/maturin#3064](https://github.com/PyO3/maturin/issues/3064)) will let the free-threaded build fold back into a single wheel per platform.

musllinux ships no PyPy wheel at all — the `musllinux_1_2` build image bundles no PyPy interpreter.

## Python Quickstart

`within`'s main user-facing function is `solve`. Provide a 2-D `uint32` array of category codes (one column per fixed-effect factor) and a response vector `y`. The solver finds x in the normal equations **D'D x = D'y**, where D is the sparse categorical design matrix.

```python
import numpy as np
from within import solve, solve_batch, LsmrOptions, PreconditionerConfig

np.random.seed(1)
n = 100_000
fe = np.asfortranarray(np.column_stack([
    np.random.randint(0, 500, n).astype(np.uint32),
    np.random.randint(0, 200, n).astype(np.uint32),
]))
y = np.random.randn(n)

# Default: additive Schwarz + LSMR
result = solve(fe, y)

# Custom tolerance / iteration cap
result = solve(fe, y, options=LsmrOptions(tol=1e-10, maxiter=2000))

# Weighted solve
result = solve(fe, y, weights=np.ones(n))

# Opt into diagonal/Jacobi preconditioning
result = solve(fe, y, preconditioner=PreconditionerConfig.Diagonal())
```

### Varying slopes

Pass a list of `Effect` terms instead of a categories array. Each term is a
factor's level codes plus an optional intercept and zero or more slope
covariates (per-level slopes, as in fixest's `f[z]` notation).

```python
from within import solve, Effect

firm = np.random.randint(0, 500, n).astype(np.uint32)
year = np.random.randint(0, 20, n).astype(np.uint32)
x = np.random.randn(n)

result = solve(
    [
        Effect(firm, intercept=True, slopes=[x]),
        Effect(year, intercept=True),
    ],
    y,
)

i = result.layout.index(0, 3, 1)
print(result.x[i])
```

## Python API

| Function | Description |
|---|---|
| `solve(design, y, weights?, options?, preconditioner?)` | Solve a single right-hand side. Returns `SolveResult`. |
| `solve_batch(design, Y, weights?, options?, preconditioner?)` | Solve multiple RHS vectors in parallel. `Y` has shape `(n_obs, k)`. Returns `BatchSolveResult`. |

`design` is either a 2-D `uint32` array or a list of `Effect` terms.

For repeated solves with the same design, `Solver` builds the preconditioner
once and reuses it.

```python
from within import Solver

solver = Solver(fe)
r = solver.solve(y)
r = solver.solve_batch(np.column_stack([y, X]))
precond = solver.preconditioner
solver2 = Solver(fe, preconditioner=precond)
```

| Property / Method | Description |
|---|---|
| `Solver(design, weights?, preconditioner?)` | Build a persistent solver. |
| `.solve(y, options?)` | Solve one right-hand side. |
| `.solve_batch(Y, options?)` | Solve multiple right-hand sides. |
| `.preconditioner` | Return the built, reusable `Preconditioner`, or `None`. |

## R quickstart

Requires R and a Rust toolchain (`cargo` on `PATH`).

From the repository root, use `devtools` to install R dependencies and build the
package:

```r
install.packages("devtools")
Sys.setenv(NOT_CRAN = "true")
devtools::install_deps("withinr/", dependencies = TRUE)
devtools::load_all("withinr/")
```

Example (FWL with two-way fixed effects):

```r
set.seed(42)
n <- 1000
n_firms <- 50L
n_years <- 20L

# 1-based fixed-effect ids in R
firm <- rep(seq_len(n_firms), each = n_years)
year <- rep(seq_len(n_years), times = n_firms)
categories <- cbind(firm, year)

beta <- 1.5
firm_fe <- rnorm(n_firms, sd = 3)[firm]
year_fe <- rnorm(n_years, sd = 1)[year]
x <- rnorm(n) + 0.3 * firm_fe
y <- beta * x + firm_fe + year_fe + rnorm(n, sd = 0.5)

res <- withinr::within_solve_batch(categories, cbind(y, x))
y_tilde <- res$demeaned[, 1]
x_tilde <- res$demeaned[, 2]
beta_hat <- sum(x_tilde * y_tilde) / sum(x_tilde^2)

print(beta_hat)
print(res$converged)
```

| Function | Description |
|---|---|
| `within_solve(design, y, options?, weights?, preconditioner?)` | Solve a single right-hand side. |
| `within_solve_batch(design, Y, options?, weights?, preconditioner?)` | Solve multiple RHS vectors in parallel. |

In R, `design` is either a 1-based integer category matrix or a list of
`Effect()` terms. For repeated solves, `Solver` builds the preconditioner once
and reuses it.

```r
solver <- withinr::Solver(categories)
r <- solver$solve(y)
r <- solver$solve_batch(cbind(y, x))

precond <- solver$preconditioner()
payload <- precond$serialize()
solver2 <- withinr::Solver(categories, preconditioner = withinr::Preconditioner(payload))
```

| Property / Method | Description |
|---|---|
| `Solver(design, weights?, preconditioner?)` | Build solver. Factorizes the preconditioner at construction. |
| `$solve(y, options?)` | Solve a single RHS with the given LSMR tuning. |
| `$solve_batch(Y, options?)` | Solve multiple RHS columns in parallel. |
| `$preconditioner()` | Return the built `Preconditioner`, or `NULL`. Reuse via `Solver(design, preconditioner=p)`. |

## Shared API concepts

### Solver configuration

| Class | Description |
|---|---|
| `LsmrOptions(tol=1e-8, maxiter=1000, local_size=None)` / `LsmrOptions(tol = 1e-8, maxiter = 1000L, local_size = NULL)` | Modified LSMR. `local_size` enables windowed reorthogonalization. |

### Preconditioners

The `preconditioner` argument accepts any of:

| Form | Meaning |
|---|---|
| `None` / `NULL` (default) | Library default — Additive Schwarz with sensible defaults. |
| `PreconditionerConfig.Off()` / `PreconditionerConfig$Off` | Explicit identity — solve unpreconditioned. |
| `PreconditionerConfig.Additive()` / `PreconditionerConfig$Additive` | Additive Schwarz shortcut, equivalent to the default. |
| `PreconditionerConfig.Diagonal()` / `PreconditionerConfig$Diagonal` | Diagonal/Jacobi preconditioner using `diag(D^T W D)^{-1}`. |
| `PreconditionerConfig.Additive(...)` / `AdditiveSchwarz(...)` | Tuned additive Schwarz configuration. |
| `Preconditioner` instance | Reuse a previously-built preconditioner across solvers. |

`PreconditionerConfig` is a tagged union: each variant is a subclass, so instances
support `match`/`case` and compare by value. A built preconditioner exposes
the configuration used to construct it as `.config` and its original build duration in seconds
as `.build_duration_seconds`; both are preserved across pickling and reuse.
R built-preconditioner environments expose the same information as `$config`
and `$build_duration_seconds` and can be serialized with `$serialize()`.

### Local solver configuration (advanced)

| Class | Description |
|---|---|
| `LocalSolverConfig(approx_chol?, schur?, dense_threshold?, scaling?, ridge?)` | Schur reduction + approximate Cholesky. Omit `schur` for the library default. |
| `Schur.approximate(config?)` / `Schur.exact()`; `Schur$Approximate(config?)` / `Schur$Exact` | Approximate or exact Schur reduction. |
| `ApproxCholConfig(seed=0, split_merge=None)` | Approximate Cholesky parameters. |
| `ApproxSchurConfig(seed=0, split=1)` | Approximate Schur complement sampling parameters. |
| `ScalingConfig(...)` | Signed-component dominance certification policy. |
| `ReductionStrategy` | `Auto` (default), `AtomicScatter`, `ParallelReduction` (class attributes, not an `Enum`). |

### Result types

**`SolveResult`**: `x` (coefficients), `unidentified` (directions the data cannot identify), `warnings`, `layout` (address-to-flat-index metadata; a lazy method-bearing handle in R), `demeaned` (residuals), `converged`, `iterations`, `residual`, `time_total`, `time_setup`, `time_solve`.

**`BatchSolveResult`**: Same fields, with `converged`, `iterations`, `residual`, and `time_solve` as lists (one entry per RHS).

Coefficients for unidentified directions are pinned to the **minimal-norm** value `0` (never NaN). This is why `x` can differ from reference tools that instead drop a reference level; the identified fit — `demeaned` — is unaffected by the choice.

## Rust API

```rust
use ndarray::Array2;
use within::{solve, LsmrOptions, PreconditionerConfig};
use within::config::{LocalSolverConfig, ReductionStrategy};

let categories = /* Array2<u32> of shape (n_obs, n_factors) */;
let y: &[f64] = /* response vector */;

// Default: LSMR + additive Schwarz (None → library default)
let r = solve(categories.view(), &y, None, &LsmrOptions::default(), None)?;
assert!(r.converged);

// Tighter tolerance with an explicit additive preconditioner
let lsmr = LsmrOptions { tol: 1e-10, ..LsmrOptions::default() };
let precond = PreconditionerConfig::Additive {
    local_solver: LocalSolverConfig::default(),
    reduction: ReductionStrategy::default(),
};
let r = solve(categories.view(), &y, None, &lsmr, &precond)?;

// Opt into diagonal/Jacobi preconditioning
let diagonal = PreconditionerConfig::Diagonal;
let r = solve(categories.view(), &y, None, &lsmr, &diagonal)?;
```

Persistent solver — build once, solve many:

```rust
use within::Solver;

let solver = Solver::new(categories.view(), None, None)?;
let r1 = solver.solve(&y, &LsmrOptions::default())?;
let r2 = solver.solve(&another_y, &LsmrOptions::default())?;  // reuses preconditioner
```

`solve` and `Solver::new` take the preconditioner as `impl Into<PreconditionerInput>`:
`None` (library default), a `&PreconditionerConfig` or owned `PreconditionerConfig`
(e.g. `PreconditionerConfig::Off` for the identity), or an owned/borrowed
`Preconditioner` for reuse. LSMR options are `impl Into<Option<&LsmrOptions>>`, so
`None` accepts the defaults and `&opts` overrides them.

| Type | Variants / Fields |
|---|---|
| `LsmrOptions` | `{ tol: f64, maxiter: usize, local_size: Option<usize> }` |
| `PreconditionerConfig` | `Off` \| `Additive { local_solver: LocalSolverConfig, reduction: ReductionStrategy }` \| `Diagonal` (`#[non_exhaustive]`) |
| `LocalSolverConfig` | `{ approx_chol, schur: SchurMode, dense_threshold, scaling, ridge }` |
| `SchurMode` | `Approximate(ApproxSchurConfig)` \| `Exact` |
| `Preconditioner` | Opaque built handle — reuse via `Solver::new(.., precond)` (owned or `&`) |
| `Effect` | `Effect::new(levels: &[u32], intercept: bool, slopes: impl IntoIterator<Item = &[f64]>) -> Result<Self, BuildError>` |
| `CoefficientAddress` | `{ channel: Channel { term, column }, level: usize }` |
| `CoefficientLayout` | `index(CoefficientAddress) -> Option<usize>`, `address(usize) -> Option<CoefficientAddress>`, `n_dofs()`, `n_terms()`, `n_levels(term)`, `n_columns(term)` |

### Varying slopes

Pass a `Vec<Effect>` instead of a categories array. Each term is a factor's
level codes plus an optional intercept and zero or more slope covariates
(per-level slopes, as in fixest's `f[z]` notation).

```rust
use within::{solve, Channel, CoefficientAddress, Effect};

let firm: &[u32] = /* level codes */;
let year: &[u32] = /* level codes */;
let x: &[f64] = /* covariate whose slope varies by firm */;

let terms = vec![
    Effect::new(firm, true, [x])?,  // firm intercept + firm-specific x slope
    Effect::new(year, true, [])?,   // year intercept
];
let r = solve(terms, &y, None, None, None)?;

// Read firm level 3's x-slope via the layout map (column 0 = intercept, 1 = first slope):
let at = CoefficientAddress { channel: Channel { term: 0, column: 1 }, level: 3 };
println!("{}", r.x[r.layout.index(at).unwrap()]);
```

`Solver::new` takes the same `Vec<Effect>`, so a slope design can be reused
across solves like any other.

### Lower-level access

| Module | Visibility | Key types |
|---|---|---|
| `within::config` | public | `LsmrOptions`, `PreconditionerConfig`, `LocalSolverConfig`, `SchurMode`, `ApproxCholConfig`, `ApproxSchurConfig`, `ScalingConfig`, `ReductionStrategy` |
| `within::observation` | public | `ObservationFrame` (columnar level-code + loading columns) |
| `within::error` | public | `WithinError`, `BuildError`, `SolveError` |
| `block_elim` / `channel` / `csr_block` / `domain` / `operator` / `solver` | `pub(crate)` | implementation layers — public items are re-exported at the crate root |

### Feature flags

| Feature | Default | Effect |
|---|---|---|
| `ndarray` | yes | Enables `from_array` constructors for `ndarray::ArrayView2` interop. |

## Project structure

```
crates/
  schwarz-precond/   Generic domain decomposition library (traits, solvers, Schwarz preconditioners)
  within/            Core fixed-effects solver (observation stores, domains, operators, orchestration)
  within-py/         PyO3 bridge (cdylib → within._within)
  within-r/          Workspace mirror for the extendr bridge
python/within/       Python package re-exporting the Rust extension
withinr/             R package pinned to within 0.4.0 and patched to the workspace in development
benchmarks/          Python benchmark framework
```

## Development

Uses [pixi](https://pixi.sh) as the task runner.

```bash
pixi run develop                         # Build Rust extension (release mode)
pixi run test                            # Rebuild + pytest
cargo test --workspace                   # Rust tests only
cargo bench -p within                    # Criterion benchmarks
pixi run python -m benchmarks run all    # Python benchmarks
```

Rust changes require rebuilding before running Python code (`pixi run develop`).

## License

MIT

## References

- Correia, Sergio. "A feasible estimator for linear models with multi-way fixed effects." *Preprint* at http://scorreia.com/research/hdfe.pdf (2016).
- Gao, Y., Kyng, R. & Spielman, D. A. (2025). AC(k): Robust Solution of Laplacian Equations by Randomized Approximate Cholesky Factorization. *SIAM Journal on Scientific Computing*.
- Toselli & Widlund (2005). *Domain Decomposition Methods — Algorithms and Theory*. Springer.
- Xu, J. (1992). Iterative Methods by Space Decomposition and Subspace Correction. *SIAM Review*, 34(4), 581--613.
