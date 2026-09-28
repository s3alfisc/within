# xhdfe / within MLSMR comparison and attribution assessment

Date of review: 18 August 2026

## Executive assessment

My conclusion is:

1. **There is strong evidence that xhdfe implements a standalone C++ adaptation of the
   `within` MLSMR--additive-Schwarz solver architecture.** The evidence is much stronger
   than the fact that both projects happen to use LSMR. The implementations share an
   unusually specific collection of domain construction rules, partition weights, local
   graph solvers, numerical defaults, and a paired modified-Golub--Kahan recurrence. xhdfe
   also exposes `within`, `within-additive`, and `within_additive` as literal aliases for its
   `mlsmr` method.

2. **I did not find evidence that xhdfe imports, links to, or calls the Rust `within`
   package.** Its implementation is self-contained C++/Eigen code with substantial
   xhdfe-specific integration. I also cannot establish byte-for-byte or line-by-line source
   copying across the two languages from the repositories alone. Accordingly, the safest
   wording is not simply "xhdfe uses our solver," which can imply a software dependency.
   The precise claim is that xhdfe **implements/adapts the solver architecture developed in
   `within`**.

3. **The paper provides meaningful credit, but the provenance wording is too weak in light
   of the code.** It calls the Fischer--Schröder solver "closely related" and says shared
   updates were "helpful." Those formulations read as parallel work or inspiration. The
   source evidence is more consistent with a direct architectural adaptation. The paper
   should say so explicitly and distinguish the inherited solver design from xhdfe's own
   C++ implementation, selection logic, threading, interface work, fallbacks, and GPU sweep
   kernels.

4. **The paper's claim that MLSMR is not run on the GPU is accurate.** xhdfe's MLSMR
   implementation is CPU-only. GPU requests are routed to separate sweep/PCG paths; a
   larger workflow may use CUDA in another phase and CPU MLSMR in a repair or recovery
   phase, but there is no CUDA or Metal MLSMR kernel.

5. **A collegial request for revision is warranted.** I would ask for (a) "adapted from" or
   "based on" language in the technical section, (b) a versioned software citation, and
   (c) disclosure of the actual absorber used in each benchmark row. I would not allege
   code copying or a license violation on the present evidence.

## Scope and snapshots reviewed

- `within`: current checkout `e5b54dc`, plus public history back to March 2026.
- `xhdfe-xfe`: current checkout `fdf8116`; the paper identifies xhdfe 2.23.1 at
  `d9958df` (6 August 2026); the first public repository commit is the squashed root
  `763c724` (3 July 2026).
- Paper: `xhdfe: Fast high-dimensional fixed-effects estimation in Stata`, pages 5--33.
  Despite the `11JUN2026` filename, the title page and PDF metadata identify an August
  2026 revision, and its appendix identifies the 6 August xhdfe commit.

The complete current `src/fe_absorption.cpp` is byte-identical to the paper's `d9958df`
snapshot; only unrelated covariance-scaling code changed elsewhere afterward. This makes
the current line links below directly representative of the code used for the paper. At the
paper snapshot, approximately 91% of the MLSMR core lines blame to the 3 July root commit;
later changes mainly added deterministic OpenMP execution, concurrent right-hand-side
solves, and numerical hardening.

The 17 March disclosure date and the 28 April pull-request opening date are supplied by the
`within` authors. I did not inspect the underlying email/chat or GitHub event record. Those
records should be preserved because they would be the best evidence of what was actually
shared and when.

## What is distinctive about the `within` solver

The distinctive object is not bare LSMR. Its individual ingredients have prior literature.
The relevant contribution for this comparison is their integration into a fixed-effects
solver:

1. A rectangular, matrix-free least-squares operator, with
   \(A=\sqrt{W}D\) and \(b=\sqrt{W}y\), solved by modified LSMR.
2. A one-level additive-Schwarz preconditioner built from every cross-factor pair.
3. One subdomain for each connected component of each factor-pair bipartite graph.
4. Symmetric two-sided partition-of-unity weights \(1/\sqrt{c_j}\) for a degree of
   freedom covered by \(c_j\) subdomains.
5. A bipartite sign flip to obtain a Laplacian/SDDM local system.
6. Elimination of the larger factor side, leaving the smaller Schur system.
7. Exact dense treatment for small retained systems and sampled clique-tree Schur
   construction otherwise.
8. Approximate-Cholesky factorization of the reduced graph.
9. A modified-Golub--Kahan MLSMR recurrence with one preconditioner application per
   iteration and paired, M-weighted reorthogonalization of \(v\) and \(\widetilde p\).

The public `within` overview states the architecture directly in
[docs/2_solver_architecture.md](/Users/admin/Documents/within/docs/2_solver_architecture.md:14).
The high-level solve constructs the rectangular operator and calls MLSMR in
[solver.rs](/Users/admin/Documents/within/crates/within/src/solver.rs:448), while the
preconditioned solver's one-\(M^{-1}\)-application contract is stated in
[lsmr.rs](/Users/admin/Documents/within/crates/schwarz-precond/src/lsmr.rs:205).

## Source-level implementation match

The following choices occur together in both projects:

| Choice | `within` evidence | xhdfe evidence |
|---|---|---|
| Matrix-free rectangular \(\sqrt W D\) operator | [solver.rs](/Users/admin/Documents/within/crates/within/src/solver.rs:448) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:13690) |
| All unordered cross-factor pairs | [factor_pairs.rs](/Users/admin/Documents/within/crates/within/src/domain/factor_pairs.rs:35) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12722) |
| Pair graph split into connected-component domains | [factor_pairs.rs](/Users/admin/Documents/within/crates/within/src/domain/factor_pairs.rs:94) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12746) |
| Two-sided \(1/\sqrt{\text{cover count}}\) partition weights | [factor_pairs.rs](/Users/admin/Documents/within/crates/within/src/domain/factor_pairs.rs:154) | weights computed at [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12871) and applied on restriction and prolongation at [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:11956) |
| Bipartite sign flip / Laplacian coordinates | [sddm.rs](/Users/admin/Documents/within/crates/within/src/domain/factor_pairs/sddm.rs:54) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:11739) |
| Eliminate the larger factor side | [sddm.rs](/Users/admin/Documents/within/crates/within/src/domain/factor_pairs/sddm.rs:88) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12143) |
| Dense retained-Schur threshold 24 | [config.rs](/Users/admin/Documents/within/crates/within/src/config.rs:17) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12714) |
| Sampled clique-tree Schur, default seed 0 and split 1 | [config.rs](/Users/admin/Documents/within/crates/within/src/config.rs:118) and [schur.rs](/Users/admin/Documents/within/crates/within/src/block_elim/schur.rs:303) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12306) |
| Approximate Cholesky, default seed 0 and split/merge 2 | [config.rs](/Users/admin/Documents/within/crates/within/src/config.rs:74) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12485) |
| Paired \(v,\widetilde p\) M-weighted reorthogonalization | [bidiag.rs](/Users/admin/Documents/within/crates/schwarz-precond/src/lsmr/bidiag.rs:192) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:14050) |
| Modified-GK update with one preconditioner application per iteration | [bidiag.rs](/Users/admin/Documents/within/crates/schwarz-precond/src/lsmr/bidiag.rs:275) | [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:14135) |

Some of these choices can follow from the literature independently. The evidentiary weight
comes from their **combination**, the matching defaults (24; Schur seed/split 0/1;
ApproxChol seed/split-merge 0/2), the direct naming evidence below, and the chronology.
Bare LSMR or generic additive Schwarz does not prescribe this bundle.

### Direct naming and self-description

- xhdfe's Python parser maps `within`, `within-additive`, and `within_additive` directly to
  `AbsorptionMethod::Mlsmr`:
  [py_hdfe_v11.cpp](/Users/admin/Documents/xhdfe-xfe/python/py_hdfe_v11.cpp:168).
  The Stata help exposes the same aliases and calls the family CPU-only:
  [xhdfe.sthlp](/Users/admin/Documents/xhdfe-xfe/stata/xhdfe.sthlp:631).
- The R showcase calls `mlsmr` the "modified/within LSMR":
  [xhdfe_r_showcase.qmd](/Users/admin/Documents/xhdfe-xfe/r/examples/xhdfe_r_showcase.qmd:437).
- The paper itself says that `mlsmr` "is also available as `within`" (page 18).
- xhdfe's source labels the block an "LSMR / MLSMR additive-Schwarz Krylov absorber"
  and says it was "ported from Tiago's ... branch, 24jun2026":
  [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:10662).
  This establishes xhdfe's own description of the code as a port from an internal branch;
  it does not by itself identify that branch's upstream source.
- The automatic selector is separately labelled "ported from Tiago branch 24jun2026":
  [hdfe_regressor_v11.cpp](/Users/admin/Documents/xhdfe-xfe/src/hdfe_regressor_v11.cpp:2637).

The aliases and full solver were already present in xhdfe's first public, squashed commit
`763c724`. The README acknowledgment of `within` was added five days later, on 8 July, in
commit `dc29273` (`docs: acknowledge within (Fischer/Schröder) ...`). The current
acknowledgment is at [README.md](/Users/admin/Documents/xhdfe-xfe/README.md:583).

## Chronology

| Date | Event | Evidentiary status |
|---|---|---|
| 4 March 2026 | `within` commit `343fc9e` adds the Rust crates, already including factor-pair connected-component Schwarz domains, Schur reduction, and approximate Cholesky. | Public git history |
| 9 March 2026 | `within` commit `fd23c16` corrects the symmetric PoU weight to \(1/\sqrt c\). | Public git history |
| 14--17 March 2026 | Public documentation/release describes the graph-preconditioner architecture; annotated tag `v0.1.0` was created 17 March. | Public git history |
| 17 March 2026 | The algorithm was shared with the xhdfe authors. | Author-supplied; preserve the underlying communication |
| 23 April 2026 | `within` commits `3da8e9c` and `e150a80` add rectangular MLSMR and connect it to the existing Schwarz preconditioner. | Public git history |
| 28 April 2026 | MLSMR PR opened; commit `430df0a` adds windowed reorthogonalization. | PR opening date author-supplied; commit public |
| 29 April 2026 | Merge commit `ef3ad04` merges PR #25, "Add LSMR/MLSMR Krylov solver alongside CG." | Public git history |
| 4 June 2026 | `within` release `v0.2.0` makes modified LSMR the principal iterative path. | Public annotated tag |
| 24 June 2026 | Date in xhdfe's source comment for the internal branch from which its solver was ported. | xhdfe self-report; private branch unavailable |
| 3 July 2026 | xhdfe's first public squashed commit `763c724` already contains the complete matching stack and `within` aliases. | Public git history |
| 8 July 2026 | xhdfe commit `dc29273` adds the README acknowledgment of Fischer, Schröder, and `within`. | Public git history |
| 6 August 2026 | xhdfe paper/release snapshot `d9958df`, version 2.23.1. | Paper appendix and public git history |

This sequence supports architectural derivation but is not, by itself, proof of source-code
copying. The private 17 March communication is important: it can establish whether the
specific MLSMR recurrence and local-solver details were transmitted before the public PR.

## Material xhdfe contributions and differences

xhdfe is not merely a wrapper around `within`. It includes substantial independent
engineering:

- A self-contained C++/Eigen implementation integrated with Stata, Python, and R rather
  than a Rust dependency or FFI call.
- A probe-based automatic selector that promotes large, slowly contracting, eligible CPU
  designs to MLSMR and otherwise selects a sweep method:
  [hdfe_regressor_v11.cpp](/Users/admin/Documents/xhdfe-xfe/src/hdfe_regressor_v11.cpp:2849).
- CPU-only MLSMR restricted to standard categorical fixed effects; heterogeneous slopes,
  FE recovery, and group/individual models use other paths.
- xhdfe-specific heuristics: pair skipping when the retained side is at most 64 levels,
  ridge regularization at `1e-10`, exact-star threshold 16, a full-dense threshold 192,
  and dense/sparse/incomplete-Cholesky/diagonal/identity fallbacks:
  [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:12708).
- Its own normalization/column-scaling policy and invalid-preconditioner fallback behavior:
  [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:13665).
- A default iteration cap of 100,000 rather than `within`'s current 1,000.
- OpenMP batching across right-hand sides and xhdfe-specific convergence, certification,
  cache, and interface machinery:
  [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:14291).
- Independent compiled sweep kernels, CPU threading, CUDA/Metal dispatch, and GPU
  certification code. These are not part of the `within` MLSMR implementation.

Conversely, current `within` has generic Rust operator/preconditioner APIs, a persistent
reusable solver, varying-slope/SDDM handling, and later warm-start/escalation machinery that
should not be projected backward onto xhdfe's July implementation.

### Test-coverage caveat

The explicit xhdfe MLSMR tests I found use factor dimensions 60 by 10 in Stata
([methods.do](/Users/admin/Documents/xhdfe-xfe/tests/stata/part1/methods.do:3)) and 300 by
40 in the R/Python parity fixture
([helper-parity.R](/Users/admin/Documents/xhdfe-xfe/r/xhdfe/tests/testthat/helper-parity.R:10)).
Because xhdfe skips a factor pair whenever its smaller side is at most 64, both tests bypass
the pairwise Schwarz domains and exercise the diagonal fallback. I found no focused
automated test that forces the sampled clique-tree/approximate-Cholesky MLSMR path, the
default auto-selector's MLSMR promotion, or GPU rejection specifically for MLSMR. This is
a robustness/coverage concern, not evidence about provenance.

## GPU question

The paper is accurate that xhdfe does not run MLSMR on the GPU:

- Low-level dispatch intercepts MLSMR before GPU dispatch, calls
  `absorb_fixed_effects_mlsmr` into a variable named `cpu_result`, and marks a requested
  GPU backend unavailable:
  [fe_absorption.cpp](/Users/admin/Documents/xhdfe-xfe/src/fe_absorption.cpp:9120).
- Automatic MLSMR selection is disabled whenever a GPU backend is requested; the request
  resolves to a sweep fallback:
  [hdfe_regressor_v11.cpp](/Users/admin/Documents/xhdfe-xfe/src/hdfe_regressor_v11.cpp:7892).
- No `Mlsmr`/`mlsmr` implementation occurs in the CUDA or Metal translation units. The
  CUDA adaptive solver is a different block-PCG/MAP path.

A composite workflow can use CUDA for one phase and then CPU MLSMR for another. For
example, a CPU MLSMR repair may follow an unavailable or non-GPU path. That does not make
MLSMR GPU-executed.

## What the paper says and how fair it is

### Credit the paper does give

The paper is not silent about `within`:

- Page 5 disclaims a new estimator.
- Pages 10--11 describe xhdfe's modified LSMR path, then identify a "closely related
  graph-preconditioned Krylov solver" developed by Fischer and Schröder, including its
  reusable factor-pair preconditioner.
- Page 18 exposes `within` as an alias for `mlsmr` and accurately says MLSMR is CPU-only
  and standard-FE-only.
- Page 30 thanks Alexander Fischer for sharing updates, names Kristof Schröder, describes
  modified LSMR plus additive Schwarz, and links the `within` repository.
- Page 31 cites the Fischer--Schröder mimeo.

That is substantive acknowledgment and should be recognized as such.

### Where the framing falls short

The difficulty is the relationship the words imply. "Closely related" and "helpful for our
ongoing work" describe adjacent or inspirational work. They do not tell the reader that
xhdfe's `mlsmr` is an implementation/adaptation of the same solver architecture. In light
of the exact code fingerprint, direct aliases, and timeline, this is materially incomplete.

The benchmark presentation is also too opaque:

- Page 26 says the default CPU absorber may select modified LSMR, but Tables 7 and 8 do
  not report `absorption_method_used` for any row.
- The repository's benchmark showcase identifies the exact
  `akm_v02_firstreg` common-seniority specification and states that its auto selector
  promotes that graph to MLSMR:
  [xhdfe_r_showcase.qmd](/Users/admin/Documents/xhdfe-xfe/r/examples/xhdfe_r_showcase.qmd:787).
  Its benchmark helper actually records the selected method
  [at line 724](/Users/admin/Documents/xhdfe-xfe/r/examples/xhdfe_r_showcase.qmd:724), but
  the paper table omits that column.
- It is therefore strongly supported that the common-seniority **CPU** rows use MLSMR,
  although the paper's underlying run logs would be needed for row-by-row confirmation.
  Those rows include the paper's comparisons in which xhdfe CPU is reported as 19--22%
  faster than pyfixest.
- The common-seniority **CUDA** rows do not use MLSMR; they use GPU sweep/PCG machinery.
- The firm-specific-seniority model has heterogeneous slopes, so all of its CPU and GPU
  rows use sweep-based absorption, not MLSMR.

This materially narrows the attribution issue. It would be unfair to imply that all xhdfe
speedups come from `within`: the GPU results and heterogeneous-slope results do not. It is
also unfair to fold the likely MLSMR-assisted CPU rows into undifferentiated headline xhdfe
results while describing the upstream architecture only as "closely related."

The paper's mechanism tables isolate sweep extrapolation, threading, and GPU execution;
they do not isolate MLSMR. There is no xhdfe-sweeps versus xhdfe-MLSMR ablation and no
xhdfe-MLSMR versus `within` comparison. Omitting `within` from a full-estimator comparison
is not automatically improper if `within` is treated as a solver library, but the selected
absorber and algorithmic provenance still need to be disclosed.

### Overall fairness judgment

- **Fair:** the authors disclaim a new estimator, acknowledge both `within` authors,
  describe the technical relationship, cite the mimeo, and correctly separate CPU MLSMR
  from GPU execution. xhdfe also contains substantial original integration and performance
  engineering.
- **Not fully fair:** the wording understates architectural dependence and the benchmark
  tables hide which CPU results use the adapted solver.
- **My bottom line:** the paper is not deceptive or attribution-free, but it should be
  revised. "Closely related" is not the clearest factual characterization of the code that
  ships under the aliases `within` and `within-additive` and reproduces the full distinctive
  architecture and defaults.

## Recommended revision

A technically accurate sentence would be:

> The xhdfe `mlsmr` absorber is a standalone C++ adaptation of the
> MLSMR--additive-Schwarz fixed-effects solver architecture developed in the `within`
> project and shared with us on 17 March 2026. It follows `within` in constructing
> factor-pair connected-component domains, applying symmetric partition-of-unity
> weighting, and solving the local bipartite systems through Schur reduction and
> approximate Cholesky. xhdfe adds its own C++ kernels, scaling and fallback policies,
> automatic method selector, multi-interface integration, and convergence/certification
> machinery.

I would ask the authors to make four changes:

1. Replace "closely related" with "based on," "adapted from," or "an independent C++
   implementation of," depending on the provenance they can substantiate.
2. Cite a `within` release/commit and repository in the technical MLSMR section, rather
   than leaving the software URL only in the acknowledgments.
3. Add an `absorber used` column or footnote to every benchmark row selected through
   `auto`, and identify the CPU common-seniority rows that used MLSMR.
4. Add an ablation separating xhdfe's sweep, MLSMR, threading, and GPU contributions; if
   feasible, include a solver-level comparison with `within` on the same demeaning task.

## Suggested outreach note

> We appreciated that the draft names both of us, cites our mimeo, and links `within`.
> After reviewing the xhdfe source, however, we think "closely related" and "helpful for
> our ongoing work" understate the implementation relationship. xhdfe's `mlsmr` path
> reproduces the distinctive `within` architecture--factor-pair connected-component
> Schwarz domains, two-sided \(1/\sqrt c\) weights, larger-side Schur elimination,
> sampled clique-tree Schur construction, approximate Cholesky with the same defaults,
> and paired M-weighted modified-GK reorthogonalization--and xhdfe exposes `within` as a
> direct alias. The chronology is also that we shared the algorithm on 17 March, opened
> the public MLSMR PR on 28 April, and merged it on 29 April; xhdfe's source dates its
> internal port to 24 June.
>
> Would you be willing to describe the xhdfe MLSMR absorber explicitly as a standalone
> C++ adaptation/implementation of the solver architecture developed in `within`, while
> separately identifying xhdfe's own C++ integration, automatic selection, fallbacks,
> threading, and GPU sweep work? We would also appreciate a versioned `within` citation
> and an indication of which benchmark rows selected MLSMR. We are happy to agree on
> precise language that credits both the upstream algorithmic work and your independent
> implementation contributions.

The tone should stay factual and collaborative. Sending a short evidence table with the
message will make the requested change concrete and give the authors an opportunity to
correct any provenance detail that is not visible in the public history.

## Licensing caveat

Both projects use MIT licensing. xhdfe's current `NOTICE` names bundled third-party
components but contains no `within` entry. That fact is relevant **only if** xhdfe translated
or copied a substantial portion of `within` source code. Architectural or algorithmic
adoption alone does not establish a software-license notice obligation. The cross-language
structural match is not enough for me to make that legal conclusion. If the authors confirm
that their C++ was translated from `within` source, ask separately that the applicable MIT
copyright/license notice be included.

## Limitations

- I did not inspect the 17 March message, attachment, meeting notes, or any private xhdfe
  branch. The first xhdfe commit is squashed, so its earlier development history is not
  auditable from the public repository.
- I did not reproduce the paper's large benchmarks; the datasets and run logs are not in
  the reviewed files. The benchmark-method inference uses the paper, public dispatch code,
  and xhdfe's own benchmark-showcase comment.
- Similar architecture is strong evidence of derivation/adaptation, not proof of literal
  source copying or intent.
- This is a technical and scholarly-attribution assessment, not legal advice.
