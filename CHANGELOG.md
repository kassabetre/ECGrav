# Changelog

All notable changes to ECGrav are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning is
semantic-ish; breaking changes are called out explicitly.

## [Unreleased]

### Changed
- **Recalibrated the equilibriation convergence criterion** in all four equilibriators
  (`GraphEquilibriate` ×2, `RandomPureSimplicialComplexMCMCEquilibriate` ×2). The criterion
  compares `sqMeanPairwiseDiff`, a squared difference of 30-point window *means*, against a
  threshold. That threshold was `meanLateVar`, the variance of the *individual* energies —
  larger than the quantity being tested by the effective sample size, so the test reduced to
  "is the autocorrelation time below about 7.5 sweeps". Measured on AR(1) tracks, the ratio it
  computed was `2·tau_int/15`; at the correlation times these chains actually show it passed
  with a 14× margin, and a systematic 1σ offset between half the replicas was declared
  converged 99% of the time.

  The threshold is now `expectedSqDiff`, twice the variance of a window mean, estimated from
  the scatter of the non-overlapping window means inside each track. That is a within-track
  quantity, so an offset between tracks cannot inflate it, and it carries the chain's own
  autocorrelation without needing a correlation time — which is not available until after
  equilibriation. The pass rate on equilibriated tracks is now flat in the autocorrelation, at
  about 50% per check, instead of running from 100% to 4%.

  Effect on the residual disagreement between replicas at the moment equilibriation exits, on
  the one test system with a real transient (`Partition[Range[20],2]`, `labelingChoise 2`):
  **1.24σ → 0.47σ**, for ten extra sweeps. Small systems are unaffected. `eqlT` grows on
  systems that were exiting early — for that system it goes from 31–32 to 59–210 — so **seeded
  runs do not carry across this change**, and any `eqlT` recorded from an earlier version is an
  underestimate.

- **The equilibriation test now reads two observables, over a longer window, once per window.**
  Four changes to the same test, all in the four equilibriators:

  - **Two scalars are watched, not one** — the energy and the vertex count (free-vertex
    complexes), the edge count (fixed-vertex), or the edge density (graphs). A single scalar
    can be flat for reasons unrelated to mixing, and one whole case was blind: the fixed-vertex
    chain at `labelingChoise -> 0` gives every state weight 1, so the energy is identically 0
    on every track. `RandomPureSimplicialComplexMCMCEquilibriate[Partition[Range[30], 3],
    True, 0]` used to return `converged -> True` at the first check for any input, however far
    from mixed. Both scalars are already carried on every replica, so the second costs nothing.
  - **A scalar that never varied gets no vote, and if none varied that is now said out loud**
    via the new `GraphEquilibriate::nosignal` / `…MCMCEquilibriate::nosignal`. It still returns
    `converged -> True` — a chain whose observables never move is stationary in them — but the
    silence is no longer passed off as evidence. This replaces the `Tr[sqMeanEMat] < 1e-5`
    metastability escape, which read only the diagonal of the matrix, i.e. each track's own
    drift, and so declared four tracks frozen at four *different* energies converged however
    far apart they sat. That also removes the `1e-5` / `1e-6` disagreement between overloads.
  - **`inWinLength` 30 → 100.** That window alone sets the resolution: an offset between
    replicas is detectable down to about `sigma*Sqrt[2g/inWinLength]`. Measured over AR(1)
    tracks, the pass rate at a 1σ offset falls from 0.18 to 0.04. It is also the batch length
    of the threshold estimate, which is what keeps the comparison exact at any autocorrelation
    — both sides average over the same number of sweeps, so any batching bias cancels.
    `outWinLength` was left alone: it only sets how many batches the threshold is estimated
    from, lengthening it was measured to buy nothing, and it is what every run pays before the
    first check.
  - **Checked once per `inWinLength` sweeps instead of every sweep.** Consecutive checks shared
    all but one of their points, so the test was re-read hundreds of times over the same data
    and stopped at the first favourable fluctuation. The acceptance is `1.5 ×` the null
    expectation rather than exactly it, chosen off the measured curve: at that tolerance 85% of
    equilibriated checks pass, so a run exits in 1.1 checks rather than 2, while the pass rate
    at a 1σ offset is still 0.04.

  `::noconv` gained a slot and now names the observable that agreed least well.

  **This is the change that costs runtime, and most of that cost is an accounting artifact.**
  Equilibriation itself got 1.0–1.6× longer, which is the intended effect. But `eqlT` grew
  3–6×, because `eqlT = numsweeps - outWinLength + inWinLength` carries a flat `+70` from the
  window change on a base of ~31, and the correlation-time stage runs `10*eqlT` sweeps. Same
  seeds, before → after: p2 M4 `labelingChoise 2` converged at sweep 401 → 501 while its
  correlation-time run went 310 → 2010 sweeps; p2 M10 `labelingChoise 2` went 480 → 760 sweeps
  while its run went 2100 → 4010. End to end on `RandomPureSimplicialComplexMCMC` at p=2, M=4,
  `labelingChoise 2`: 0.9 s → 16 s.

  The root cause is that `eqlT` does double duty — a burn-in report *and* the length scale for
  measuring autocorrelation — and those want different numbers: measurement length should scale
  with the correlation time, not with burn-in. `10*eqlT` only looked reasonable while the window
  arithmetic pinned `eqlT` near 31. See the next entry. Seeded runs do not carry across this one.

- **Shortened and capped the correlation-time measurement run**, in all six routines that set
  one: `GraphComputeCorrelationTime` (×4, which used `30*eqlT`) and
  `RandomPureSimplicialComplexMCMCCorrelationTime` (×2, `10*eqlT`). Both now use
  `Min[5*eqlT, $ECGravMaxCorrelationSweeps]`, the cap being a new global with a default of
  1000, distributed to subkernels alongside `$ECGravMaxEquilibriationSweeps`.

  The multiple came down because it was never measuring what it was scaled to: at the
  correlation times these chains show, roughly 2 to 4, resolving one needs on the order of 20
  correlation times of run, i.e. under a hundred sweeps, and `30*eqlT` was spending 930. The
  cap then stops the growth the previous commit introduced — at the larger `eqlT` the new
  equilibriation test reports, the complex side would have run 4010 sweeps and the graph side
  3030. Together: the graph side goes 930 → 155 sweeps at its old `eqlT`, the complex side
  310 → 155, and neither can now exceed 1000. End to end on
  `RandomPureSimplicialComplexMCMC` at p=2, M=4, `labelingChoise 2`: **16 s → 6.7 s**, against
  0.9 s before any of this work, the remainder being genuine extra equilibriation.

  Note `lastLag = numsweeps - 10`, so at `5*eqlT` an `eqlT` below 3 leaves no lag range and
  every reported correlation time falls to the floor of 2. Nothing clamps this.

  A cap that silently truncates would be worse than the coupling, so the routines now check
  what they found against how long they watched: if the reported correlation time exceeds a
  twentieth of the run — the usual minimum for resolving one — they say so through the new
  `::shortrun`, naming the run length and the cap. It does not fire on ordinary runs at any
  `labelingChoise`, or on the graph driver.

  Still heuristic: the run length should be set by the measurement need, `numsweeps >= 20*tau`
  discovered adaptively, not by the burn-in at all. At the correlation times these chains show
  that would usually run shorter than the cap.

## [1.6.0] - 2026-08-10

Fully unlabeled pure complexes. The counting trio is completed by
`NumUnlabeledPureComplexes`, which counts the isomorphism classes -- neither vertices nor
facets labelled -- and the unlabeled random generator, which had been rejecting with an
acceptance rate that collapsed factorially, is now exact and rejection-free: 11.4 s to 0.82 ms
per sample at `{3,4,10}`. Two top-level specifications, `UnlabeledCount.md` and
`UnlabeledSampler.md`, record the derivations alongside `Theory.md`. No public function was
removed and no argument order changed.

Note for anyone relying on reproducibility: `RandomUniformUnlabeledPureSimplicialComplex`
returns a *different* complex for a given `SeedRandom` than it did at 1.5.0. The distribution
is the same -- uniform over isomorphism classes -- but the algorithm is not, so seeded
sequences do not carry across the version.

### Changed
- **`RandomUniformUnlabeledPureSimplicialComplex` is now exact and rejection-free.** It drew
  from the vertex-labeled generator and accepted with probability `|Aut|/n!`, which collapses
  factorially. Measured against the previous implementation: **146 → 0.32 ms** per sample at
  `{3,4,6}` (458×), **3.31 s → 0.68 ms** at `{3,4,8}` (4,800×), and **11.4 s → 0.82 ms** at
  `{3,4,10}` (**14,000×**). Same distribution — uniform over isomorphism classes.

  It samples a fixed pair `(sigma, S)` of the `S_n` action, where `S` is a `sigma`-invariant
  covering `M`-set of facets: every orbit contributes exactly `n!` such pairs whatever its size,
  so dropping `sigma` leaves the isomorphism class uniform, with no automorphism weighting. The
  draw picks a cycle type, then a *size profile* — how many `<sigma>`-orbits of each size the
  facet set uses — then the orbits themselves, each weighted by its number of completions.

  Sampling the profile up front is what makes the last step tractable. Read the count as orbits
  of separating incidence tableaux under `S_M` rather than as `S_n` orbits of facet sets, and the
  acting permutation's cycle type *is* the size profile. Without that it would be hidden state:
  the completion count depends on how many orbits **of each size** are already used, not merely
  how many — which is where this differs from the facet-labeled sampler, whose objects each cost
  one unit of the budget.

  Covering is imposed during the draw rather than by rejecting afterwards. Rejecting is exact and
  much simpler, but its expected trial count is exactly `A(p,M,n)/U(p,M,n)` — fine mid-range,
  hopeless at the top of the vertex range where the facets must nearly partition `[n]`: 66 trials
  at `{3,4,12}`, 4,279 at `{3,6,18}`, over 10⁷ at `{4,7,28}`.

  Verified by three exact identities rather than by sampling alone: the cycle-type weights sum to
  `n!·NumUnlabeledPureComplexes[p,M,n]`; the profile marginals sum back to the per-cycle-type
  covering count; and at the identity cycle type that covering count reproduces
  `NumVertexLabeledPureComplexes`. Together these force the class probability to `1/U` exactly,
  which was also computed directly — not sampled — and came out exactly `1/U` for every class on
  seven parameter sets. Plus chi-square uniformity and validity across the range including
  `n = pM`.

  With the vertex count free, `{p,M}` now draws `n` from `NumUnlabeledPureComplexes` directly.
  The old body drew it from the **vertex-labeled** counts and let rejection repair the
  difference. `$RandomUnlabeledParallelThreshold` is 500, measured: nearly all the work is in
  memo tables keyed on the cycle type and each subkernel rebuilds its own, so parallel pays much
  later here than for the other generators.

### Added
- `NumUnlabeledPureComplexes[p, M, n]` and `NumUnlabeledPureComplexes[p, M]` — the number of
  fully unlabeled pure complexes: `M` pairwise distinct `p`-subsets of an `n`-element vertex
  set covering it, up to relabelling the vertices *and* permuting the facets. These are the
  isomorphism classes, the count both labelled ones refine. Completes the trio alongside
  `NumVertexLabeledPureComplexes` and `NumFacetLabeledPureComplexes`, same argument order.

  The `S_n × S_M` action collapses to `S_n` alone — forgetting the facet order maps tuples
  onto sets, and `S_M` acts simply transitively on each fibre because the facets are
  distinct — so this is one Burnside over the vertex permutations, with `|Fix(sigma)|` the
  number of `sigma`-invariant covering `M`-sets. Those are unions of `<sigma>`-orbits of
  `p`-subsets, giving `[z^M] Product_d (1 + z^d)^n_d`, and covering is restored by the same
  isolated-vertex differencing `U(n) = A(n) - A(n-1)` the facet-labeled count uses.

  The orbit counts `n_d` are never formed. Taking the logarithmic derivative of that product
  turns it into a recurrence whose coefficients are the fixed-subset counts directly —
  `m c(m) = Sum_k l(k) c(m-k)` with `l(k) = f(k)` for odd `k` and `f(k) - 2 f(k/2)` for even
  — which is Newton's identity, using `Sum_{d|m} d n_d = f(m)`. That removes the Möbius
  inversion and, with it, binomial expansions over enormous `n_d`. Same values everywhere
  (61,172 helper-level comparisons against the Möbius route, zero mismatches), and faster
  where `n` dominates, which is the regime that limits this function: 0.28 → 0.21 s at
  `{2,10,20}`, 1.94 → 1.57 s at `{4,8,28}`. It costs a few ms more when `M` greatly exceeds
  `n` (0.05 → 0.08 s at `{3,100,12}`), where the O(M²) convolution outweighs the handful of
  polynomial multiplies it replaces.

  Verified against a brute-force oracle that enumerates the covering `M`-sets and groups
  them into `S_n` orbits explicitly (29 parameter sets, which reproduce the vertex- and
  facet-labeled counts from the same enumeration); against the Burnside average computed
  from explicitly enumerated fixed sets rather than from the orbit polynomial; and
  externally at `p = 2`, where these are the unlabelled graphs with no isolated vertex —
  `1, 0, 1, 2, 7, 23, 122, 888` for `n = 0..7`.

  Unlike the facet-labeled count there is no restriction to cycle types with parts `<= p`:
  with the facets permutable a long cycle can be covered by a facet whose orbit walks around
  it. Every cycle type contributes, so cost grows like `PartitionsP[n]` — `n`, not `p` or
  `M`, is the limiting argument (about 1 s at `n = 25`, 3 s at `n = 30`). Results are
  memoised for the session behind the shared `NumPCClearCache[]`.

## [1.5.0] - 2026-08-06

Pure-complex counting and random generation. Three of the random generators were silently
returning malformed output above a sample-count threshold; that is fixed, and the
facet-labeled generator no longer uses rejection at all. The counting functions are renamed
to say which labelling they count, and a facet-labeled count is added. No public function was
removed, and `NumPureComplexes` still works, with the same arguments in the same order.

### Added
- `NumFacetLabeledPureComplexes[p, M, n]` and `NumFacetLabeledPureComplexes[p, M]` — the
  number of facet-labeled pure complexes: `M` labeled, pairwise distinct `p`-subsets of an
  `n`-element vertex set covering it, counted up to permutation of the unlabeled vertices.
  Equivalently the separating `(p, n, M)` incidence tableaux. Verified three independent ways
  — against the author's `Facet-labeled-count.nb` (11 values plus its 92-digit
  `sF[3,10,50]`), against direct enumeration of the defining description on 12 parameter
  sets, and against the Stirling identity `B = Sum_k StirlingS2[M,k] F` on brute-forced
  counts of *all* tableaux.
- `NumVertexLabeledPureComplexes` — the new name of `NumPureComplexes` (see Changed).

### Changed
- **`NumPureComplexes` → `NumVertexLabeledPureComplexes`.** The old name said nothing about
  which of the two labellings it counted, which matters now that both exist. Arguments are
  unchanged — `(purity, facet order, vertex count)`, the facet order now written `M` rather
  than `q` — and `NumPureComplexes` remains as a forwarding alias, so existing code keeps
  returning exactly what it did.
- **`RandomUniformFacetLabeledPureSimplicialComplex` is now exact and rejection-free.** It
  drew from the vertex-labeled generator and accepted with a stabiliser-ratio weight, whose
  acceptance rate collapses with `n`. It now samples a fixed pair `(sigma, x)` of the `S_n`
  action: every orbit contributes exactly `n!` such pairs whatever its size, so dropping
  `sigma` leaves the orbit uniform. Same distribution, verified by chi-square over the
  enumerated isomorphism classes. Per sample: 2.09 → 0.48 ms at `{3,4,6}`, 50.9 → 0.86 ms at
  `{3,4,8}`, **516 → 0.99 ms at `{3,4,10}`** — nearly flat in `n` where the old cost grew
  factorially. With the vertex count free, it now draws `n` from
  `NumFacetLabeledPureComplexes` directly rather than from the vertex-labeled counts.
- `NumVertexLabeledPureComplexes` computes one whole row of `n` at a time, bottom-up over
  `M`, and memoises the rows for the session: **68–250×** on a cold two-argument call, 4–6×
  on a single three-argument one. It also no longer hits `$RecursionLimit` — the old form
  died past `M = 700`, the new one handles `M = 2000`. The row cache is capped
  (`$NumPCMaxCachedQ`, default 150) and clearable.
- `RandomVertexLabeledPureSimplicialComplex` is **4.9–8.1×** faster: the composition weights
  depend only on the facet index and the covered-vertex count, so they are built once per
  pair instead of once per sample (82–89% of the old runtime), a provably dead feasibility
  re-draw is gone, and one random permutation replaces the covered/uncovered bookkeeping.
- Random generators now parallelise only when kernels are already running. `ParallelTable`
  otherwise launches them itself, which cost ~3.7 s — far more than the whole serial draw at
  any relevant sample count.

### Fixed
- **`RandomVertexLabeledPureSimplicialComplex`, `RandomUniformUnlabeledPureSimplicialComplex`
  and `RandomUniformFacetLabeledPureSimplicialComplex` returned malformed output above their
  parallel thresholds** — facet lists containing unevaluated `RandomSample[...]` expressions,
  for 100% of samples, silently. Their `ParallelTable` distributed only `ECGrav\`Private\``,
  so the subkernels had no `NumPureComplexes` and the composition weights never evaluated.
  The two uniform generators tripped at 21 samples; the vertex-labeled one above 10^4.
- Empty sample spaces now give a message and `$Failed` instead of malformed complexes, and
  the one-complex overloads propagate that rather than `First[]`-ing into it.
- `NumVertexLabeledPureComplexes[p, 0, n]` returned `Indeterminate` for every `p` and `n`
  (there was no `M == 0` branch, so it divided by `M`). `M = 0` is the empty complex: the
  count is 1 at `n = 0` and 0 otherwise — also the seed that reproduces the `M = 1` row
  through the recurrence.
- The two-argument `NumVertexLabeledPureComplexes[p, M]` returned an unevaluated
  `Total[Table[...]]` when no vertex count carried enough distinct `p`-subsets; that sum is
  empty, so it returns 0.
- `RandomUniformFacetLabeledPureSimplicialComplex`'s usage message described its output as
  "unlabeled" complexes. It produces facet-labeled ones.

### Documentation
- Reference pages for `NumVertexLabeledPureComplexes` and `NumFacetLabeledPureComplexes`;
  the `NumPureComplexes` page now documents the alias. The old page's usage text was stale
  independently of the rename — it still described a recursive implementation replaced in
  this release.
- 115 public symbols, 115 reference pages, verified one-to-one.

### Tests
- 87 → **132**. New coverage for both counting functions (including independent
  brute-force oracles and the Stirling identity), the alias, the parallel branches of all
  three fixed generators, uniformity of the vertex-labeled, unlabeled and facet-labeled
  samplers, the Burnside weight identity behind the facet-labeled sampler, and the
  empty-sample-space guards.

## [1.4.0] - 2026-08-02

Performance release. The Monte Carlo hot path no longer constructs `Graph` objects: the
three functions it leans on are computed straight from the adjacency matrix, which is
8–22× faster at 8 vertices and takes a `GraphEquilibriate` run there from 324 s to 21 s on
an identical trajectory. No public function was added, removed or renamed, and every
changed function was verified to return identical results. The one behavioural change is
the new progress report from `GraphEquilibriate`, called out under Added.

### Changed
- `HyperDeg[Amat_List, clq]` and `delH2dCombManifold` are now computed entirely from the
  adjacency matrix and build no `Graph` object. Both are pure optimisations: outputs are
  unchanged (verified `SameQ`-identical over 8142 `HyperDeg` and 3074 `delH2dCombManifold`
  calls across 95 graphs — random over sizes 3–9 and four densities, plus empty, complete,
  cycle, path, star, wheel, bipartite, disjoint-union and singleton graphs, both signs of
  `J`, every vertex pair, and degenerate/duplicate cliques), and the 87-test suite passes.

  - `HyperDeg[Amat, clq]` tests the clique condition as a sum over the induced submatrix
    and counts common neighbours as `Total[Times@@Amat[[clqV]]]`. It previously ran
    `CompleteGraphQ[Subgraph[...]]` and then delegated to the `Graph` overload, which ran
    that same test a second time. **0.88 ms → 0.040 ms** (22×) at 8 vertices.
  - `delH2dCombManifold` replaces `Sph`/`Subgraph` with adjacency-row intersections,
    `FVector` (which called `FindClique`) with vertex/edge/triangle counts read off the
    submatrix — only the first three f-vector entries survived its `PadRight[...,3]` — and
    the `HyperDeg` sum over `Tuples` with a single dot product. **4.18 ms → 0.32 ms** (13×).
  - `H2dCombManifold` speeds up as a side effect, since it calls the matrix `HyperDeg`:
    **12.4 ms → 1.47 ms** (8.4×).

  End to end this takes `GraphEquilibriate` at 8 vertices from 331 s to 122 s with an
  identical trajectory (same `eqlT`, `minEnergy` and final energy). Callers that pass
  `AdjacencyGraph[am]` to `HyperDeg` rather than `am` still take the slower `Graph` path.

- `EulerChi[Amat]` no longer builds a `Graph` or calls `FindClique`. It accumulates the
  alternating clique sum directly: with `g[C]` the signed count of cliques inside `C`,
  `g[C] = g[C without v] - g[C intersect N[v]]` and `chi = 1 - g[V]`, so one memoised
  recursion over vertex subsets (as bitmasks) gives `chi` without listing a clique.
  **0.66 → 0.071 ms** at 8 vertices, **3.86 → 1.22 ms** at 30. Output is unchanged —
  verified identical over 226 adjacency matrices (sizes 1–10 across seven densities
  including 0.0 and 1.0, plus complete, cycle, path, star, wheel, bipartite, grid,
  Petersen, disjoint unions, all-zero, `{{0}}`, `{{1}}`, `{{1,1},{1,1}}`), with `chi`
  spanning −5 to 10; the facet-list and `Graph` overloads are untouched.

  The `Amat == {{1}}` special case was removed as redundant: the chosen vertex is cleared
  from the candidate set before intersecting, so a diagonal self-loop cannot make the
  recursion revisit it and any 1×1 matrix yields 1.

  This also gives callers a cheap edge-toggle delta. Toggling `{i,j}` adds or removes
  exactly the cliques containing both `i` and `j` — each of them `{i,j}` plus a clique of
  the link — so the f-vector shifts by two and the alternating sum collapses to
  `s (EulerChi[link] - 1)` with `s = -(2 Amat[[i,j]] - 1)`, replacing an `FVector` of a
  `Subgraph` with `EulerChi` of a submatrix.

### Added
- Both `GraphEquilibriate` overloads now emit a `PrintTemporary` progress line reporting
  the inverse temperature and Hamiltonian parameters they were called with. This is the
  one user-visible behavioural change in the release: the cell is removed automatically
  when the evaluation finishes, so it leaves nothing in a saved notebook, and it does not
  display at all from the parallel subkernels that `GraphParallelTempering` and the
  temperature schedulers run it in.
- Private helper `EulerChiAM` — the adjacency-matrix Euler-characteristic recursion behind
  the `EulerChi[Amat]` overload above.
- Private helpers `ConnectedComponentCountAM` and `LinkComponentCountAM` — component
  counts of a graph and of a vertex link straight from an adjacency matrix, by boolean
  transitive closure (repeated squaring), matching
  `Length[ConnectedComponents[AdjacencyGraph[...]]]` including isolated vertices.

## [1.3.1] - 2026-07-29

Documentation and tooling only — no changes to the package's functions or public API,
and the installable paclet is unchanged in content from 1.3.0.

### Added
- `Theory.md` — a specification of the model (configuration space, observables with
  exact `SpectralDim` / `HausdorffDim` / `EulerChi` / Betti-number definitions, the
  Hamiltonians, and the statistical mechanics), plus a marked outline for the physics
  background. Linked from the README.
- `CITATION.cff` — citation metadata (enables GitHub's "Cite this repository"); the
  project is registered with Zenodo for a DOI.
- `build.wls` — a reproducible, self-verifying paclet build script
  (`wolframscript -file build.wls`) that scripts the manual build steps.

### Changed
- README: restructured Installation into three options (load source / install a
  released paclet / build from source) and refreshed the version references.

## [1.3.0] - 2026-07-28

Phase 3 (quality & CI). The three known bugs are fixed, test coverage expanded
(63 → 87 tests), GitHub Actions CI runs the suite on every push/PR, all ~120
`ReturnAmbiguous` code-analysis warnings are cleared, and the subpackage context
convention was simplified to a single `ECGrav`` public namespace. No public API
was removed; the two `CombinatorialSphereQ` / `DSphereQ` changes below are the
only behavioural (breaking) changes.

### Added
- `ClosedCombinatorialManifoldQ` — tests whether a pure simplicial complex (given as a
  facet list) is a combinatorial manifold *without boundary*: connected, every
  codimension-1 face contained in exactly two facets, and every vertex link recursively
  a closed combinatorial manifold. This is the predicate the old `CombinatorialSphereQ`
  actually implemented — a torus, Klein bottle, and sphere all return `True`.
- `ClosedDGraphQ` — the graph/clique-complex analogue: tests whether a graph's clique
  complex is a closed combinatorial manifold (the predicate the old `DSphereQ` actually
  implemented). Reference pages added for both new symbols.
- Regression tests: `CombinatorialSphereQ` / `ClosedCombinatorialManifoldQ` and `DSphereQ`
  / `ClosedDGraphQ` over the reference surfaces and their clique graphs, plus a
  `GraphParallelTempering` smoke test exercising the beta-history boundary (a previously
  untested function). Suite is now 73 tests.

### Changed
- **Breaking:** `CombinatorialSphereQ` now tests genuine sphere-ness. It returns `True`
  only for a combinatorial sphere — a closed combinatorial manifold with the Euler
  characteristic of a sphere of its dimension (χ = 1 + (−1)^d). Previously it returned
  `True` for *any* closed combinatorial manifold, so a torus, Klein bottle, or projective
  plane wrongly returned `True`. Callers that want the old behaviour should use
  `ClosedCombinatorialManifoldQ`. For surfaces the new test is exact; in dimension ≥ 3 the
  Euler-characteristic condition is a necessary (homology-level) filter only, since
  recognizing PL spheres is not algorithmically decidable in general.
- **Breaking:** `DSphereQ` now tests whether a graph's clique complex is a combinatorial
  manifold *homeomorphic to a sphere* — a closed d-graph (see `ClosedDGraphQ`) with the
  Euler characteristic of a sphere. Previously it returned `True` for any closed d-graph,
  so a graph whose clique complex is a torus wrongly returned `True`. `DGraphBoundary`'s
  interior-point test was repointed to `ClosedDGraphQ` to preserve its behaviour (the two
  agree on all genuine geometric graphs).
- Internal refactor, no behaviour change: the kernel's early-exit `Return`s are now
  explicit tagged `Catch`/`Throw`, clearing all ~120 `ReturnAmbiguous` code-analysis
  warnings (116 in `PureComplexes.wl`, 4 in `MCSims.wl`). `Return[expr]` inside
  `If`/`With`/`While` has a dynamically-scoped exit point that the code analyzer flags as
  ambiguous; each affected function body is now wrapped in `Catch[…, "tag"]` with its
  guards throwing that unique tag, and seven whole-body `If[c, Return[a], Return[b]]`
  predicates were simplified to `If[c, a, b]`. The full 87-test suite passes unchanged.
- Internal refactor, no behaviour change: the two subpackages (`MCSims.wl`,
  `PureComplexes.wl`) no longer open their own `BeginPackage["ECGrav`Sub`"]` context.
  They are now loaded as plain includes into the single `ECGrav`` public context that
  `ECGrav.wl` establishes, and the ~890 hand-qualified `ECGrav`Foo` prefixes were dropped.
  Public functions resolve through their up-front `::usage` declarations; any undeclared
  symbol falls into `ECGrav`Private`` and is automatically private, so adding a function
  without a usage message can no longer make it silently vanish from the API. All 113
  public symbols remain present, documented, and `Protected`; the suite stays at 87/87.

### Fixed
- `CombinatorialSphereQ[torus]` (and Klein bottle, projective plane) no longer returns
  `True` — see the split above.
- `DSphereQ` no longer returns `True` for a graph whose clique complex is a torus (or any
  other closed non-sphere manifold) — see the split above.
- `ExactExpectationValue` no longer emits a stray `Part::partd` when passed a list of
  observable *functions*. The numeric-observable overload's pattern test now uses
  `MatrixQ[obsValues, NumericQ]` instead of `NumericQ[obsValues[[1,1]]]`, so dispatch no
  longer evaluates `{func,…}[[1,1]]`. Returned values are unchanged.
- `GraphParallelTempering` and `GraphMultiHistogram` no longer emit
  `Part::take: Cannot take positions -N through -1 …` when a replica's accepted-swap count
  reaches the sweep count `NN`: the rolling beta-history buffer trim now clamps the slice
  to the buffer length. Returned results are unchanged.
- The Monte Carlo drivers (`GraphParallelTempering`, `GraphMultiHistogram`,
  `SimulatedAnnealing`, `SGradDescent`, `GradDescent`) no longer emit `Mod::indet` /
  `Mod[_, 0]` on short runs: the progress-print interval `Floor[NN/5]` is now
  `Max[Floor[NN/5], 1]`.

## [1.2.0] - 2026-07-17

First openly published, MIT-licensed release — a cleaned-up, tested, and
documented line of development. The earlier codebase is preserved on the
`v1.1.0` and `v1.0.0` tags (whose history is unrelated to this one).

### Added
- MIT `LICENSE`.
- `README.md` with an overview, install instructions, verified quickstart
  examples, and an API map by theme.
- `ROADMAP.md` tracking documentation and release work.
- `Tests/` — a `.wlt` regression suite (63 tests) runnable via `TestReport`,
  loading the package from source. See `Tests/README.md`.
- `.gitignore` and `.gitattributes` (build tree ignored; the built `.paclet`
  kept and marked binary; Wolfram source tagged for GitHub language stats).
- Installable `build/ECGrav-1.2.0.paclet` artifact.

### Changed
- **Breaking:** `GraphParallelTempering` — the `minEtoBeat` (minimum-energy-to-beat)
  argument moved to the **last** position in the overloads that take it
  explicitly, making the argument order consistent across overloads. Update any
  calls that previously passed it as the third argument.
- `HIsing` reimplemented to compute its sums directly with `Total` instead of
  building an all-ones matrix — faster, with identical results.
- Monte Carlo progress messages now use `PrintTemporary` (self-clearing) instead
  of `Print`.
- Warning conditions now raise proper `Message[...]` calls instead of `Print`.

### Fixed
- `Branchedness` no longer misreads a facet list as an adjacency matrix (which
  leaked a `Null` into the result); facet-list input is routed through
  `GraphFromCliques`, matching how `RmsPurity` handles its input.
- `LowEnergyStates` no longer divides by `$KernelCount`, so it no longer fails
  with `Power::infy` when called before `LaunchKernels[]` (now `Max[$KernelCount, 1]`).
- Restored the diagnostic plots (`ListLinePlot` / `Plot`) emitted by
  `GraphEquilibriate`, `GraphComputeCorrelationTime`, `GraphMultiHistogram`,
  `GraphCTLSchedule`, `ConstrainedProbConjugateField`, and the density-plot
  functions, which had been removed during an earlier debug-print cleanup.
- Fixed a variable-shadowing bug in `GraphMetropolis`: a local `N` shadowed the
  built-in numeric-evaluation symbol and has been renamed to `vCount`.

### Removed
- Dead `LandauFreeEnergy` declaration (usage/error messages with no
  implementation behind them).
- Extensive debug `Print` scaffolding throughout the Monte Carlo code.

### Known issues
- `CombinatorialSphereQ[torus]` returns `True` (a torus is not a sphere).
- `ExactExpectationValue` emits a stray `Part::partd` message (the result is
  still correct). See `Tests/README.md`.

## [1.1.0] — preserved as tag `v1.1.0`

Pre-1.2.0 codebase (unrelated history; reconstructed from its commit log):

- Added `Visualize2DComplex` and `OrientableCombinatorialManifoldQ`.
- Added parallel-tempering capability with an optional `delH`.
- Bug fixes in `Sph`, `Bll`, and other complex/graph functions.

## [1.0.0] — preserved as tag `v1.0.0`

- Initial version.

[Unreleased]: https://github.com/kassabetre/ECGrav/compare/v1.3.1...HEAD
[1.3.1]: https://github.com/kassabetre/ECGrav/releases/tag/v1.3.1
[1.3.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.3.0
[1.2.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.2.0
[1.1.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.1.0
[1.0.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.0.0
