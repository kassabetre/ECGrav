# Changelog

All notable changes to ECGrav are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning is
semantic-ish; breaking changes are called out explicitly.

## [Unreleased]

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
