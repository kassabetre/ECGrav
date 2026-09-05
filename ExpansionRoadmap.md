# ECGrav Expansion Roadmap

Planning document for the expansion arc that begins after 1.13.0. `ROADMAP.md` Phases 0–6 record
what the paclet has been; this records where it is going. Phase 6 of that document — the loose ends
from the tempering/plotting cycle — is not superseded: several of its open items are prerequisites
here and are cross-referenced by name.

_Written 2026-09-04, at `2e262bb` (four commits past `v1.13.0`)._

Two papers drive everything, in this order:

| | Paper | Draft due | State space |
|---|---|---|---|
| **P1** | Random pure simplicial complexes at fixed facet count: typical topology and connectivity under vertex-, facet-, and unlabeled equivalence | ~2026-09-25 (2–3 weeks) | pure complexes, d = 1…4 |
| **P2** | Emergent combinatorial manifolds from `H2dCombManifold` with source terms: phase map, order of the manifold-emergence transition, continuum limit by finite-size scaling | ~2026-11-04 (2 months) | **graphs**, 2-complex derived as the clique complex |

---

## 0. Decisions taken

Settled in the 2026-09-04 planning session; recorded so they are not relitigated.

1. **P1 is a probability paper, not an algorithms paper.** It reports how typical topology and
   connectivity change as the facet count *M* grows, separately for each of the three labelings.
   It does **not** attempt rigorous thresholds or asymptotic formulae. A draft already exists.
2. **P1 covers dimensions d = 1, 2, 3, 4** — purity `p = d+1 = 2, 3, 4, 5` in the package's argument
   convention.
3. **P2 keeps the graph state space.** `H2dCombManifold` stays as it is structurally: the state is
   an adjacency matrix, moves are single-edge toggles, and the 2-complex is the clique complex of
   the graph. Its one coupling `J` is **fixed**, and phases are mapped by varying *source terms*
   instead.
4. **P2's source terms couple to three operators**: Euler characteristic, triangle count, and
   boundary-edge count.
5. **Wang–Landau first** for the order of the manifold-emergence transition. Whether a further
   scheme (Berg–Neuhaus recursion, replica-exchange WL, transition-matrix MC) is needed is deferred
   until WL has been run on the real model.
6. **The quantum expansion targets D > 2, k > 2**, which subsumes everything currently implemented.
   The immediate need inside that is **D = 2, k = 3**. It follows both papers.
7. **P2's manifold model is LABELED.** No automorphism weight in the acceptance ratio
   (`UnlabeledVerticesYes = 0`). This simplifies Track C considerably — see §C0.
8. **P2's finite-size scaling runs to N ≤ 21**, larger being cost-prohibitive.
9. **`FacetLabeledCountExamples.md` stays untracked.** Decided 2026-09-04; do not re-raise.
10. **The quantum state is an association from level to facet set** — the user's design, adopted
    with two amendments to the move set. See §D0.

---

## 1. What the audit found

Two facts reshaped this plan and are worth stating up front, because both invert the intuitive
expectation about where the risk lies.

### 1.1 For P1, the samplers are solid and the *observables* are not

The functions that would worry you are the best-covered code in the package. The functions the
paper's figures actually consist of are at zero.

| Function | Tests | Role in P1 |
|---|---|---|
| `NumVertexLabeledPureComplexes` | 30 | ensemble sizes |
| `NumUnlabeledPureComplexes` | 17 | ensemble sizes |
| `NumFacetLabeledPureComplexes` | 16 | ensemble sizes |
| the three exact samplers | 8 / 6 / 9 | draws |
| **`SpectralDim`** | **0** | headline: dimension |
| **`HausdorffDim`** | **0** | headline: dimension |
| **`ConnectedComplexComponents`** | **0** | headline: connectivity |
| **`KpathConnectedComponents`** | **0** | headline: connectivity |
| **`FractionInLargestKPathComponent`** | **0** | order parameter |
| **`PureComplexAutomorphismGroupOrder`** | **0** | ensemble weights |
| `CountHoles`, `FVector`, `RmsPurity`, `AvgKDim`, `FractionInLargestComponent` | 1 each | topology |

Connectivity and dimension are P1's two headline quantities and every function computing them is at
zero or one test. These are deterministic functions with known values on small complexes, so this is
roughly a day of work — and it is the difference between a defensible figure and an undefended one.

**Resolved 2026-09-04 (§A1).** 31 tests added, `Tests/PureComplexes.wlt` now 172/172 green. The
functions themselves turned out to be correct on every known-topology reference; what the exercise
actually found was **three silent wrong-answer paths in how they are called** — see §A1. That is the
better outcome: a wrong function fails loudly on a test, a wrong calling convention does not fail at
all.

### 1.2 The campaign is cheap, and the cost that exists is in the sampler, not the counters

Measured 2026-09-04. The full grid **p = 2…5 (d = 1…4), n = 8…18**, drawing 50 unlabeled complexes
at `M = min(C(n,p)/2, 50)` and measuring `EulerChi`, `FVector`, `ConnectedComplexComponents` on each
plus `PureComplexAutomorphismGroupOrder` on five, completed with **every cell under 0.01 s** — the
counters included. `NumUnlabeledPureComplexes[5, 50, 16]` returns a 105-digit integer instantly.
p = 2 and p = 3 at n = 20, M = 60 are likewise instant. **The whole d = 1…4 campaign region is
affordable.**

Verified non-vacuous: the counts are big integers, samples carry 50 facets, `EulerChi` returned −30
and `PureComplexAutomorphismGroupOrder` returned 1 on a real draw. (A call that silently returns
unevaluated times as 0.00 s — the trap this project has hit before, so the assert matters.)

**The cost that does exist is in `RandomUniformUnlabeledPureSimplicialComplex`, and it depends
sharply on how it is called.** Drawing 5 samples at p = 2 *without* first calling the counter for
that cell:

| n | M = 30 | M = 45 | M = 60 |
|---|---|---|---|
| 16 | 1.16 s | 1.38 s | 1.71 s |
| 18 | 2.43 s | 3.00 s | 3.49 s |
| 20 | 5.14 s | 6.01 s | 7.23 s |

Cost roughly doubles per +2 vertices and grows only weakly with M. In the grid above, the same calls
were free — there the counter for each cell ran first, and the grid walked n upward in one kernel,
so the memoized count rows (`NumUnlabeledPureComplexes` "computes one whole row of n at a time …
and memoizes the rows for the session") were already built. **Campaign-design consequence: call the
counter for a cell before sampling it, walk n upward within a session, and take many draws per cell
rather than few.** That is the difference between seconds and free.

**Correction to an earlier reading.** A first pass over `M = C(n,p)/3` capped at 60, n up to 20, ran
more than ten minutes without output and was initially recorded here as a *cliff above n ≈ 16 at
p ≥ 4*. That was wrong. The output was block-buffered so nothing could be seen until exit, and that
exact configuration — **p = 2…5 at n = 20, M = 60, up to C(20,5) = 15504 candidate facets** — has
since been re-run with each purity in its own process and the counter called first, and **every call
completed at 0.000 s**: counter, 50 draws, observables, and 10 automorphism orders alike. There is
**no cliff** in the campaign region. Whatever the original run was doing (four purities accumulating
memo tables in one kernel, or an unlucky complex in its larger automorphism sample), it is not a
property of the parameters. Process isolation is cheap insurance for a long unattended run —
`row.wls` and `probe20.wls` in this session's scratchpad are working templates — but it is a
convenience, not a response to a measured failure.

---

## 2. Dependency map

```
      A1 observable hardening ─┐
      A2 ensemble spec ────────┼─► A3 campaign runner ─► A4 measurements ─► A5 figures ─► P1 DRAFT
      A0 frontier sizing ──────┘

      B1 H2dCombManifold tests ─┐
      B2 land Homogeneous ──────┼─► B3 source operators ─► B5 phase map ──┐
      B4 companion tests ───────┘                                         ├─► P2 DRAFT
                                 C1 WL spec ─► C2 impl ─► C3 validate ─► C4 order
      B6 tempering defects ──────────────────────────────────────────────┘

      D0 data structure ─► D1 D=2,k=3 ─► D2 operators ─► D3 H-builder ─► D4 D>2,k>2
```

A and B are independent and can overlap. C depends on B1–B3. D follows both papers, though **D0
should be captured now** while the design is fresh.

---

## Track A — P1, the combinatorics paper (weeks 1–3)

**Scope fixed 2026-09-04: the vertex count n is FREE.** P1 samples from `S_V(p,M)`, `S_F(p,M)`,
`S_U(p,M)` — vertex-labeled, facet-labeled and unlabeled p-pure complexes with M facets, over all
vertex counts. So the `[p, M]` counter overloads and the `{p, M}` sampler overloads are the primary
API, and **n itself becomes the first random variable**, not a parameter. The paper covers the
counts, the exact samplers and the MC samplers, then the distributions of: number of vertices,
number of connected components, Betti numbers, and automorphism group order — compared across the
three labelings.

### A-theory. The three ensembles are related exactly — *this is a result, not a measurement*

Each isomorphism class C contributes `n(C)!/|Γ(C)|` vertex-labeled objects and `M!/|Γ_F(C)|`
facet-labeled ones, where `Γ` is the vertex automorphism group and `Γ_F` its image on facets. So as
weights over classes:

| ensemble | weight on class C |
|---|---|
| `S_U` | 1 |
| `S_F` | `1/|Γ_F(C)|` — bounded, and → 1 as complexes become asymmetric |
| `S_V` | `n(C)!/|Γ(C)|` — **factorial in a varying n** |

**This explains the observed F↔U convergence and predicts that V never joins them.** F differs from
U only by a bounded symmetry correction; V is weighted factorially toward large n and so
concentrates near `n_max = pM`. Confirmed in 8 draws at p=3, M=5: V gave n ∈ {12…15} (max 15), F
gave {6…11}, U gave {6…9}.

It also yields an **exact, sampling-free convergence diagnostic**:

> **E_U[1/|Γ_F|] = NumFacetLabeled(p,M) / (M! · NumUnlabeled(p,M))**

computable from two counters, no Monte Carlo. It equals 1 exactly when every complex is
facet-asymmetric, i.e. when F and U coincide. Measured (and validated against 2000 draws from U at
p=3, M=5: sampled 0.5236 vs exact 0.5308):

| p | M=2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|
| 2 | .500 | .300 | .265 | .254 | **.246** | .261 | .277 | .297 |
| 3 | .500 | **.403** | .420 | .531 | .628 | .721 | .788 | .835 |
| 4 | .500 | **.453** | .546 | .707 | .830 | .903 | .940 | .960 |

It dips, then climbs toward 1 — and **the rate depends sharply on purity**: by M = 9 it is 0.96 at
p = 4, 0.84 at p = 3, and only 0.30 at p = 2. *Convergence is dramatically slower for graphs than
for higher-dimensional complexes.* The M = 1 and M = 2 columns are exactly 1 and 1/2 for every p,
which is provable by hand (any two distinct p-subsets admit a swap), so the framework self-checks.

- [ ] Write this up as `EnsembleWeights.md` — it reframes P1's central claim from an empirical
      observation into a theorem plus an exactly computed rate.
- [ ] Extend the table to larger M and to p = 5, and find where each purity crosses, say, 0.99.

### A-stats. Kolmogorov-type testing — **one trap to avoid**

- [ ] **The classical two-sample KS test is conservative under ties, and every statistic here is
      discrete** (vertex count, component count, Betti numbers, and |Aut|, which takes divisor
      values and is violently lumpy). Conservative means inflated p-values, means under-rejection —
      and since the desired conclusion is *"F and U agree"*, i.e. failure to reject, **a
      conservative test is biased toward the answer you want.** A referee will find this. Use
      **permutation-based p-values** for the KS (or Cramér–von Mises / Anderson–Darling) statistic:
      exact under the null regardless of discreteness, and about ten lines of code.
- [ ] **Report the effect size, not the p-value.** A test that fails to reject at each M says
      little — power depends on sample size, so it can be manufactured. The evidence for
      convergence is `D_M → 0`: plot the KS statistic against M with bootstrap intervals, with the
      exact `E_U[1/|Γ_F|]` curve above as the theoretical prediction.
- [ ] Consider reporting **total variation distance** between F and U directly. Since the
      normalizer is known exactly, `d_TV = ½ E_U[|1/(|Γ_F| E_U[1/|Γ_F|]) − 1|]` is an
      importance-sampling estimate with an exact constant — far better conditioned than any
      two-sample test, and it measures the whole distribution rather than one marginal.
- [ ] For |Aut|, compare **log|Aut|**: group orders are multiplicative and heavy-tailed.

### A0. Size the frontier — *mostly done; §1.2 has the numbers*

d = 1…4 at n ≤ 18 is confirmed affordable, and there is no known cliff. What remains:

- [x] **Settled: the whole region is free on the warm path.** p = 2…5 at n = 20, M = 60 — up to
      C(20,5) = 15504 candidate facets — every call (counter, 50 draws, observables, 10 automorphism
      orders) at 0.000 s, each purity in its own process with the counter called first. `n_max` is
      not the binding constraint anywhere in d = 1…4; pick the grid on scientific grounds.
- [ ] If the grid is pushed past n = 20, re-time rather than extrapolate — the cold-path table in
      §1.2 shows cost doubling per +2 vertices, so the warm path is doing real work that could
      eventually bite.
- [ ] Confirm the warm path holds on the real grid — counter first, n ascending, many draws per
      cell — since that is what makes sampling free rather than seconds (§1.2).
- [ ] Write the resulting grid into A2's spec.

### A1. Observable hardening — **DONE 2026-09-04 for P1's four statistics**

**31 tests added to `Tests/PureComplexes.wlt`; the file now runs 172/172 green** (was 141). Every
expected value is known mathematics — S₄ on the tetrahedron boundary, the octahedral group,
b = (1,2,1) on the torus, b = (1,1,0) on the Klein bottle — not recorded behaviour. Coverage now
spans Betti numbers, both automorphism-order functions, `ConnectedComplexComponents`,
`KpathConnectedComponents`, `FractionInLargestComponent`, `SpectralDim`, and the exact
counter identity from §A-theory.

**The functions are correct. The input conventions are the hazard.** Three ways to get a plausible
wrong number, all found by writing these tests, all live in a P1 sweep:

1. **`CountHoles[c, k]` returns b_{k−1}, not b_k.** It is `BettiNumbers[c][[k]]` and that vector
   starts at b₀, so `CountHoles[c,1]` is the *component count*. The usage message calls it "the
   number of k-holes, i.e. the k-th Betti number", which disagrees. A figure labelled b₁ would be
   plotting b₀, and both are small plausible integers. **Use the one-argument `CountHoles[c]`,
   which returns the whole vector.** Out of range the two-argument form leaks an unevaluated
   `{1,1}[[3]]` with a `Part::partw`.
2. **`FractionInLargestComponent` silently misreads a square facet list.** It takes a graph or an
   adjacency matrix; a facet list with M = p is square, so it is accepted as a weighted adjacency
   matrix. `FractionInLargestComponent[{{1,2,3},{4,5,6},{7,8,9}}]` returns **1**; the correct
   answer is 1/3. Non-square facet lists error safely — **only the M = p case is silent, and a
   sweep over M hits it every time.** Go through `GraphFromCliques` always.
3. **`FacetAdjacencyMatrix` is not an adjacency matrix** (§A1 trap box below) — facet sizes on the
   diagonal, intersection sizes off it.

- [ ] **Guard (1) and (2) at the source.** Both are one-line pattern guards: reject a facet list in
      `FractionInLargestComponent` via the existing `PureComplexQ` convention, and either fix
      `CountHoles`'s indexing or fix its usage message. Deferred as a behaviour change rather than
      a test addition — decide before P1's figures, since the tests currently pin the *actual*
      convention so a change to either code or message is caught.
- [ ] Raise `FVector`, `RmsPurity`, `AvgKDim` above a single test each (not P1 headline
      statistics, so lower priority than the four above).
- [ ] **`SpectralDim` at a fixed step s is a LOCAL probe, not an asymptotic dimension.** On a cycle
      it returns 1.2983… identically for n = 20, 60 and 150 — a 6-step walk cannot see global
      structure. That is correct behaviour, but it means P1 must report the step s alongside any
      dimension and must not call the result "the spectral dimension". Decide the s convention in
      A2.
- [ ] **Resolve the dimension-function input question — and mind the trap below.**
      `SpectralDim[Amat, s]` and `HausdorffDim[Amat, Dmat, s]` accept only an adjacency matrix or a
      `Graph`; verified 2026-09-04 that there is **no complex-accepting overload**. So for d ≥ 2
      something must decide *which graph represents the complex*, and the choice changes what
      "dimension" means in the figures. Candidates already in the package: the **1-skeleton**
      (vertices, joined when they share a facet — the identity at d = 1); the **dual graph**
      (facets as nodes, adjacent when they share a (p−1)-face — usually what "spectral dimension"
      means for random walks in the quantum-gravity literature); and **`KFaceDistanceMatrix`**,
      which supplies a distance matrix over k-faces and so pairs directly with `HausdorffDim`'s
      optional `Dmat` argument.

      > **Trap — `FacetAdjacencyMatrix` is not an adjacency matrix.** Despite the name, it returns
      > facet *sizes* on the diagonal (`DiagonalMatrix[Length /@ facets]`) and
      > `|facet_i ∩ facet_j|` off it. Passing it to `SpectralDim` does **not** error: on the
      > tetrahedron boundary it returns 2.516…, a plausible number that means nothing, because the
      > diagonal is read as self-loops and the intersection sizes as edge weights. (A bare facet
      > list *does* error, so that route is safe.) To get the dual graph, threshold it:
      > `Boole[Map[# == p - 1 &, fam - DiagonalMatrix[Diagonal[fam]], {2}]]`.
      > This is precisely the failure the zero test coverage is hiding, and precisely the figure
      > that would have survived review.

### A2. Ensemble and campaign specification — `RandomComplexEnsembles.md`

Following the repo's convention that new work ships with a spec.

- [ ] For each d ∈ {1,2,3,4} (p = 2…5) and each labeling ∈ {vertex, facet, unlabeled}: the M sweep,
      the n values, sample counts, and the error-bar method.
- [ ] The observable list and, for each, exactly what is computed (see A1's signature question).
- [ ] **Record that d = 1 is the control.** Pure 1-complexes are graphs, so the d = 1 row is the
      fixed-edge-count (microcanonical) slice of the ensemble in *Dynamical Quantum Multigraphs*.
      It should reproduce the labeled/unlabeled split established there. That is a known answer to
      check the whole pipeline against before trusting d = 2, 3, 4 — and a reviewer-facing argument
      that the method works.
- [ ] Note that `FractionInLargestComponent` is the same order parameter `s₁` as in that paper.

### A3. Campaign runner

- [ ] Order the sweep to stay on the warm path: counter before sampler per cell, n ascending within
      a session, many draws per cell (§1.2). This is the single biggest constant-factor decision in
      the campaign.
- [ ] Process isolation per grid point is optional insurance for long unattended runs, not a
      response to a measured failure (§1.2) — but note it forfeits the cross-cell memoization that
      makes the warm path free, so isolate per *n*, not per cell.
- [ ] Thin: all real logic stays in tested package
      functions, so the runner itself carries no untested computation. This is the standing lesson
      from `TemperingPlots.wls`/`TemperingMixing.wls`, which the `.wlt` suite still cannot reach.
- [ ] Results to a durable on-disk format that a re-run appends to rather than regenerates — the
      campaign will be run more than once.

### A4. Measurements
- [ ] Run the grid; record provenance (package version, seed, date) alongside the numbers.

### A5. Figures and draft integration
- [ ] Figures against the existing draft.
- [x] **`FacetLabeledCountExamples.md`: leave untracked.** Decided 2026-09-04. It stays a local
      working document; do not re-raise this in a later session.

---

## Track B — P2 foundations (weeks 1–4, overlapping A)

### B1. Test `H2dCombManifold` / `delH2dCombManifold`

The Hamiltonian the research runs on, at **zero tests**. Non-negotiable before it carries a paper.

- [ ] Differential test: for random graphs and every single-edge toggle,
      `delH2dCombManifold == H2dCombManifold[new] - H2dCombManifold[old]` within tolerance.
- [ ] **Include the vacuity assert.** A differential test that compares nothing passes silently;
      this project has been bitten twice. Assert that the comparison actually ran and that the
      energies are numeric.
- [ ] Fixed-value tests on hand-checkable states (the V=12 torus reference).

### B2. Land the homogeneous form

- [ ] `HomogeneousHamiltonian.wls` is **written and verified but lives outside the repo.** Land it,
      with tests. `HomogeneousHamiltonian.md` is already the spec.
- [ ] This is the source-term mechanism (`H = c·O`, `bt = 1`, `beta = c₀`) and therefore P2's
      critical path. It is *also* the natural API for Track D's general Hamiltonian (§D3), so it
      pays for itself twice.

### B3. The three source operators — **the main technical risk in P2**

Each needs an operator *and* an incremental delta, because the Metropolis inner loop consumes `delH`,
not `H`. Recomputing the operator per proposal would make the phase map unaffordable.

- [ ] **Triangle count.** `T = Tr(A³)/6`. Toggling {a,b} changes it by `±(A²)_{ab}` — the number of
      common neighbours. Exact, O(n), trivial.
- [ ] **Boundary edges.** Toggling {a,b} changes the triangle-membership count of {a,b} itself and of
      every edge {a,c}, {b,c} with c a common neighbour, so the delta touches O(deg) edges.
      Tractable; `CombinatorialBoundary` and `DGraphBoundary` already exist to check against.
- [ ] **Euler characteristic — the hard one, and the item to scope first.** χ of the clique complex
      is an alternating sum over *all* clique sizes. Toggling {a,b} changes the number of k-cliques
      by ±(the number of (k−2)-cliques inside `N(a) ∩ N(b)`), so

      `δχ = Σ_k (−1)^(k+1) · #{(k−2)-cliques in N(a) ∩ N(b)}`

      i.e. the delta is the (signed) clique-count vector of the induced subgraph on the common
      neighbourhood. That is far cheaper than recomputing χ, but it is not O(1) and its cost scales
      with the common-neighbourhood size — which is exactly where the model spends its time in the
      dense phase. **Benchmark this before committing to the phase-map grid.**
      Note `HWeightedFaceCounts` counts cliques only to size 5, so it cannot supply full χ as-is.
- [ ] Spec: `SourceTerms.md`, with the three derivations and their verification.

### B4. Companion scripts under test

- [ ] `TemperingPlots.wls` and `TemperingMixing.wls` produce P2's figures and the `.wlt` suite does
      not touch them (`ROADMAP.md` Phase 6). Either have the suite `Get` them, or accept hand
      verification — but decide deliberately, before they carry a paper's plots.

### B5. Phase map
- [ ] Sweep the three source couplings; locate the manifold-emergence transition.
- [ ] Finite-size scaling **across N up to 21** (§0.8) for the continuum limit. That is 5–6 usable
      sizes if the sequence is spread (e.g. 9, 12, 15, 18, 21), which is enough to fit a scaling
      exponent but leaves little slack — so the small-N end must be clean rather than merely
      cheap.
- [ ] Note for B3's cost: `N = 21` gives 210 possible edges, and the Euler-characteristic delta's
      expense scales with the common-neighbourhood size. **In the manifold phase the graph is
      sparse** (a 2-manifold triangulation has ~3N edges, and a manifold edge has exactly two
      common neighbours), so the delta is cheap exactly where the physics is. The cost lands in the
      disordered high-temperature phase, which matters for equilibration and for WL's traversal of
      the full energy range — not for the measurements themselves.

### B6. Tempering defects from Phase 6

Diagnosed on real runs, never acted on; each will distort a phase map if it is still live:

- [ ] Field CTL bootstrap's swap graph is disconnected (3 components, {8,4,4}) — round trips are
      impossible by construction, so a zero count means a severed ladder, not slow sampling.
- [ ] That run's schedule put 7 of 16 replicas outside the sampled fan; `HomogeneousHamiltonian.md`
      §8 has the inscribed-box recipe and the schedule was never rebuilt with it.
- [ ] Beta CTL bootstrap has a 42× acceptance spread, bottlenecked at the hot end.

---

## Track C — Wang–Landau (weeks 4–7)

Purpose: decide whether the manifold-emergence transition is **first or second order**. This is
precisely the question parallel tempering answers badly and flat-histogram sampling answers well.

### C0. What the labeled decision buys

P2 is labeled (§0.7), which removes the hardest part of this track before it starts:

- `g(E)` is a **plain density of states**, not a Γ-weighted one. Every standard WL formula applies
  as written.
- **No automorphism computation in the MC inner loop.** The unlabeled path calls
  `GroupOrder[GraphAutomorphismGroup[…]]` per proposal — the call whose own source comment flags a
  memory leak. Labeled runs never touch it.
- Validation against exact enumeration (§C3) counts labeled graphs, which is the straightforward
  count.

One thing to keep in view rather than act on: *Dynamical Quantum Multigraphs* found labeled graph
systems showing **no** thermodynamic phase transition for both the free and Ising Hamiltonians,
while the unlabeled ones did. `H2dCombManifold` is a different kind of object — nonlocal
connectivity penalties, manifold ground states rather than `K_N` — so that finding does not
transfer, but a referee will ask. Worth a sentence in P2 saying why the labeled ensemble is the
right one here. A secondary consideration if the manifold phase looks marginal: highly symmetric
states (a torus triangulation has a large `Γ`) have multiplicity `N!/|Γ|` among labeled states, so
the labeled ensemble entropically suppresses exactly the ordered states P2 is hunting — by an O(1)
factor, not an exponential one, so this is a caveat rather than a problem.

### C1. Spec — `WangLandau.md`
- [ ] Energy binning. With `J` fixed and integer-weighted source terms the spectrum may be discrete
      enough to use exact levels; otherwise bin width is a real parameter and must be justified.
- [ ] Modification-factor schedule (`f → √f` with a flatness criterion) versus **1/t-WL**, which
      avoids the error saturation the naive schedule suffers. Recommend 1/t.
- [ ] Move set: single edge toggle, matching `GraphSweepReplica`, so the existing driver is reusable.
- [ ] **Labeled or unlabeled** — see Open Question 1. If unlabeled, WL accumulates a Γ-weighted
      density of states, not a plain one, and every downstream formula changes. Settle before coding.

### C2. Implementation
- [ ] `GraphWangLandau[...]` returning `g(E)` on the energy grid plus the usual diagnostics.
- [ ] Reuse the existing `delH` convention and the log-space arithmetic already in the package
      (the log-space Metropolis work from 1.9.1 applies directly — `g(E)` spans many orders of
      magnitude and must be accumulated as `log g`).

### C3. Validation — *the reason to build this in ECGrav rather than ad hoc*
- [ ] Compare `g(E)` against **exact enumeration** at small N. The package can already enumerate
      the state space, so this is an exact test, not a statistical one — the same standard the
      MCMC generators were held to.
- [ ] Cross-check thermodynamics reconstructed from `g(E)` against the existing multi-histogram /
      MBAR path on a case where both are valid. Two independent routes to the same curve.

### C4. Order of the transition
- [ ] From `g(E)`, form `P(E) ∝ g(E)e^{−βE}` at the transition and look for a **double peak**.
- [ ] First order: peak separation ΔE grows with system size (latent heat) and the histogram
      minimum deepens with the interface free energy. Second order: single peak broadening, and the
      Binder energy cumulant minimum → 2/3.
- [ ] Report against several N — one N cannot distinguish the two.

### C5. Integration
- [ ] Feed `g(E)` into the existing reweighting path so surfaces come out of the same pipeline as
      the tempering runs.

---

## Track D — Quantum multigraphs (after P1 and P2)

From *Dynamical Quantum Multigraphs* (arXiv:2509.08296). The shipped package **is** the D = 2, k = 2
corner: `HIsing[Am,J,L] = (J/2)Σ_{i≠j}A²_{ij} + (L/2)Σ_{i≠j}A_{ij}` is the paper's angle-counting
Ising Hamiltonian, and `GraphSweepReplica`'s `UnlabeledVerticesYes` flag is the |Γ(G′)|/|Γ(G)|
acceptance ratio. The expansion is the other three quadrants.

### D0. Data structure — **decided**

The state of a k-uniform (p = k) complex over n vertices with D single-facet levels is an
association from level to facet set:

```
s = <| 0 -> X_0, 1 -> X_1, …, D-1 -> X_{D-1} |>
```

with each `X_j ⊆ C([n],k)` and `(X_0, …, X_{D-1})` a **weak ordered partition** of `C([n],k)` —
empty blocks allowed, union is everything, blocks disjoint. Example, D = 3, k = 3, n = 4:
`<|0 -> {}, 1 -> {}, 2 -> {{1,2,3},{1,2,4},{1,3,4},{2,3,4}}|>`.

**Adopted.** It is a direct transcription of the paper's occupation graph basis `|G_0,…,G_{D-1}⟩`,
so code will read like the formalism; it generalizes uniformly in D and k with no special cases; and
at D = 2 each block *is* an ordinary facet list, so **every existing complex observable
(`EulerChi`, `FVector`, `ConnectedComplexComponents`, …) applies to `s[j]` unchanged.** That last
point is worth a lot — it means the whole `PureComplexes.wl` observable library comes along for free.

Three amendments:

- [ ] **Move set — use the union variant; it removes a Hastings correction.** As described (pick an
      ordered pair `(a,b)`, pick a facet uniformly from `X_a`, move it to `X_b`) the proposal is
      *not* symmetric: the forward probability carries `1/|X_a|` and the reverse `1/(|X_b|+1)`, so
      acceptance must include a factor

      `|X_a| / (|X_b| + 1)`

      on top of the Boltzmann ratio. Omitting it biases the chain toward whichever level is
      emptier, and it will not show up as an obvious failure — it will show up as a wrong
      distribution. **Instead: pick an unordered pair `{a,b}`, pick a facet uniformly from
      `X_a ∪ X_b`, and move it to whichever side it is not on.** The union is invariant under the
      move, so the forward and reverse selection probabilities are identically
      `1/C(D,2) · 1/(|X_a|+|X_b|)`, they cancel, and the acceptance is pure Metropolis
      `min{1, e^{-βΔE}}`. It also handles empty blocks gracefully and reduces at D = 2 to exactly
      the single-facet toggle the existing drivers already use.
- [ ] **Keep the association as the canonical form; add a packed level vector as the working
      representation.** The association answers "give me a random facet at level j" in O(1), which
      the move needs; it answers "what level is facet f at?" only by search, which ΔE needs for every
      neighbouring facet. Carry alongside it an integer vector `ℓ` of length `C(n,k)` with
      `ℓ[[RankComb[f]]] ∈ {0,…,D-1}` — O(1) both ways, O(1) to update. **`RankComb` / `UnrankComb`
      are already public in the package** and are exactly this index. Converting between the two
      representations per MC step would dominate the inner loop, so build both once and update both
      incrementally. This mirrors what the package already does, where public Hamiltonians take a
      packed `Am_List` while complexes live as facet lists elsewhere.
- [ ] **Add a validator.** The weak-partition invariant is not enforced by the type — an association
      can hold overlapping or incomplete blocks. Follow the existing `PureComplexQ` / `DGraphQ`
      convention with a predicate that checks disjointness, coverage, and that every element is a
      k-subset of `[n]`.

*Note, so it is not "corrected" later:* this move set deliberately does **not** mirror the ladder
operators `L^±`, which shift a facet one level at a time. For a Hamiltonian diagonal in the
occupation basis, any ergodic reversible chain samples the same Gibbs measure, and letting a facet
jump between arbitrary levels mixes strictly better than a nearest-level walk. The operator algebra
is for expressing observables, not for constraining the sampler.

### D1. D = 2, k = 3 — the immediate need
- [ ] Quantum 3-uniform hypergraphs: the state is a set of triples, which is exactly a pure
      2-complex's facet set. `PureComplexes.wl` already supplies the combinatorial substrate, and
      P1's campaign will have exercised it hard.

### D2. The antisymmetric sector — settle before writing two samplers

Derived and verified this session (exhaustively over all of S_N for N ≤ 7, sampled at N = 8, for
k = 2, 3, 4): the phase an antisymmetric k-relation ket picks up under σ is

**χ(σ) = sgn(σ)^C(N−2, k−2)**

- At **k = 2** the exponent is 1, so χ ≡ sgn. The symmetric and antisymmetric sectors have identical
  automorphism groups and identical thermodynamics. This is why the paper can enumerate unlabeled
  quantum graphs as plain unlabeled graphs, and it is correct.
- At **k ≥ 3 this can fail.** When C(N−2, k−2) is even, χ ≡ 1, the stabilizer sum collapses to
  `Σ_{σ∈Γ(G)} sgn(σ)`, and that vanishes whenever Aut(G) contains an odd permutation — those states
  are annihilated, `𝒜|G⟩ = 0`. For k = 3 the exponent is N−2, so the A-sector matches the S-sector
  at odd N and shows exclusion at even N.

- [ ] Have the author check the derivation against his own conventions.
- [ ] **Consequence for the code: k = 2 needs one path; k ≥ 3 needs two, selected by the parity of
      C(N−2, k−2).** Discovering this mid-sampler would be expensive.
- [ ] If it holds, it is an exclusion principle for antisymmetric quantum uniform hypergraphs that
      is invisible at k = 2 — plausibly a section of the follow-up paper, not just an implementation
      note.

### D3. Operators and the Hamiltonian builder
- [ ] Ladder `L±_{ij}`, indicator `I^k_{ij}`, occupation `N^k`, subgraph occupation `N^k_g`.
- [ ] The general diagonal Hamiltonian `H = Σ_g f(g) Σ_k E^k_g N^k_g`. **This is the same object as
      the homogeneous form `H = c·O`** from B2 — a Hamiltonian built from a table of (subgraph,
      level, energy) weights. Landing B2 early makes this mostly assembly.
- [ ] The incremental `delH` for a subgraph-weight Hamiltonian is the real work, and it is the same
      local-subgraph-count problem as B3's Euler characteristic delta.

### D4. Aut for edge-coloured states
- [ ] For D > 2, Γ(|G⟩) = ⋂_k Γ(G_k) — the automorphism group of an edge-coloured graph. The
      current code calls `GroupOrder[GraphAutomorphismGroup[...]]`, which **its own source comment
      flags as leaking memory.** `ComplexAutomorphismIncidenceGraph` already showed how to get Aut
      from an auxiliary graph with no group enumeration (197× on pathological states) — the same
      technique applies here. **Carry over the warning that the gadget can silently overcount:** a
      colour-class construction must be checked against brute force before it is trusted.

### D5. D > 2, k > 2
- [ ] The full quadrant, once D1–D4 are in place.

---

## Open questions

Resolved 2026-09-04: the labeled/unlabeled choice (§0.7, §C0), the data structure (§D0), the
`FacetLabeledCountExamples.md` call (§A5), and the N range (§0.8, §B5). Remaining:

1. **Which graph represents a complex** for `SpectralDim` / `HausdorffDim` — 1-skeleton, dual
   graph, or `KFaceDistanceMatrix` (§A1). Still open because it is a *modelling* choice about what
   P1 means by "dimension", not a coding one. It must be fixed before any A-track figure, and
   §A1 now carries the trap that makes getting it wrong silent.
2. **Does P1 report one dimension notion or compare several?** Falls out of (1): if the three
   candidates disagree interestingly on random complexes, that comparison may be a result rather
   than a methodological footnote.
3. **What is the energy spectrum of `H2dCombManifold` + source terms** — discrete levels or
   effectively continuous? Decides whether WL bins or enumerates (§C1). Answerable by inspection
   once B3 lands, not by discussion.

## Risks

| Risk | Track | Mitigation |
|---|---|---|
| Euler-characteristic delta too slow in the dense phase | B3 | Benchmark before fixing the phase-map grid; fall back to a truncated clique sum if χ's tail is negligible |
| P1's figures rest on untested observables | A1 | Hardening is day one, not a follow-up |
| Campaign run 10–100× slower than it needs to be | A0/A3 | Stay on the warm path — counter before sampler, n ascending, many draws per cell (§1.2). Cold sampling costs seconds per draw and rises steeply with n |
| WL never converges on a rugged landscape | C | 1/t-WL; fall back to replica-exchange WL, reusing the existing tempering ladder |
| Companion scripts carry paper figures untested | B4 | Decide deliberately before the figures are made |
| Installed paclet lags the repo | all | Below 1.11.0 the private MBAR helpers do not exist, and an undefined private symbol returns *unevaluated, silently*. Check first whenever something quietly does nothing |
