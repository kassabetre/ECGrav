# Porting MCSims to Quantum Multigraphs (D > 2, k > 2)

Staged build plan for **Track D** of `ExpansionRoadmap.md`, started **2026-09-21** at `9df4e87`
(released state: 1.16.0). This document plans one thing only: bringing every Monte-Carlo routine in
`Kernel/Subpackages/MCSims.wl` from the D = 2, k = 2 corner to general (D, k). The mathematics of
the model stays in `ExpansionRoadmap.md` §D0–D5 and in the spec this plan produces
(`QuantumMultigraphs.md`, Stage 1).

**Sequencing note.** `ExpansionRoadmap.md` placed Track D *after* P1 and P2. Starting it in parallel
is a deliberate override, and it has one hard consequence: **§B2 (the homogeneous Hamiltonian) stops
being a P2 convenience and becomes a blocking prerequisite**, because §D3's general diagonal
Hamiltonian is the same object. `HomogeneousHamiltonian.wls` is written and verified but is **not in
the repo** — confirmed by inspection, nothing in `Kernel/` mentions it. See §5.

---

## 0. Decisions taken

Settled 2026-09-21; recorded so they are not relitigated.

1. **One code path, refactored in place.** The shipped `Graph*` drivers are repointed at a *state
   backend*; they keep their exact signatures and behaviour. Quantum entry points pass a different
   backend to the same drivers. The rejected alternative was a parallel `QMCSims.wl` stack, which
   would have meant reimplementing ~7 000 lines of equilibration, schedule and swap logic and
   letting two paths drift.
2. **The refactor must be provably behaviour-preserving, not plausibly so.** The adjacency backend
   reproduces the current proposal *exactly*, including its non-uniformity (§4.3), so a seeded
   golden-trajectory harness can assert bit-identity. That harness is Stage 0 and is built before
   anything is touched.
3. **`mag` generalises to mean level.** Today `mag = Total[Flatten[Am]]/(v(v−1))`, which is the edge
   density — i.e. exactly the occupation fraction of level 1. The generalisation
   `Σ_j j·|X_j| / C(n,k)` reduces to it identically at D = 2, k = 2 and **stays a scalar**, so the
   replica record keeps its shape and the 153 `Key["graph"]` / `Key["mag"]` sites need no
   restructuring.
4. **General (D, k) from the start.** The §D0 data structure has no special cases, so there is no
   reason to build D = 2, k = 3 first and generalise later. D = 2, k = 2 is the *validation corner*;
   D = 2, k = 3 is the first physics milestone (§0.6 of the roadmap).
5. **The antisymmetric sector (§D2) is out of scope for this plan.** The roadmap's first D2 item —
   have the author check the derivation of χ(σ) = sgn(σ)^C(N−2,k−2) against his own conventions — is
   still open. The samplers built here are the symmetric sector plus the optional Aut weight. Adding
   the A-sector later is a change to *one* backend primitive (`logSelectionProb`), not to the
   drivers.
6. **Labeled by default.** Matching §0.7, `UnlabeledVerticesYes = 0` is the path that must work
   first. The edge-coloured Aut machinery (§D4) is Stage 7 and is not on the critical path to a
   running quantum chain. §6 gives the measured reason: unlabeled costs ~2 500× labeled per step.

### Added 2026-09-21, from the design session

7. **The move is facet-first.** *Pick a hyperedge uniformly from all `C(n,k)` slots, then pick its
   new level uniformly from the `D−1` levels it is not on.* This is **not** a departure from §D0's
   union move — it is that move with the level pair drawn in proportion to `|X_a| + |X_b|` rather
   than uniformly, which collapses to picking the facet first. Both give transition probability
   `1/((D−1)·C(n,k))`; they are the same kernel. See §4.2 for why the correction matters.

   Verified exhaustively over all 3⁶ = 729 states at D = 3, k = 2, n = 4 against a structureless
   random energy: proposal symmetry `q(s→s′) − q(s′→s) = 0` **exactly**, and detailed balance
   against the Gibbs measure to 3 × 10⁻¹⁶. No Hastings factor, as §D0 requires.

8. **The working representation is the colouring; the association is the canonical form.** The state
   is a function `ℓ : C([n],k) → {0,…,D−1}`, stored as an integer vector of length `C(n,k)`, plus a
   length-D vector of block counts. The §D0 association is its fibre decomposition, materialised by a
   `GroupBy` **once per measurement, not per step**. This supersedes §D0's second amendment, which
   asked for both structures maintained incrementally: the facet-first move never needs "a random
   facet at level j", so the association is never required inside the loop. See §4.1.

9. **`RankComb` / `UnrankComb` define the canonical order; they do not implement it.**
   `Subsets[Range[n],{k}]` emits exactly `RankComb`'s order — verified over every slot at all twelve
   combinations of n ∈ {6,8,10,12}, k ∈ {2,3,4}. So both directions are lookups. The package
   functions stay in the **test suite as the independent oracle** (§4.4).

10. **The level count cannot be called `D` in code.** `D` is `Protected` — it is the derivative
    operator, and `Set::wrsym` is what you get. The spec and the paper keep `D`; the code uses
    `nLev`, fixed from Stage 1 before it is threaded through forty overloads.

---

## 1. The inventory — what actually needs porting

The decisive finding of the survey: **the port is not 9 214 lines of work.** MCSims splits cleanly,
and the split is verifiable rather than impressionistic — everything above `HIsing`
(`Kernel/Subpackages/MCSims.wl:966`) contains no occurrence of `Amat`, `seedGraph`,
`AdjacencyGraph`, or `Subsets[Range[...]]`.

| Layer | Lines | Public symbols | Port work |
|---|---|---|---|
| Settings (`$ECGravMax…`, `$ECGravCorrelationRunMultiplier`) | — | 3 | **none** |
| Statistics & MBAR reweighting (`LogSumExp` … `ConstrainedProbConjugateField`) | 45–949 | 13 | **none** — operates on energy/observable *numbers* |
| `ExactExpectationValue` | 1963–2090 | 1 | **none** — maps `Hamiltonian[#]` over an `ensemble_List`; a list of Associations matches |
| Hamiltonians + deltas | 966–1240 | 10 | new quantum members (Stage 4) |
| Optimisers (`GradDescent`, `SGradDescent`, `SimulatedAnnealing`, `LowEnergyStates`) | 1246–1960 | 4 | backend (Stage 5) |
| MC drivers | 2093–9214 | 8 symbols, 36 overloads | backend (Stages 3, 6) |

**17 of the 39 public MCSims symbols need no work at all.** That includes the entire MBAR stack —
`ComputeMinusBetaTimesFreeEnergy`, `MBARWeightBasis`, `NegativeBetaTimesFreeEnergy`,
`ExtrapolatedExpectationValue`, `InternalEnergy`, `CvOverT`, `MBAREffectiveSampleSize`,
`ConstrainedProbConjugateField` — so reweighting, free energies, conjugate-field densities and the
CTL metric all come along for free the moment a quantum chart exists. Three of
`GraphCTLSchedule`'s nine overloads (`MCSims.wl:6096`, `:6376`, `:6937`) are likewise pure
post-processing on `minusbetaF` and take no state, so **33 of 36 driver overloads** are the real
surface.

### 1.1 Where the driver overloads are

| Symbol | Overloads | State-bound |
|---|---|---|
| `GraphMetropolis` | 1 | 1 |
| `GraphSweepReplica` | 2 | 2 |
| `GraphEquilibriate` | 2 | 2 |
| `GraphComputeCorrelationTime` | 4 | 4 |
| `GraphMultiHistogram` | 6 | 6 |
| `GraphCEITempSchedule` | 4 | 4 |
| `GraphCTLSchedule` | 9 | 6 |
| `GraphParallelTempering` | 8 | 8 |

---

## 2. The backend interface

Every driver reaches the state through exactly nine operations. Enumerated from the inner loops of
`GraphSweepReplica` (`MCSims.wl:2190`, `:2363`) and `GraphEquilibriate` (`:2612`, `:2793`), which
between them contain every distinct state access in the file.

| Primitive | Adjacency backend (unchanged behaviour) | Multigraph backend |
|---|---|---|
| `size` | `vCount = Length[seedGraph]` | `C(n,k)` slots, plus `n`, `k`, `nLev` |
| `seedEmpty` | `Table[0,{vCount},{vCount}]` | `ConstantArray[0, C(n,k)]` |
| `seedRandom` | random symmetric 0/1 matrix | `RandomInteger[{0,nLev−1}, C(n,k)]` |
| `propose` | `row=RandomInteger[{1,v−1}]; col=RandomInteger[{row+1,v}]` | `i = RandomInteger[{1,C(n,k)}]`, then a level ≠ `ℓ[[i]]` (§0.7) |
| `deltaE` | `delH[s, params, row, col]` | `delH[s, params, i, from, to]` — **arity break, §4.5** |
| `apply` | `Mod[Am[[r,c]]+1, 2]` on both triangles | `ℓ[[i]] = to`; two ±1s on the block counts |
| `magDelta` | `(2/(v(v−1)))·(2·Am[[r,c]]−1)` | `(to − from)/C(n,k)` |
| `logSelectionProb` | `Log[N[|Γ(G′)|/|Γ(G)|]]` or `0` | coloured-incidence Aut (Stage 7); `0` while labeled |
| `canonicalise` | `ChooseNonIsomorphicGraphs` for `minStates` | `ℓ` equality while labeled (§4.1); Stage 7 otherwise |

The record the drivers pass around — `<|"graph"→s, "energy"→e, "mag"→m|>` — is unchanged in shape.
Only `"graph"` changes type, and decision §0.3 keeps `"mag"` scalar so that nothing downstream has to
learn a new shape.

**`HamiltonianUsableQ` / `DelHUsableQ` (`MCSims.wl:2559`, `:2563`) become backend-supplied.**
`DelHUsableQ` currently probes with `dh[amat, dhpar, 1, 2]` — two move coordinates, hard-coded. The
quantum probe needs three. This gate is exactly the one that caused the 1.8.1 bug, where a check
added across many call sites rejected correct Hamiltonians; the lesson recorded then was that such a
gate needs **a test per calling convention**, which here means one per backend.

---

## 3. Stages

Each stage ends in something runnable and asserted. Sizes are relative (S/M/L), not calendar.

### Stage 0 — The safety net *(S, blocks everything)*

The 61 MCSims tests are, by their own header, mostly *smoke tests* for the drivers: they assert a
well-formed return, not a particular chain. A refactor could silently change what the chain does and
they would stay green. Fix that before touching anything.

- [ ] `Tests/QuantumMCSims.wlt` — CI auto-discovers `*.wlt` (`Tests/ci-run.wls:66`), so a new file
      needs no workflow change.
- [ ] Golden-trajectory harness: `SeedRandom` a fixed seed, run each state-bound driver overload,
      record the energy/mag chain, assert bit-identity against stored reference data.
- [ ] **Establish first whether the parallel drivers are reproducible under a fixed seed.**
      `ParallelTable` is not deterministic across kernel counts in general. If they are not, the
      parallel overloads get a distributional assertion (two-sample on the energy histogram) instead
      of bit-identity, and that distinction is recorded per overload rather than assumed.
- **Exit:** every state-bound overload has a chain-level assertion; suite green.

### Stage 1 — The state type and its spec *(M)*

Implements `ExpansionRoadmap.md` §D0 as amended by decisions §0.7–§0.10.

**The representation, concretely.** At n = 4, k = 3 there are C(4,3) = 4 slots in `Subsets` order —
`{1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}`. The state

```
<|0 -> {{1,2,3},{1,2,4}}, 1 -> {{1,3,4}}, 2 -> {{2,3,4}}|>      association (canonical)
```

reads off as `ℓ = {0, 0, 1, 2}` with block counts `{2, 1, 1}`. Moving `{2,3,4}` from level 2 to
level 0 is `ℓ[[4]] : 2 → 0`, giving `ℓ′ = {0, 0, 1, 0}` and counts `{3, 1, 0}` — one integer write
and two counter bumps, and the empty block needs no special handling. Converting back reproduces
`<|0 -> {{1,2,3},{1,2,4},{2,3,4}}, 1 -> {{1,3,4}}, 2 -> {}|>` exactly.

At D = 2, k = 2 the colouring **is** the shipped representation: `ℓ` is literally the upper triangle
of the adjacency matrix read in `Subsets` order — verified identical. Nothing new is being adopted;
the matrix was always a two-colour colouring, and the colour count is becoming a parameter.

- [ ] The two lookup tables, built once per (n, k):
      `slots = Subsets[Range[n],{k}]` (position → facet) and
      `slotOf = AssociationThread[slots -> Range[S]]` (facet → position).
- [ ] `QMultigraphStateQ` — validates the weak-**ordered**-partition invariant: disjoint blocks,
      union is all of `C([n],k)`, every element a `k`-subset of `[n]`. Follows the existing
      `PureComplexQ` / `DGraphQ` convention.
- [ ] Constructors: ground (`ConstantArray[0,S]`), uniform random, from a list of facet sets.
- [ ] **The D = 2 bridge** — `AdjacencyMatrixFromQState` / `QStateFromAdjacencyMatrix`, exact and
      round-tripping. This is the oracle every later stage validates against.
- [ ] **Precompute the neighbour stencil**, which is what makes the inner loop conversion-free
      (§4.4). It is static in (n, k), never in the state.
- [ ] `QuantumMultigraphs.md` — the spec, opened here and grown through Stage 6, per the repo
      convention that new work ships with one. It states the order via `RankComb`; the code does not.
- **Tests:** association ↔ colouring round trip; validator negatives (overlapping blocks, missing
  facet, wrong arity, out-of-range vertex); and the oracle check
  `slots[[i]] === UnrankComb[i−1,n,k]+1` and `slotOf[f] === RankComb[f−1,n]+1` over the **whole**
  slot set — the assertion that catches a ±1 shift instead of letting it produce a plausible
  permutation of the levels.
- **Exit:** state algebra tested; the D = 2 bridge proven exact.

### Stage 2 — The backend records *(M)*

- [ ] Both backends as associations of the nine §2 primitives, plus the `…UsableQ` probes.
- [ ] Adjacency backend reproduces current behaviour *exactly*, non-uniform proposal included.
- [ ] Multigraph backend implements the facet-first move (§0.7).
- **Tests:** unit tests per primitive; the proposal-symmetry and detailed-balance assertions of §0.7,
  run exhaustively on a small state space rather than sampled.
- **Exit:** backends exist and are tested; **no driver has been changed yet**.

### Stage 3 — Repoint the chain core *(M)*

`GraphMetropolis` (1), `GraphSweepReplica` (2), `GraphEquilibriate` (2). Everything else calls these.

- [ ] Repoint the five overloads at the backend.
- [ ] **Stage 0's golden chains must stay bit-identical.** This is the gate; nothing proceeds past a
      red one.
- [ ] Add quantum entry points. Dispatch needs care: `GraphParallelTempering` already has
      `inputReplicas_Association` overloads (`MCSims.wl:7676`, `:7859`), so a bare
      `state_Association` would be ambiguous. Use a wrapper head or a `QMultigraphStateQ` condition.
- **Exit:** a quantum chain runs at D = 2, k = 2 and reproduces the shipped equilibrium energy
  histogram. Note this is a *distributional* comparison, not a trajectory one — the facet-first move
  is uniform over slots where the shipped move is not (§4.3), so the two chains sample the same
  Gibbs measure by different routes.

### Stage 4 — Hamiltonians and deltas *(L)*

- [ ] **Land §B2 first.** `HomogeneousHamiltonian.wls` → paclet, with tests. It is the prerequisite,
      it is not in the repo, and §D3 is mostly assembly on top of it.
- [ ] The general diagonal Hamiltonian `H = Σ_g f(g) Σ_k E^k_g N^k_g` — a table of
      (subgraph, level, energy) weights, which is the homogeneous `H = c·O` form.
- [ ] The incremental `delH`. **This is the real work of the stage**, and it is the same
      local-subgraph-count problem as §B3's Euler-characteristic delta: moving one facet between
      levels changes only the subgraph counts that contain it, so ΔE is a gather over the
      precomputed stencil — `ell[[nbr[[i]]]]`, measured at 0.5 μs.
- [ ] **Benchmark the stencil width before fixing any phase-map grid.** At n = 21, k = 3 a
      Hamiltonian coupling facets that share k−1 vertices has 54 neighbours; one coupling all facets
      that share *any* vertex has 513. That is the cost question for this stage, not the indexing.
- [ ] Operators as observables: ladder `L±_{ij}`, indicator `I^k_{ij}`, occupation `N^k`, subgraph
      occupation `N^k_g`.
- **Tests:** the delta identity, in the shape the suite already uses for
  `delHIsing-matches-HIsing-difference` — `delH == H[after] − H[before]` over random states and
  random moves, at several (n, k, D). Plus the reduction: at D = 2, k = 2 the quantum Hamiltonian
  must equal `HIsing` exactly on the bridged state.
- **Exit:** an H/delH pair whose delta identity holds and whose D = 2 corner reproduces `HIsing`.

### Stage 5 — Correlation time and the optimisers *(M)*

`GraphComputeCorrelationTime` (4), `GradDescent`, `SGradDescent`, `SimulatedAnnealing`,
`LowEnergyStates`.

- [ ] Repoint at the backend. `LazyCorrelationTime` (`MCSims.wl:3013`) is already numeric and needs
      nothing.
- [ ] `LowEnergyStates` and `SweepReplica`'s `minStates` dedup call `ChooseNonIsomorphicGraphs`.
      While labeled (§0.6) `canonicalise` is the identity and dedup is exact equality on `ℓ`, which
      is canonical for free (§4.1); the real version is Stage 7.
- [ ] The five `edgeList = Subsets[Range[nn],{2}]` sites in the optimisers become
      `backend["moveSet"]`.
- **Exit:** corrT measured on a quantum chain. The frozen-chain signature — a linear ACF falling to
  zero then flat, meaning `corrT ≈ s/2` for a single jump at sweep `s` — applies unchanged and
  should be expected at cold β on a large slot count.

### Stage 6 — Tempering and schedules *(L)*

`GraphParallelTempering` (8), `GraphMultiHistogram` (6), `GraphCTLSchedule` (6 state-bound),
`GraphCEITempSchedule` (4). The largest overload count, but the *least* state-coupled work: these
manage replicas, β ladders, swaps and chart layouts, and touch the state almost only through
`Key["graph"]` plumbing.

- [ ] Repoint; `ChooseRandomIndependentEdgeSet` is about the *replica swap graph*, not the state, and
      needs nothing.
- [ ] **Energy must stay at chart column 3** in every layout. That is the only column the package
      reads out of a β chart, and `TemperingPlots.wls` / `TemperingMixing.wls` depend on it.
- [ ] Confirm the MBAR stack consumes a quantum chart unchanged — it should, since it only ever sees
      `chart[[All,All,3]]`.
- **Exit:** parallel tempering on quantum states, with MBAR reweighting and a CTL schedule built from
  the result.

### Stage 7 — Aut and the unlabeled sector *(M)*

Implements §D4. Only needed for `UnlabeledVerticesYes = 1`. **The construction below was built and
validated on 2026-09-21; the first attempt at it was wrong, which is the stage's main lesson.**

**Do not intersect groups.** `Γ(|G⟩) = ⋂_j Aut(X_j)` is by definition the subgroup of `S_n`
preserving the *colouring*, so there is one automorphism group to compute, not D groups plus an
intersection — and intersecting permutation groups has no known polynomial algorithm, while one
coloured-graph Aut is what canonical labelling does. Verified on the n = 4 example: the intersection
route and the colour-preserving route return the identical group `{e, (34)}`, order 2, while the
individual `|Aut(X_j)|` are 4, 6 and 6 and are never needed.

**No group elements are enumerated anywhere.** `GroupOrder` on a permutation group given by
generators runs Schreier–Sims — a stabiliser chain, `|G|` as a product of basic orbit lengths.
`GroupOrder[SymmetricGroup[60]]` returns an 82-digit number in 1.7 μs. Finding the generators is
canonical labelling, which is backtrack search with partition refinement, also enumeration-free.

- [ ] **The validated encoding.** `Options[GraphAutomorphismGroup]` is `{}` — no colour support — and
      IGraphM is not installed, so the colouring goes into the graph's shape:
      a node per vertex, a node per k-subset, incidence edges between them, and **a pendant path of
      `ℓ(s)+1` nodes hanging off each facet-node.** Distinct levels give distinct tail lengths;
      a tail attached to a node of degree ≥ 3 cannot reverse and so contributes no factor.
      Validated **152/152** against brute force: 72 cases at n ∈ {5,6}, k ∈ {2,3}, D ∈ {2,3,4},
      plus 80 at the degenerate end (n,k) ∈ {(3,2),(4,2),(4,3),(5,4),(6,5)} where the vertex-node
      degree `C(n−1,k−1)` falls to 2 — plus the n = 4 example → 2 and uniform colourings →
      `6! = 720`, `7! = 5040` exactly.
- [ ] **Understand why it holds, because the guarantee is narrow.** In this construction **only
      facet-nodes carry tails**, so a vertex-node never carries a leaf and cannot be confused with a
      tail-interior node. Colouring *both* sides with tails does collapse — it failed **34 of 40**
      sparse pure complexes, exactly as `ComplexAutomorphismIncidenceGraph`'s own source comment
      predicts ("a factor of 2 per facet too large"). The quantum case is safe because every
      k-subset is present, so the incidence graph is never sparse; do not carry this construction
      over to P1's sparse complexes, which already have the package's rigid-triangle gadget.
- [ ] **Colours are mandatory, not an optimisation.** An uncoloured incidence graph admits
      side-swapping automorphisms: the triangle (n=3, k=2) gives 12 against a true 6, and the Fano
      plane (n=7, k=3) gives 336 against a true 168 — both self-dual, both 2× over.
- [ ] **Keep the brute-force comparison as a permanent test at n ≤ 7.** The obvious construction —
      one node per *level*, facets joined to their level-node, distinct pendant paths to stop levels
      swapping — **overcounts**: 7 of 72 cases wrong, always by 2×, and a uniform colouring at
      n = 6, k = 3 gives 2880 instead of 720. An *empty* level leaves its level-node dangling as a
      bare path endpoint, the component becomes a free-floating path, and the path reverses. Empty
      levels are common at cold β, which is exactly where the measurements are. This is §D4's
      "the gadget can silently overcount" warning, realised — a plausible rigidity argument was
      wrong and only brute force caught it.
- [ ] Replace the four `GroupOrder[GraphAutomorphismGroup[…]]` calls, whose own source comment flags
      a memory leak.
- [ ] **Isolated-vertex convention: resolved 2026-09-21.** `ComplexAutomorphismIncidenceGraph`
      (`PureComplexes.wl:2634`) builds its vertex set as `Union @@ facetsLst`, so the package
      **drops isolated vertices** — `Aut` is taken on the union of the facets, not on `[n]`. This
      never bites the quantum case, where every vertex sits in `C(n−1,k−1)` facets, but any P1
      reuse must account for it. Note the orders coincide on the n = 4 example, so no test built
      from that example would have caught a mismatch.
- **Exit:** unlabeled quantum chains run, with the Aut order verified against brute force.

### Stage 8 — Spec audit and release *(S)*

- [ ] Audit `QuantumMultigraphs.md` by **extracting its checkable identities and brute-forcing them
      against the definitions** — not by re-reading it for plausibility. That method found five
      defects in the pure-complex specs that no test caught; re-reading found none of them. The
      Stage 7 gadget failure is the same lesson in code rather than prose.
- [ ] `md2pdf.py` render, `CHANGELOG.md`, `PacletInfo.wl`, `README.md` (four version sites plus the
      tag bullet), `ReleaseNotes/vX.Y.Z.md`, `.gitignore` negation swap.
- **Exit:** released and asset-verified via the release API's `digest` field.

---

## 4. Traps carried forward

### 4.1 Do not carry the association as the working state

Two states that are equal can be different expressions: `<|0→{{1,2,3},{1,2,4}}|>` and
`<|0→{{1,2,4},{1,2,3}}|>` differ structurally but denote the same state. Anything using `===`,
`Union` or `DeleteDuplicates` — `minStates` dedup in `GraphSweepReplica`, `tracksAgreeQ` in
`GraphEquilibriate` — would silently count them as distinct and never raise anything. With `ℓ` as
the working representation, equality is structural and canonical for free. Materialise the
association per *measurement*, by `GroupBy` over `ℓ`, at O(C(n,k)) against C(n,k) proposals per
sweep — free.

### 4.2 Why the move was corrected, so it is not "simplified" back

The pair-first draw (pick `{a,b}` uniformly, then a facet from `X_a ∪ X_b`) is exactly symmetric and
needs no Hastings factor — that part of §D0 is right and was verified. But the probability of
proposing a *given* facet at level ℓ is `(1/C(D,2)) Σ_{m≠ℓ} 1/(|X_ℓ|+|X_m|)`, which depends on the
block sizes. Measured:

| occupancy | per-facet proposal rate, by level | max/min | pair draws hitting an empty union |
|---|---|---|---|
| D = 3, `(1, 1, 2000)` | 1.67e-1 / 1.67e-1 / **3.33e-4** | **501×** | 0 % |
| D = 3, `(70, 70, 70)` | uniform | 1.0× | 0 % |
| D = 10, `(5, 5, 200, 0×7)` | 3.34e-2 / 3.34e-2 / **9.95e-4** | 34× | **47 %** |

A facet in a nearly-empty block is proposed up to 500× more often than one in the bulk, so in a
condensed phase the chain shuffles two lonely facets between two sparse levels and starves the slow
mode. The last row is the one that matters for D > 2 generally: with most levels empty, nearly half
of all pair draws hit an empty union and do nothing. At D = 2 the two moves are **identical**, so
none of this touches the validation corner.

### 4.3 The current D = 2 proposal is non-uniform — and that is not a bug

`row = RandomInteger[{1,v−1}]; col = RandomInteger[{row+1,v}]` gives
`P({row,col}) = 1/((v−1)(v−row))`, so pairs with larger `row` are proposed more often. It is a proper
distribution and it is **independent of the state**, hence symmetric, hence Metropolis-valid. It
affects mixing and what "one sweep" covers, not the stationary distribution. **Do not tidy it during
the refactor** — that would break every golden chain for no correctness gain. Record it in the
adjacency backend with a comment saying why it is preserved.

### 4.4 `RankComb` is a specification device, not a runtime one

It is 0-based twice over — it takes subsets of `{0,…,n−1}` and returns a rank in `0 … C(n,k)−1` —
so the conversion is `RankComb[Sort[f]−1, n] + 1` going in and `UnrankComb[i−1,n,k] + 1` coming
back. Both shifts, both directions; get one wrong and you get a valid-looking permutation of the
levels rather than an error.

It is also slow, and it is never needed at runtime. At n = 21, k = 3:

| operation | cost |
|---|---|
| build both lookup tables (`Subsets` + `AssociationThread`) | **0.28 ms**, once |
| same tables via `RankComb` / `UnrankComb` over all slots | 32 ms / 48 ms |
| per lookup, position → facet (`Part`) | 0.074 μs |
| per lookup, facet → position (Association) | 0.25 μs |
| per call, `RankComb` | 24 μs — **96×** slower |

And the inner loop converts nothing at all, because the neighbour stencil is static in (n, k): build
`nbr` once (61 ms at n = 21, k = 3, 54 neighbours per facet), then ΔE is `ell[[nbr[[i]]]]` at 0.5 μs.
Ranking per step instead would cost 1.3 ms per ΔE and **1.73 s per sweep of ranking alone**. The
table can never itself be the bottleneck: it has C(n,k) entries, which is exactly one sweep's worth
of steps.

`RankComb`/`UnrankComb` earn their keep in the **test suite**, as an independently written oracle for
the table's order (§ Stage 1). 32 ms in a test file is nothing.

### 4.5 The `delH` arity break

Every driver calls `delH[state, params, row, col]`. The quantum move has three coordinates
(facet, from-level, to-level). This is the one interface change that cannot be hidden, which is why
it goes through the backend rather than through the drivers. It also propagates to `DelHUsableQ`'s
probe (§2) and to the test-file wrapper convention: `Tests/MCSims.wlt` documents at length why the
inert `dH[am, jj, ll, i, j]` wrapper **must** be a plain multi-argument definition and not a curried
one, since `hamiltonian_[hparams___]` binds a `Function` otherwise and the failure surfaces only
later, as an unevaluated energy delta. The quantum wrappers inherit that constraint with one more
argument.

### 4.6 `D` is `Protected`

It is the derivative operator; `nLev = 3` works, `D = 3` gives `Set::wrsym` and then silent
downstream nonsense — `Table[…,{lev,0,D-1}]` simply does not evaluate. Decided in §0.10, repeated
here because the spec, the roadmap and the paper all call it `D`.

### 4.7 PureComplexes observables have sharp edges

Quantum observables will want the complex library on `s[j]`, which §D0 correctly notes comes along
for free. Three known silent-wrong-answer paths apply:
`CountHoles[c,k]` returns `b_{k−1}`, not `b_k` (use the one-argument form);
`FractionInLargestComponent` misreads a *square* facet list as a weighted adjacency matrix (go
through `GraphFromCliques`); and `FacetAdjacencyMatrix` is not an adjacency matrix despite the name
(facet sizes on the diagonal, intersection sizes off it) — `SpectralDim` accepts it silently and
returns a plausible, meaningless number.

### 4.8 `SpecificHeat` / `Susceptibility` normalisation is a reporting convention

Both take `NN` = "number of sites" and divide (`MCSims.wl:120`, `:152`). For a k-uniform state the
site count is `C(n,k)`, not `n`. No code change — but fix the convention in the spec before any
figure is made, because a per-site `C_v` computed with the wrong `NN` is off by a constant factor
that nothing will flag.

### 4.9 Check the installed paclet version first

Below 1.11.0 the private MBAR helpers do not exist, and **an undefined private symbol returns
unevaluated, silently**. Whenever something quietly does nothing during this work, check
`PacletFind["ECGrav"]` before investigating anything else. This has cost a session before.

---

## 5. Running in parallel with the rest of Phase 7

### 5.1 Track interactions

1. **§B2 is now blocking, and it was not before.** §D3 assembles on the homogeneous form. Landing it
   serves P2's critical path *and* Stage 4, which is the strongest argument for doing it first.
   **Scope discovered 2026-09-21: B2 cannot land usefully without B3.** `HomogeneousHamiltonian.wls`
   does not exist on disk — the code lives only as listings in `HomogeneousHamiltonian.md` §4 and in
   `HomogeneousHamiltonianBug.nb`, as notebook scratch that captures globals (`e0`, `gN0`,
   `gNBdry0`) and calls five helpers. Of those, **only `EulerChi` is in the package**;
   `delEulerChi`, `numTriangles`, `delnumTriangles`, `numBdryEdges` and `delnumBdryEdges` are not.
   Those are §B3's three source operators, so B2's wrapper is a few lines and B3's operators are the
   substance. They land together or P2 gets nothing.
2. **Split the port by additive vs. mutating, not by stage number.** B2/B3, Stage 0, Stage 1 and
   Stage 2 all *add* symbols and change no existing behaviour, so they are safe on `main` and cannot
   disturb a running campaign. **Stage 3 is the first mutating stage** and is where a branch becomes
   necessary — by which point the harness that protects the refactor already exists.
3. **Wang–Landau is the backend's second customer.** §C2 needs the same nine primitives — propose,
   apply, ΔE, size — against a flat-histogram acceptance rule instead of Metropolis. Doing the
   extraction for Track D makes Track C substantially cheaper, and it means WL works on quantum
   states the day it is written rather than needing its own port later. Track B (§B3) and Track C
   both edit `MCSims.wl`, so landing the backend early means they build *on* it rather than against
   it.

Stages 1, 4 and 7 have no dependency on the other tracks and can be interleaved freely. Stages 3, 5
and 6 are strictly ordered after Stage 2.

### 5.2 Running P2's campaign alongside the port

P2's simulations run against the *installed* paclet while this work changes the *repo*, so the two
are already decoupled — but only until something is rebuilt. The rules that keep them decoupled:

- [ ] **Stage 0's harness is the contract between the two threads.** It proves the refactor leaves
      the chain bit-identical, which is exactly what licenses P2 to adopt a post-refactor release
      **without re-running anything**. Without it, any mid-campaign upgrade means re-running, or
      arguing in prose that results across two versions are comparable.
- [ ] **Never run `build.wls` or install from a development branch.** Releases are cut on `main`
      only. Then the paclet the simulations load cannot move while the refactor is in progress. This
      is free isolation and already matches the release mechanics; it needs to be a rule, not a
      habit.
- [ ] **Pin the campaign and make every run self-describing.** Assert the version at the top of each
      campaign notebook and stamp it into the saved chart beside the data:
      `Assert[First[PacletFind["ECGrav"]]["Version"] === "1.16.0"]`. A result that cannot name its
      own build is not reproducible, and this campaign spans a period in which the package changes.
- [ ] **Two paclets are installed as of 2026-09-21 — 1.16.0 and 1.5.0.** 1.16.0 resolves first, so
      this is latent rather than active, but 1.5.0 predates the entire MBAR layer, the chart tag
      column and `corrTMeasured`. If anything ever resolves to it, undefined private symbols return
      *unevaluated and silent*. Uninstall it.
- [ ] **Keep one real P2 chart as a test fixture.** §Stage 6 requires energy to stay at chart column
      3; make that provable rather than aspirational by asserting the post-refactor code still reads
      a stored campaign chart. 1.13.0 already moved chart observables once, and a campaign
      accumulates data over weeks.
- [ ] **Run the golden harness single-kernel** (`$KernelCount = 0`). Determinism probably requires it
      anyway (§7.1), and it keeps the suite off the cores the campaign is using — this machine has
      11.
- [ ] **Use two worktrees** — `git worktree add ../ECGrav-p2 main` — so the repo being consulted
      while debugging a long run is never in a mid-refactor state.

**The rule when they conflict: P2 figures win.** A refactor step that would invalidate a stored chart
or change a number already in a draft waits for a phase boundary. Adopt a new release into the
campaign only when it delivers something P2 needs — B2/B3's source terms is the one item on the list
that qualifies — and otherwise finish the campaign on a single version.

---

## 6. Sizing, and what it rules out

The slot count is `C(n,k)` and a sweep is one attempt per slot, so cost per sweep rises steeply in
k while the state space rises as `D^C(n,k)`.

| n | C(n,2) | C(n,3) | C(n,4) | states at D = 3, k = 3 |
|---|---|---|---|---|
| 10 | 45 | 120 | 210 | ~10^57 |
| 15 | 105 | 455 | 1 365 | ~10^217 |
| 18 | 153 | 816 | 3 060 | ~10^389 |
| 21 | 210 | 1 330 | 5 985 | ~10^634 |

At P2's ceiling of n = 21 a k = 3 sweep is ~6.3× a k = 2 sweep at the same n, and k = 4 is ~28×.
**Exhaustive enumeration is out beyond about n = 8, k = 3** — which matters because it is the only
route to an independent check of the Gibbs distribution. Plan the exact-enumeration cross-validation
at n ≤ 7 and treat everything above it as sampled.

**Unlabeled sampling is ~2 500× labeled, per step.** Measured at n = 12, k = 3: the Aut computation
itself is cheap (0.1–0.2 ms up to n = 14) and the graph *construction* dominates, so incremental tail
swapping brings a proposed move from 2.83 ms to 1.28 ms — against 0.5 μs for a labeled ΔE. Caching
`|Γ(s)|` for the current state saves nothing, since the ratio needs `|Γ(s′)|` for the proposal every
step. This is the measured justification for §0.6 and it means unlabeled statistics need a
deliberate design: smaller n, fewer measurements, or a colour-aware canonical labeller.

**The one lever with real headroom is IGraph/M**, a third-party Wolfram package wrapping the igraph
C library, which bundles the Bliss canonical-labelling engine. `IGBlissAutomorphismCount[{g, colours}]`
takes vertex colours natively and returns the count directly. That removes the gadget entirely and
shrinks the graph to `n + C(n,k)` nodes — **1 351 at n = 21, k = 3**, against ~4 011 for the tail
construction and 6 671 for the package's own triangle gadget (`n + 5F`). Better, the graph becomes
*static* in (n,k): per proposed move only the colour vector changes, which is where the 1.28 ms
actually goes (the Aut call is 0.1–0.2 ms; the rest is `EdgeAdd`/`EdgeDelete`).

**Licence caveat — this is why it cannot simply be adopted.** IGraph/M is **GPLv3**; ECGrav is MIT.
A paclet that `Needs["IGraphM`"]` at load is the combined-work case copyleft targets. It is also not
on either Wolfram paclet server (`PacletFindRemote["IGraphM"]` → `{}`), installing instead via
`Get["https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m"]`; it needs
Mathematica 11.0+, and this machine is on 15.0.1. **If adopted, it must be an OPTIONAL dependency**:
the validated gadget stays the default MIT path, Bliss is used only when present, users install it
themselves, and CI keeps testing the gadget path rather than adding a second flake source beside the
Wolfram Engine activation flake. Validate `IGBlissAutomorphismCount` against the same brute-force
suite before trusting it — a colour-convention mismatch would look exactly like the gadget bug did.

---

## 7. Open questions

1. **Are the parallel drivers seed-reproducible?** Decides whether Stage 0's harness asserts
   bit-identity or a distribution for roughly half the overloads. Answerable by experiment in the
   first session, not by discussion.
2. **What is the natural "empty" seed at D > 2?** `GraphEquilibriate` runs three tracks — user seed,
   empty, random — and uses their convergence as the equilibration test (`tracksAgreeQ`).
   All-at-level-0 is the obvious analogue of the empty graph, but at D > 2 it is no longer the
   *extreme* of anything; all-at-level-`nLev−1` is equally extreme. Whether the test needs a third
   bracketing track falls out of the first runs.
3. **Should IGraph/M be adopted as an optional Aut backend?** GPLv3 against ECGrav's MIT, so only
   ever optional (§6). Open pending the user's call; the gadget path is validated either way.
4. **Does `SpectralDim` on a level's 1-skeleton mean anything here?** Independent of this port, but
   it will be asked of the first quantum results. It is a *local* probe at a fixed step — identical
   on a cycle at n = 20, 60 and 150 — so it is not an asymptotic dimension, and saying so once here
   is cheaper than rediscovering it.

---

## 8. Risks

| Risk | Stage | Mitigation |
|---|---|---|
| Refactor silently changes released chain behaviour | 3, 5, 6 | Stage 0's golden harness is the gate; it is built before any edit |
| Parallel drivers not seed-reproducible, so half the harness is weak | 0 | Fall back to distributional assertions, recorded per overload rather than assumed |
| `delH` for subgraph weights too slow to be usable | 4 | The cost is the *stencil width* (54 vs 513 neighbours at n = 21, k = 3), not the indexing; benchmark it before committing to a grid, and share whatever locality trick §B3 finds |
| Aut gadget silently overcounts | 7 | Already happened once — the level-node construction fails on empty levels. Brute-force `S_n` check at n ≤ 7 stays a permanent test, not a one-off |
| Unlabeled runs too slow to give statistics | 6, 7 | 2 500× is measured, not feared: scope unlabeled work to small n, or install IGraphM and drop the gadget |
| Merge conflicts in `MCSims.wl` with Tracks B and C | all | Land Stages 0 and 2 first so the other tracks build on the backend |
| Wrong `NN` makes every per-site `C_v` off by a constant | 6, 8 | Pin the convention in the spec before the first figure |
