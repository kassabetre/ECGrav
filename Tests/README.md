# ECGrav test suite

Automated regression tests for the ECGrav paclet, derived from the manual
checks in `BuildingAndInstallingPaclets.nb`.

## Layout

| File | Contents |
| --- | --- |
| `TestPrelude.wl` | Loads the package **from source**, defines shared fixtures and float-comparison helpers. |
| `PureComplexes.wlt` | 115 tests: complex constructions, observables, dimensions, geometric predicates, and structural checks on the random-complex generators. |
| `MCSims.wlt` | 34 tests: exact assertions for Hamiltonians / data & free-energy helpers, internal-consistency and smoke tests for the stochastic MC drivers, and ground-state search. |

## Running

Via the Wolfram MCP `TestReport` tool, or from a kernel:

```wl
Needs["MUnit`"];
TestReport["/Users/012759760/Desktop/Research/ECGravMathematicaPackage/Paclet/ECGrav/Tests"]
```

Expected: **148 tests, 148 passing** (~60 s; the MC smoke tests and the unlabeled-sampler uniformity check dominate).

## Design notes

**Tests load from source, not from the installed paclet.** `TestPrelude.wl`
does `Get[.../Kernel/ECGrav.wl]` so the suite always exercises the current
working tree. Running the notebook instead tests whatever paclet is installed,
which can silently lag behind the source until `PacletBuild` is re-run.

**Paths are absolute.** Under `TestReport`, `$InputFileName` does not point at
the `.wlt` being run, so relative-path resolution is unreliable. Each file
prefers a `$InputFileName`-derived path only when it actually resolves, and
otherwise falls back to the absolute root in `TestPrelude.wl`. If the repo
moves, update `$ECGravRoot` there.

**`LaunchKernels[]` in the prelude.** Much of MCSims runs on `ParallelTable`,
and `LowEnergyStates` divides by `$KernelCount`, so it fails outright when no
parallel kernels exist. The notebook always calls `LaunchKernels[]` before use;
the prelude mirrors that.

**Float comparisons** use `SameTest -> floatEq` (1e-6) or `looseEq` (1e-3)
rather than exact equality.

**Three classes of test.** Hamiltonians and data-analysis helpers are pure
functions of their inputs, so exact values are asserted. The Monte Carlo
drivers are stochastic, so they are smoke-tested: the result is well-formed,
and the diagnostic-plot functions still emit their plots (this guards the
plot-restoration fix — these are the tests that would have caught the plots
being stripped).

Where a driver is stochastic but its output is *self-referential*, a third
option beats a smoke test: assert **internal consistency**. `GraphMultiHistogram`
and `GraphParallelTempering` are not seed-reproducible — two runs under the same
`SeedRandom` give different numbers — so no value can be pinned. But every
replica records both a graph and that graph's energy, so recomputing the
Hamiltonian from the recorded graph must reproduce the recorded energy, and the
recorded ground-state energy must be the energy of the recorded ground states.
Those hold for any correct run and fail on malformed output, which a shape check
does not. See `mcReplicasConsistentQ` and `groundStatesConsistentQ` in the
prelude, used by `GraphMultiHistogram-invariants` and
`GraphParallelTempering-invariants`.

One trap worth knowing: `minEstates` is legitimately empty when nothing beats the
`minEToBeat` threshold, and the ground-state check then passes vacuously. The
invariant test passes a reachable threshold (`0.0`, not the `-100.0` used by the
`Part::take` smoke test) precisely so the check has something to compare.

**Oracles.** Where possible, expected values are independent of the
implementation: known surface topology (`EulerChi[tetrahedron] == 2`,
`EulerChi[torus] == 0`, octahedron Betti numbers `{1, 0, 1}`, f-vector
`{6, 12, 8}`, orientability of the seven surfaces) and values stated in the
notebook (`NumVertexLabeledPureComplexes[3, 3] == 2649`, automorphism order 24). The
remainder are characterization tests that lock in current behaviour.

## Fixed bugs (found while building the suite, now regression-tested)

- **`Branchedness[facetList]` used to leak a `Null`** (e.g.
  `-1/4 + (-3/2 + Null)^2`) because a facet list was misinterpreted as an
  adjacency matrix. It now routes facet lists through `GraphFromCliques`, with
  the adjacency-matrix overloads guarded (mirroring `RmsPurity`). Guarded by
  `Branchedness-*` tests in `PureComplexes.wlt`.
- **`LowEnergyStates` divided by `$KernelCount`** and failed with `Power::infy`
  if called before `LaunchKernels[]`. Now uses `Max[$KernelCount, 1]`. The
  `LowEnergyStates-energies-no-kernels` test calls `CloseKernels[]` first to
  exercise the `$KernelCount == 0` path directly.
- **`CombinatorialSphereQ` conflated spheres with closed manifolds** — it returned
  `True` for a torus, Klein bottle, and projective plane. Split into
  `ClosedCombinatorialManifoldQ` (the closed-manifold test the old code actually
  implemented) and a corrected `CombinatorialSphereQ` (closed manifold **and** the Euler
  characteristic of a sphere, χ = 1 + (−1)^d; exact for surfaces). Guarded by
  `CombinatorialSphereQ-surfaces` and `ClosedCombinatorialManifoldQ-surfaces`.
- **`ExactExpectationValue` emitted a stray `Part::partd`** while probing its observable
  overload. The numeric-observable pattern test now uses `MatrixQ[obsValues, NumericQ]`.
  The `ExactExpectationValue-AvgDeg-highbeta` test now declares **no** expected messages,
  so it fails if the message reappears.
- **`GraphParallelTempering` / `GraphMultiHistogram` emitted `Part::take`** when a
  replica's accepted-swap count reached `NN`; the beta-history trim now clamps to the
  buffer length. A `Mod[_, 0]` progress-print glitch on short runs (`Floor[NN/5] → 0`)
  was fixed at the same time. Guarded by `GraphParallelTempering-smoke-no-Part-take`,
  which asserts a clean, message-free run at the triggering boundary (`NN == swapAccept`).

## Known gaps / suspected bugs

None outstanding — the three suspected-wrong-behaviour bugs previously flagged here have
been fixed at the source (see above) and are regression-tested.

## Coverage

70 of the package's 116 public symbols are exercised by the suite. Recently added:

- **Free-energy / reweighting helpers** — `EmpCorrelationTime`, `DSpecificHeat`,
  `NegativeBetaTimesFreeEnergy`, `CvOverT`, `ExtrapolatedExpectationValue`,
  `ComputeMinusBetaTimesFreeEnergy`, as deterministic oracles. The reweighting family is
  tested at the **single measured beta**, where the free-energy factor cancels and each
  reduces to a closed form (−βF equals the supplied value, the extrapolated observable
  equals its sample mean, Cv/T equals β²(⟨E²⟩−⟨E⟩²), the self-consistent free energy is 0).
- **Ground-state search** — `GradDescent` (deterministic: returns a valid adjacency matrix
  and never increases the energy), `SGradDescent` and `SimulatedAnnealing` (seeded smoke).
- **Random-complex generators** — structural checks on
  `RandomFacetLabeledPureSimplicialComplex`, `RandomUnlabeledPseudoManifold`, and
  `RandomPureSimplicialComplexMCMCSweep` (the output is a well-formed pure complex of the
  requested shape, whatever the draw).
- **Uniform samplers** — all three are covered, both branches. An earlier note excluded them
  for a "~30 s one-time parallel-distribution cost"; that cost was a **bug**, not an inherent
  price — the subkernels never received the counting functions, so every parallel sample came
  back malformed. Fixed in 1.5.0, and the parallel branch is now asserted directly by forcing
  the thresholds down. `RandomUniformUnlabeledPureSimplicialComplex` additionally has exact
  identity checks on its Burnside weights and profile marginals, which is what pins its
  uniformity; a wrong weight would still emit plausible complexes.

Deliberately **still untested**, each for a concrete reason:

- **Temperature schedulers** (`GraphCTLSchedule`, `GraphCEITempSchedule`) — the CTL/CEI
  algorithms run MC + reweighting at every schedule point and take minutes even on a tiny
  input, too slow for the suite.
- **Non-reproducible parallel functions** (`ErrorBootstrap`, the parallel
  `RandomPureSimplicialComplexMCMC` pipeline) — no `SeedRandom` reproducibility, so no
  stable assertion on values. Note this is not the same as untestable: where the output
  cross-checks itself, assert that instead (see the internal-consistency tests above).
