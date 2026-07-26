# ECGrav test suite

Automated regression tests for the ECGrav paclet, derived from the manual
checks in `BuildingAndInstallingPaclets.nb`.

## Layout

| File | Contents |
| --- | --- |
| `TestPrelude.wl` | Loads the package **from source**, defines shared fixtures and float-comparison helpers. |
| `PureComplexes.wlt` | 54 tests: complex constructions, observables, dimensions, geometric predicates, and structural checks on the random-complex generators. |
| `MCSims.wlt` | 33 tests: exact assertions for Hamiltonians / data & free-energy helpers, smoke tests for the stochastic MC drivers and ground-state search. |

## Running

Via the Wolfram MCP `TestReport` tool, or from a kernel:

```wl
Needs["MUnit`"];
TestReport["/Users/012759760/Desktop/Research/ECGravMathematicaPackage/Paclet/ECGrav/Tests"]
```

Expected: **87 tests, 87 passing** (~22 s, dominated by the MC smoke tests).

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

**Two classes of test.** Hamiltonians and data-analysis helpers are pure
functions of their inputs, so exact values are asserted. The Monte Carlo
drivers are stochastic, so they are smoke-tested: the result is well-formed,
and the diagnostic-plot functions still emit their plots (this guards the
plot-restoration fix — these are the tests that would have caught the plots
being stripped).

**Oracles.** Where possible, expected values are independent of the
implementation: known surface topology (`EulerChi[tetrahedron] == 2`,
`EulerChi[torus] == 0`, octahedron Betti numbers `{1, 0, 1}`, f-vector
`{6, 12, 8}`, orientability of the seven surfaces) and values stated in the
notebook (`NumPureComplexes[3, 3] == 2649`, automorphism order 24). The
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

Roughly 58 of the package's public symbols. Recently added:

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

Deliberately **still untested**, each for a concrete reason:

- **Temperature schedulers** (`GraphCTLSchedule`, `GraphCEITempSchedule`) — the CTL/CEI
  algorithms run MC + reweighting at every schedule point and take minutes even on a tiny
  input, too slow for the suite.
- **Parallel uniform samplers** (`RandomUniformUnlabeledPureSimplicialComplex`,
  `RandomUniformFacetLabeledPureSimplicialComplex`, `RandomVertexLabeledPureSimplicialComplex`)
  — verified interactively, but the first call in a fresh kernel pays a ~30 s one-time
  parallel-distribution cost, not worth a structural smoke test.
- **Non-reproducible parallel functions** (`ErrorBootstrap`, the parallel
  `RandomPureSimplicialComplexMCMC` pipeline) — no `SeedRandom` reproducibility, so no
  stable assertion.
