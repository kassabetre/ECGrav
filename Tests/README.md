# ECGrav test suite

Automated regression tests for the ECGrav paclet, derived from the manual
checks in `BuildingAndInstallingPaclets.nb`.

## Layout

| File | Contents |
| --- | --- |
| `TestPrelude.wl` | Loads the package **from source**, defines shared fixtures and float-comparison helpers. |
| `PureComplexes.wlt` | 38 deterministic assertions: complex constructions, observables, dimensions, geometric predicates. |
| `MCSims.wlt` | 20 tests: exact assertions for Hamiltonians / data helpers, smoke tests for the stochastic MC drivers. |

## Running

Via the Wolfram MCP `TestReport` tool, or from a kernel:

```wl
Needs["MUnit`"];
TestReport["/Users/012759760/Desktop/Research/ECGravMathematicaPackage/Paclet/ECGrav/Tests"]
```

Expected: **58 tests, 58 passing** (~20 s, dominated by the MC smoke tests).

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

## Known gaps / suspected bugs (not encoded as passing tests)

These were found while building the suite and are deliberately **not** asserted,
because doing so would bake in probably-wrong behaviour:

- **`CombinatorialSphereQ[torus]` returns `True`.** A torus (χ = 0) is not a
  sphere. Only the correct cases (tetrahedron, octahedron) are asserted.
- **`Branchedness[facetList]` leaks a `Null`**, returning e.g.
  `-1/4 + (-3/2 + Null)^2`. It works on `Graph` input (only that form is
  tested). It should either handle facet lists or reject them via `::argerr`.
- **`LowEnergyStates` divides by `$KernelCount`** and fails with `Power::infy`
  if called before `LaunchKernels[]`. Worked around in the prelude.
- **`ExactExpectationValue` emits `Part::partd`** while probing its observable
  argument. Harmless (result is correct) but declared as an expected message.

## Coverage

Not exhaustive: ~40 of the package's 116 public symbols. Untested areas include
parallel tempering (`GraphParallelTempering`, `GraphCTLSchedule`,
`GraphCEITempSchedule`), ground-state search (`SimulatedAnnealing`,
`GradDescent`, `SGradDescent`), the random-complex generators, and the
free-energy extrapolation family.
