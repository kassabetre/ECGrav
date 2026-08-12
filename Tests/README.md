# ECGrav test suite

Automated regression tests for the ECGrav paclet, derived from the manual
checks in `BuildingAndInstallingPaclets.nb`.

## Layout

| File | Contents |
| --- | --- |
| `TestPrelude.wl` | Loads the package **from source**, defines shared fixtures, float-comparison helpers, a message-argument collector and the correlation-time oracle. |
| `PureComplexes.wlt` | 135 tests: complex constructions, observables, dimensions, geometric predicates, structural checks on the random-complex generators, and the complex-space MCMC driver with its two preparatory stages. |
| `MCSims.wlt` | 34 tests: exact assertions for Hamiltonians / data & free-energy helpers, internal-consistency and smoke tests for the stochastic MC drivers, and ground-state search. |

## Running

Via the Wolfram MCP `TestReport` tool, or from a kernel:

```wl
Needs["MUnit`"];
TestReport["/Users/012759760/Desktop/Research/ECGravMathematicaPackage/Paclet/ECGrav/Tests"]
```

Expected: **169 tests, 169 passing** (~60 s; the MC smoke tests and the unlabeled-sampler uniformity check dominate).

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

**Stochastic is not the same as unreproducible.** The whole complex-space MCMC family —
`RandomPureSimplicialComplexMCMC` and its two stages — contains no `ParallelTable` in any of
its ten definitions, so two runs under the same `SeedRandom` agree exactly and results *can*
be pinned: `eqlT`, `converged`, `corrT` and `corrTValues` are all asserted as values. Check
for parallelism before falling back to a smoke test.

**Asserting messages.** `messageArguments["Sym::tag", expr]` (prelude) returns the filled-in
arguments of every matching message, not just the fact that one fired. It is built on
`Internal`HandlerBlock`, which sees messages *before* `Quiet` and before the `General::stop`
repeat limit — so the call under test can be wrapped in `Quiet`, keeping the run's output
clean, while the test still sees all of them, including the fourth and later copies of the
same message. `Check` and `$MessageList` give you neither the arguments nor the copies past
the limit. Every message assertion is paired with a negative control that the message does
*not* fire on the healthy path.

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

The three suspected-wrong-behaviour bugs previously flagged here have been fixed at the source
(see above) and are regression-tested. What follows came out of an audit of the equilibriation
criterion; the scale error it found has since been fixed, the rest has not.

All four equilibriators (`GraphEquilibriate` ×2 on three tracks, the two complex-space ones on
four) share one criterion: build `sqMeanEMat[[i,j]] = (mean of track i's newest 30 energies −
mean of track j's oldest 30)²`, and stop when its mean, `sqMeanPairwiseDiff`, falls below a
threshold. A second criterion stops when `Tr[sqMeanEMat]` is below `1e-5`.

**Fixed — the threshold was off by the effective sample size.** It used to be `meanLateVar`,
the variance of the *individual* energies, against a left-hand side that is a squared
difference of 30-point *means*. The criterion therefore reduced to a test of whether the
autocorrelation time is short: for equilibrated tracks the ratio it computed is `2·tau_int/15`,
measured across AR(1) tracks at `rho` 0 to 0.95 and matching that prediction, with the 50%
crossing at `tau_int ≈ 7.5` for both three and four tracks. At the `corrT` of 2 these chains
show, it passed with a 14× margin. The cost was power: two of four tracks displaced by 1σ were
declared converged 99% of the time, 50% detection needing ≈1.35σ. On the real chain,
`Partition[Range[20],2]` at `labelingChoise 2` correctly failed its first check at a 1.59σ
cross-track gap, then passed at sweep 430 with 1.04σ still outstanding, against 0.45σ once
genuinely mixed.

The threshold is now `expectedSqDiff`, twice the variance of a window mean, estimated from the
scatter of the non-overlapping window means inside each track. Its pass rate on equilibriated
tracks is flat in the autocorrelation at ≈50% per check. The residual gap on that same system
drops to **0.47σ** for ten extra sweeps. See the CHANGELOG for the reproducibility note.

Still open:

- **The second criterion can declare a stuck chain converged.** `Tr[sqMeanEMat]` is the
  diagonal only — each track's own drift — so four tracks that are each individually frozen at
  *different* energies give a trace of exactly 0 and pass, however far apart they are.
  Demonstrated at energies 1, 4, 9, 16: the first criterion correctly refuses, the second fires
  anyway. Not a knife edge — jitter of 1e-4 leaves the trace at 1e-9. It never fired in 4000
  real checks across five systems, so this is a latent path rather than an observed failure,
  and unlike the budget-exhausted case it reports `converged -> True` with no message. Its
  threshold also disagrees between overloads: `1e-5` in three equilibriators, `1e-6` in the
  fixed-vertex complex-space one.
- **The criterion is applied every sweep and exits on the first pass**, so it stops at the
  first favourable fluctuation rather than at the first sweep where the chain is genuinely
  mixed. Requiring it to hold for 30 consecutive checks was measured and takes the residual gap
  on that system further, 0.47σ → 0.21σ, but costs a great deal on systems that were already
  fine (`p2 M4 ch1`: sweep 401 → 875), so it was not adopted.
- `Abs` on `sqMeanPairwiseDiff` and on `Tr[sqMeanEMat]` is redundant — both are sums of squares.

`expectedSqDiff`'s predecessor was *not* the noisy part, contrary to the note this section used
to carry: its level varied by only 9–19% across checks on real runs, and consecutive checks
share 29 of their 30 points, so the median step was 2%. The scale error dominated it entirely.

Verified sound, so as not to leave these in doubt: the ring buffer is written at index 1 and
`RotateRight`ed each sweep, so at the first check (sweep 401) it holds the last 400 energies
with no initial zeros left, index 1 newest; `Entable[[i]][[1;;30]]` really is the recent window
and `[[371;;400]]` the old one. Both windows are 30 long. The criterion is scale-free in the
energy unit.

## The complex-space MCMC family

`RandomPureSimplicialComplexMCMC` promises independent uniform draws, and that promise rests
entirely on its two preparatory stages. The driver and both stages are covered in
`PureComplexes.wlt`, each in both its free-vertex and its fixed-vertex overload.

**Equilibriation.** `outWinLength = 400` and `inWinLength = 30`, so a converged run reports
`eqlT = numsweeps - 370`. The smallest systems still pass at the first check after the window
fills, which pins `eqlT` to 31 — do not read a repeated 31 as a stuck value. Anything larger now
takes a few more checks, so the pinned values vary (31 to 55 across the suite).

Failure is forced deterministically by setting `$ECGravMaxEquilibriationSweeps` **below 400**:
the convergence criterion is inside `If[numsweeps > outWinLength, …]`, so it never runs and
non-convergence does not depend on the draw. `eqlT` is then the lower bound
`Max[budget - 400 + 30, 1]`; budgets of 100 and 380 exercise the clamped and unclamped
branches. A *realistic* non-convergence — one where the criterion runs, fails, and puts real
numbers rather than `Indeterminate` into the `::noconv` diagnostics — needs a large system
whose replicas start far apart: `Partition[Range[20], 2]` with `labelingChoise -> 2` and a
budget of 402 fails on every seed tried, but costs ~19 s, so it is not in the suite. Nothing
smaller was reliable; at `{3,4}` and `{3,6}` the same setup converges on most seeds.

**Correlation time.** The reported values are checked against `mcmcCorrT` in the prelude, an
independent implementation of the documented integration rule. The driver measures its
observables internally, so the test recovers them: one of the supplied operators records the
complex it is handed — the same post-sweep complex the driver measures in that iteration —
and the energy and the vertex/edge count are deterministic functions of that complex, rebuilt
with the sweep's own arithmetic so the reconstruction is bit-for-bit.

Two things make that test non-vacuous. Real MCMC observables at small `(p, M)` decorrelate
within a handful of lags, so every value sits on the floor of 2 and an oracle comparison
proves almost nothing; a second operator that just counts sweeps supplies a monotone ramp
whose autocorrelation never turns over, which drives the lag loop to its end and produces a
value well off the floor. And the same oracle on shuffled series must give a *different*
answer, so the agreement is not an artifact of everything defaulting to 2.

**The end of the lag range needs a short run to be visible at all.** The integral runs to
`lastLag = numsweeps - 10`, includes its terminating non-positive term, and stops one lag
short when the autocorrelation never turns over — conventions carried over from the
pre-`f454c39` tabulated form. For a ramp the normalized autocorrelation at lag `t` is
`((N-t)/N)^2`, so over 310 sweeps the last ten lags come to about 0.02 in total and *no*
change to the lag budget can move `Ceiling`: the long-run test above is blind to it, and a
`lastLag -> numsweeps-20` mutation slips through it untouched. Over 20 sweeps those same ten
lags are half the integral, and both conventions become first-order.

That is what `-short-run-lag-range` is for. At `eqlT -> 2` the linear sweep counter reports 7
where `lastLag = numsweeps-20` reports 2, and a quadratic counter reports 6 where dropping
the one-lag-short convention reports 7 — the only place either is observable in the output.
Both mutations were run and both now fail the test. Note the probes have to be synthetic:
a real observable turns over within a handful of lags and never reaches the end of the range.

**The driver.** Four definitions: from a seed complex, free-vertex and fixed-vertex, and two
continuation overloads that take a previous run's association and skip straight to sampling.
Each returns `{measurements, data}`, a measurement row being
`{sweep index, energy, vertex or edge count, operator values…, complex}`.

The rows are checked for **internal consistency** rather than pinned: recomputing the energy,
the count and the operators from each row's own complex must give back that row's numbers
(`mcmcMeasurementsConsistentQ`, which recomputes through the package's own sweep at `NN = 0` —
no steps run, so what comes back is derived from the complex it was handed). A stale state or
an operator column out of step fails there where a shape check would not; a corrupted copy is
the negative control on the predicate.

Three behaviours needed a probe rather than an assertion on a healthy run:

- **The abort on unequilibriated chains** passes a call-counting operator. It has to stay at
  zero — `$Failed` on its own does not distinguish aborting from sampling and then giving up.
- **The continuation really skips both stages**, shown by the diagnostic plot: equilibriation
  is the only stage that `Print`s one, so the full pipeline emits at least one and the
  continuation must emit none.
- **`corrT` is the spacing between samples**, not just a reported number. Feeding the
  continuation `corrT -> 0` makes the sweep run no steps and freezes the chain, so every
  sample must come back identical to the state it started from; `corrT -> 6` is the control
  that the same call moves when allowed to. Zero is a probe, not a value the stages produce.

## Coverage

74 of the package's 117 public symbols are exercised by the suite. Recently added:

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
- **The complex-space MCMC family** — `RandomPureSimplicialComplexMCMCEquilibriate` and
  `RandomPureSimplicialComplexMCMCCorrelationTime` with pinned values, the `::noconv` /
  `::stuck` / `::alldefault` paths and their negative controls, and the correlation-time
  oracle; then `RandomPureSimplicialComplexMCMC` itself, all four definitions, by internal
  consistency of its measurement rows plus probes for the abort, the continuation and the
  sample spacing. See the section above.

Deliberately **still untested**, each for a concrete reason:

- **Temperature schedulers** (`GraphCTLSchedule`, `GraphCEITempSchedule`) — the CTL/CEI
  algorithms run MC + reweighting at every schedule point and take minutes even on a tiny
  input, too slow for the suite.
- **Non-reproducible parallel functions** (`ErrorBootstrap`, the graph MC drivers) — no
  `SeedRandom` reproducibility, so no stable assertion on values. Note this is not the same
  as untestable: where the output cross-checks itself, assert that instead (see the
  internal-consistency tests above).
- **Uniformity of the MCMC draws** — that the family samples uniformly is settled by exact
  enumeration (whole state space, orbit sizes against the per-class weights, detailed balance
  on the transition matrix in exact arithmetic), not by the suite. A goodness-of-fit test on
  the drawn samples would be both slow and weak here, since consecutive draws are correlated
  by construction; the suite checks that the machinery runs and reports honestly, and the
  distribution argument lives outside it.
