# Changelog

All notable changes to ECGrav are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning is
semantic-ish; breaking changes are called out explicitly.

## [Unreleased]

## [1.12.0] - 2026-08-24

### Changed
- **The MBAR free-energy solve is rebuilt on one shared matrix.** `ComputeMinusBetaTimesFreeEnergy`
  is the self-consistency iteration that *produces* the `minusBetaF` every other MBAR routine
  consumes, and it was the most expensive thing left in the package: `O(K^3 N)` interpreted work per
  iteration, iterated to convergence. 1.11.0 did not touch it.

  Writing `A[[s,j]]` for the reduced potential of sample `s` at state `j`, each pass is
  `L[s] = LogSumExp_j (logn_j - F_j - A[s,j])` then `F_k = LogSumExp_s (-L[s] - A[s,k])`. The term
  carrying the target is `A[s,k]`, which has no `j` in it, so it factors out of the inner
  `LogSumExp` and the outer loop over targets collapses into one column-wise reduction. `A` does not
  depend on `F`, so it is built once and reused across every iteration. All three overloads now
  build `A` and call a shared `MBARSolveFreeEnergies`.

  Measured on a real 11-beta, 999-sample `GraphCTLSchedule` output: **218.2 s on eleven kernels
  became 0.378 s on one.** The two external-field overloads agree to exactly `0.0` on the same
  one-dimensional data written as scalars and as one-element vectors.
- **The convergence test is replaced, and this changes results.** The old
  `Sqrt[Sum[(1 - F_k/Fnew_k)^2]]` had three defects. Free energies are defined only up to an
  additive constant, so dividing by `F_k` made the test depend on a gauge carrying no physical
  meaning — on one converged fixture it read `2.4e-15` in the mean-zero gauge, `2.9e-16` shifted by
  ten, and **`Indeterminate`** in the `F_1 = 0` gauge, where it divides by zero. It compared the
  already mean-centred `F` against an uncentred `Fnew`, so a gauge shift leaked in. And it stopped
  about **twelve times short** of the accuracy it appeared to promise.

  It is now `Max[Abs[Fnew - Mean[Fnew] - F]]` — gauge-invariant, no division, like compared with
  like — converged against a geometric extrapolation of the remaining error rather than the raw
  step, so the tolerance means what it says: asked for `1e-8` it returns `9.8e-9` where the raw step
  returns `1.2e-7`. The default tightens from `1e-5`/`1e-6` to **`1e-8`**, which the rewrite makes
  free: `1e-6` costs 0.27 s and `1e-12` costs 0.56 s.

  **Free energies therefore come out about four orders of magnitude more converged, and values
  computed by earlier releases are off by ~1e-4.** See the release notes before upgrading.

### Tests
- The fast form is pinned against the **literal nested-sum definition**, written out independently
  and iterated a fixed number of times rather than to a tolerance. Its first draft was vacuous — it
  compared the two field overloads at `bt = 1`, where deleting the `bt` factor is the identity, and
  the suite passed with that break in place. It now compares at `bt = 1.7` with explicit
  `bt`-sensitivity legs.

## [1.11.1] - 2026-08-23

### Fixed
- **`GraphParallelTempering` rejected correct hamiltonians when handed a schedule from
  `GraphCTLSchedule`.** Its two schedule-driven overloads (`MCSims.wl:7880` with `delH`, `8178`
  without) returned `$Failed` with `ECGrav::badham` before doing any work, so

  ```wl
  GraphParallelTempering[seed, bt, h[], dH[], sched[[2 ;; 3]], conjObs, obs, NN, unlabeled, minE]
  ```

  failed every time with the `h[am_List, f_Real]` form those drivers require.

  This is the 1.8.1 bug in two sites that fix did not reach. These drivers take the callback
  parameters from a replica's field **label** — `Apply[hamiltonian,
  replicaSwapEdgesLabels[[2, Key[i]]]]` — which is a *third* calling convention alongside
  `hparams` (beta tempering) and the field table (`GraphMultiHistogram`, `GraphCTLSchedule`).
  1.8.1 fixed the other two. With the hamiltonian passed as `h[]`, `hparams` is empty, so the
  probe evaluated `h[graph]`, which does not match `h[am_List, f_Real]` and comes back
  unevaluated. Both gates now probe with `Sequence@@replicaSwapEdgesLabels[[2,1]]`, matching how
  the continue-from-replicas overloads already read `inputReplicas[[1,"externalField"]]`.

  All 24 gate sites were audited against how their own driver invokes the callbacks; exactly these
  two disagreed.

### Tests
- **The callback gate is now pinned at the call site, not only at the probe.**
  `EnergyCallbackProbe-external-field-arity` checks that the probe behaves correctly *given* the
  right arguments — deliberately, since the drivers equilibrate before returning — and that is the
  gap this bug fell through: nothing checked that a driver *hands* it the right arguments. The new
  test stubs both probes so the gate fails on purpose and the driver returns before equilibrating,
  then asserts the field value reached the probe. Against the unfixed code all three of its legs
  return `False`.

## [1.11.0] - 2026-08-22

### Changed
- **The CTL thermodynamic metric is rebuilt from one shared MBAR weight vector.** Building it was
  the dominant cost of a multi-field `GraphCTLSchedule`, and it grew as the *square* of the
  bootstrap table: a few seconds with one external field, over seven hours with two. It is now a
  fraction of a second. Three reductions, all of the same observation — writing `D_s(x)` for the
  log MBAR denominator of sample `s` at target field `x`, every one of
  `NegativeBetaTimesFreeEnergy`, `ExtrapolatedExpectationValue` and `MBAREffectiveSampleSize` is
  built from that one vector.

  The free energy **cancels**: an expectation value is `(v.u)/Total[u]` with `u_s = Exp[-D_s]`, an
  ordinary weighted mean, so the free energy was only ever the normaliser and divides back out. It
  is no longer computed. The weight vector does not depend on which observable is averaged against
  it, so the six interpolation grids at two tempered components collapse to **one pass**. And
  `D_s(x) = beta x.O_s + Log[Sum_j Exp[Log[n_j] - F_j - beta h_j.O_s]]`, whose second term — the
  only part quadratic in the bootstrap table — carries no `x` and is therefore hoisted out of the
  grid entirely.

  Measured against 1.10.1, same fixture, same code path, serial: 24.18 s -> 0.072 s at K=16,N=120;
  767.50 s -> 1.215 s at K=100,N=100. **The schedules are unchanged**; this is an algebraic
  identity, verified to 2e-14 against the form it replaces, on targets inside and far outside the
  sampled box. Full derivation, cost tables and implementation map in `MBARWeights.md`.
- **The MBAR weights are formed in log space with a maximum shift.** `ExtrapolatedExpectationValue`
  divided by raw denominators and emitted `General::munfl` under extrapolation — which is the
  regime the interpolation box is built for. Every output is a ratio of the weights, so they can be
  rescaled freely.
- **The metric grid is filled serially.** The pass is now sub-second, so distributing the sample
  matrix to subkernels costs more than it saves. It was `ParallelTable`'s coarse-grained chunking
  that left the front end sitting on an unchanging "Distributing definitions" panel for the length
  of a build.
- **The grid-size message reports the real count.** It claimed `numInterpolatingSamples^numVars`
  = 400; the box is widened 20% each way at the unwidened step, so the grid is 29^numVars = 841.

### Removed
- **The effective-sample-size interpolant.** `rho` now evaluates the MBAR effective sample size
  exactly, which costs 0.107 ms and only became affordable with the shared weight basis. **No
  schedule changes** — `densityPoints` steps by the same increment as the interpolation grid and
  starts four steps inside it, so its points coincided with interpolation *nodes* and the
  interpolant was only ever asked for values it stored exactly (0 of 441 mask verdicts differ, at
  every bootstrap size tried). It goes because that correctness rested on an unrecorded coincidence
  between two grid definitions; evaluated off-node the interpolant disagrees about the mask at ~6%
  of points and returns negative values, meaningless for a count bounded below by 1.

### Documentation
- **`MBARWeights.md`** — the algebra above, the reference implementation, the cost and verification
  tables, and the implementation map.
- **`HomogeneousHamiltonian.md`** — how to temper in inverse temperature *and* an external field at
  once using only the shipped API, by writing the hamiltonian in homogeneous form so that beta
  becomes one of the tempered couplings. Includes the units caveat, which is the part most easily
  missed: recorded energies become `c.O = beta*H_phys`, meaningless to compare across replicas once
  beta varies.
- **`FacetLabeledCount.md`** — the specification of `NumFacetLabeledPureComplexes`, the last of the
  four pure-complex counters and samplers without one. Every numeric claim in it was checked
  against the shipped code.
- The pure-complex specification line anchors are refreshed, and `UnlabeledCount.md` §2.1 now gives
  a `p = 2` counterexample for the non-faithfulness of the `Aut` action, where the old text gave a
  `p = 3` one next to a claim that the kernel appears only from `p >= 3`.

## [1.10.1] - 2026-08-19

### Tests
- **The sweep's dependence on beta is pinned.** `GraphSweepReplica` reads `beta` in exactly one
  place — `logRatio = Log[selectionProb] - delE*beta` — and nothing else in the driver touches it.
  That is what makes a *homogeneous reparameterisation* sound: rewrite a hamiltonian so every
  coupling is explicit, `H = Sum_i c_i O_i`, run at `bt = 1` with `c` scaled by beta, and the chain
  is not merely similar but identical, because the code only ever forms the product `beta*delH`.

  This matters because the multi-field tempering path tempers the field vector `c` at *fixed* `bt`,
  so beta is not one of its coordinates. Promoting the beta-conjugate part of a hamiltonian to a
  conjugate observable with its own coupling makes it one — which is how a single ladder can
  temper in temperature and external field together — and that construction is only valid while
  this identity holds. If a second use of beta is ever added to the acceptance, the new test
  notices: injecting one makes it fail on both the trajectory and the energy-scaling row.

### Fixed
- **`GraphCTLSchedule` permuted the components of a replica's external-field vector.** The
  replica centers were ordered with `LexicographicSort[Sort/@centers, NumericalOrder]`, and the
  inner `Sort` sorted the components *inside* each center rather than ordering the list of
  centers. Each center is a field vector whose component `i` is the coupling of conjugate
  observable `i`, so permuting them hands the replica different physical parameters.

  With a single external field the inner `Sort` is the identity, which is why this survived —
  and the three multi-field overloads, which require `Length[externalFieldTable[[1]]] > 1`, were
  the only callers that could reach it. There it swapped the components of any center whose
  entries were not already ascending, **and only those**, so a schedule came back with some
  replicas correct and others carrying permuted couplings: the result is not even a permutation
  of the centers the clustering produced. On a hamiltonian written in the homogeneous form
  `c·O`, where the first component is beta and the rest are beta times the physical fields, it
  turns a positive beta into whichever coupling happened to be smallest.

  Single-field schedules are unaffected: the inner `Sort` is the identity on length-1 vectors, so
  no existing one-field result changes. Untested — `centers` is a `Module` local, so reaching it
  needs a full schedule run (an MC bootstrap, 400 MBAR extrapolations, 40,000 kernel-density
  draws), and the comment now carried at each of the three sites is the guard against the inner
  `Sort` coming back.
- **The CTL thermodynamic metric was built from second moments, not covariances.** `sigmaMat` was
  assembled from `ExtrapolatedExpectationValue[... O_i*O_j ...]` — and that function is a plain
  MBAR-reweighted expectation with no mean subtraction, despite the comment above `sigmaMat`
  reading `S_{a,b}=b^2Cov(M^a,M^b)`. The metric that places replicas is the Fisher information,
  which for an exponential family in these natural parameters is the covariance of the sufficient
  statistics; `E[O_i O_j]` overstates it by `E[O_i]E[O_j]`.

  With one external field the error is a smooth positive inflation of a 1×1 metric, so schedules
  still came out plausible — which is why this survived. With several fields, or with a hamiltonian
  written in homogeneous `c·O` form where one conjugate observable carries the whole energy, the
  mean term dominates and replica placement stops tracking fluctuations at all. The first moments
  are now interpolated over the same grid through the same MBAR call and subtracted.

  **This changes the schedules single-field runs produce.** The metric feeds both the sampling
  density (`rho = Sqrt[Det[sigmaMat]]^pSoft`) and the K-means distance function, so an inflated
  1×1 metric placed replicas differently — the ordering was sane but the spacing was not the
  constant-thermodynamic-length spacing it claimed to be. Any schedule generated before this fix
  should be regenerated; schedules are inputs to a tempering run, not results, so nothing measured
  becomes wrong, but the acceptance uniformity CTL exists to deliver was not being delivered.
- **`Max[0.0, …]` zeroed the metric's off-diagonal entries.** The same expression floored *every*
  entry at zero. That is correct for a variance and wrong for a covariance: two anti-correlated
  observables have a legitimately negative cross term, and clipping it discarded precisely the
  off-diagonal structure that makes one multi-field schedule better than one schedule per axis.
  Only the diagonal is floored now — a negative variance is MBAR noise and zero is the right
  floor — while off-diagonals keep their sign. Residual indefiniteness is absorbed downstream,
  where `Det` and the quadratic form in `metricDistance` are both already floored at zero.

  Never fired at one external field: a second moment `E[O^2]` cannot be negative.

  Both fixes live in the new private `CTLMetricFromMoments`, called from all three multi-field
  overloads. Unlike the ordering fix above this one is real arithmetic with a wrong answer to
  catch, so it is extracted and tested — against a population covariance computed directly from
  samples, on a fixture whose means are nonzero *and* whose cross second moment is negative, so a
  single case exercises both defects. Cost: the metric interpolation now also builds `numVars`
  first-moment interpolants alongside the `numVars(numVars+1)/2` second-moment ones.
- **`GraphCTLSchedule` placed replicas on an extrapolated metric.** The metric is interpolated
  over the bounding *box* of the bootstrap field table, but a box is not the region the bootstrap
  run actually covered. The corners are reached purely by extrapolation — the function really is
  named `ExtrapolatedExpectationValue` — and a metric invented there can attract replicas to
  parts of parameter space where nothing was ever sampled. The existing `pSoft` exponent, which
  "softens extreme spikes so one region doesn't eat all points", was treating the symptom.

  The replica density is now masked by the **Kish effective sample size** of the MBAR reweighting
  at each point: `rho` is zero wherever fewer than 1% of the samples effectively stand behind the
  estimate. Measured on a 12-field, 1200-sample synthetic reweighting, the effective count runs
  127–300 across the sampled range and drops to exactly 1.0 one step outside it, so the floor sits
  in a wide gap rather than cutting through the distribution.

  This is a hard restriction of the domain, not a reweighting. Inside the supported region the
  density is still exactly `Sqrt[Det[g]]^pSoft`, which is what makes the spacing
  constant-thermodynamic-length; weighting by effective sample size instead would have tilted
  schedules towards well-sampled regions and stopped being CTL.

  It is geometry-free, which is deliberate. The unsupported region is an unvisited box corner for
  ordinary multi-field tempering and a fan through the origin for a hamiltonian in homogeneous
  `c·O` form, and a support test handles both without either case needing its own grid shape.

  **Limitation, stated because it is easy to over-trust:** effective sample size measures weight
  *concentration*, not distributional gaps. A region can carry a respectable effective count and
  still be badly extrapolated if the true distribution there lives outside what was sampled. This
  is a necessary condition for trusting a point, not a sufficient one.

  The new private `MBAREffectiveSampleSize` computes it in log space with a maximum shift. The
  direct form does not fail under extrapolation — WL escapes to arbitrary precision rather than
  overflowing — but that escape is expensive: 34x slower one field-width outside the sampled range
  on a 60-sample fixture, and unfinished after ten minutes on a 720-sample one, against a log form
  flat in the target. Cost is one additional interpolant per schedule.

### Notes
- **`LogSumExp` shifts by the mean, not the maximum.** Mathematically identical, but the mean
  shift does not survive a wide spread of exponents — the largest term can still overflow. This is
  why `MBAREffectiveSampleSize` does its own maximum shift locally rather than calling it. Not
  changed here; flagged because the 1.9.1 work fixed exactly this pattern in `SGradDescent`'s
  softmax and did not reach `LogSumExp` itself.

## [1.10.0] - 2026-08-17

The correlation-time routines sized their measurement run as `Min[5*eqlT,
$ECGravMaxCorrelationSweeps]`. The cap was a setting; the multiplier was a literal. Since the
smaller of the two wins, and `5*eqlT` is the smaller for every `eqlT` below 200 at the default
cap of 1000, a user on a short run who raised the cap to lengthen it saw no change at all and had
no knob that worked. The multiplier is now a setting too.

No measured value moves: the default is the 5 it was hard-coded to.

### Added
- **`$ECGravCorrelationRunMultiplier`**, defaulting to 5, sets how many equilibriation times long
  a correlation-time measurement run is. The run length is now
  `Ceiling[Min[$ECGravCorrelationRunMultiplier*eqlT, $ECGravMaxCorrelationSweeps]]` in all six
  routines — the four `GraphComputeCorrelationTime` overloads and both
  `RandomPureSimplicialComplexMCMCCorrelationTime` overloads.

  It is the half of the pair that binds on short runs, which is exactly the case where a run is
  too short to resolve what it found. To tell which of the two is binding, read the run length the
  `::shortrun` message reports: equal to the cap means capped, anything else means the multiplier
  set it.

  Fractional values are allowed and the run length is rounded up. That is not cosmetic — the lag
  range is `numsweeps - 10` and `LazyCorrelationTime` takes `lastLag_Integer`, so a non-integral
  run length would leave it unevaluated and then be indexed as though it were an association.

  Settable by plain assignment, like the other two settings: it survives a package re-load, stays
  assignable after the blanket `Protect`, and `SyncParallelSettings` pushes it to the subkernels
  alongside `$ECGravMaxCorrelationSweeps` so a parallel driver does not silently fall back to the
  subkernel's own default.

### Changed
- **Both `::shortrun` messages ended in "raise the latter"**, meaning the cap — which was the only
  advice that could be given when the multiplier was a literal, and is now wrong whenever the
  multiplier is what set the run. They name both settings and tell you how to work out which one
  is binding. The multiplier is message slot 3, so it also now reports the value actually in
  effect rather than a hard-coded 5.

### Tests
- `GraphComputeCorrelationTime-run-length-multiplier` and
  `RandomPureSimplicialComplexMCMCCorrelationTime-run-length-multiplier` count the sweeps a run
  actually takes at raised, lowered, fractional and cap-bound multipliers. The complex-side test
  also checks the `::shortrun` argument, since the multiplier appears twice per routine — once
  setting the run, once in the message — and the two moving apart would make the message
  misreport the setting it is telling you to raise.
- `SyncParallelSettings-reaches-subkernels` now covers the new setting, checked through
  `OwnValues` on the subkernel rather than `ParallelEvaluate[sym]`, which inlines the master's
  value and would report success either way.

## [1.9.1] - 2026-08-16

Three fixes, all of them about a quantity leaving the range a machine number can hold.
`EulerChi[{}]` returned `$Failed` instead of 0, which affects you if any energy callback takes the
**link of an edge** — that slice is `{}` whenever the two endpoints share no neighbour. The
Metropolis acceptance test now compares in log space, so high-beta runs no longer flood with
`General::munfl`. And `SGradDescent`'s softmax is shifted by its maximum before exponentiating,
which is the one that was failing outright rather than merely complaining.

Bug fixes only. No public function was removed, no argument order changed, and no result that
previously computed changes value — the acceptance decisions and the softmax distribution are
both identical where the old code produced an answer at all, verified against it.

### Fixed
- **The empty complex now has an Euler characteristic of 0.** Both `List` overloads gate on
  `Depth == 3`, and `Depth[{}]` is 2, so `{}` fell through to the catch-all — even though
  `EulerChiAM`, which the adjacency overload delegates to, already had an `n == 0` branch
  returning 0. The convention was implemented and unreachable.

  This is not an exotic input. Any callback computing how a quantity changes when edge `(i,j)`
  is toggled evaluates `am[[sphInt, sphInt]]`, and that slice is `{}` exactly when `i` and `j`
  have no common neighbour. On sparse graphs that is most pairs — on 12 vertices with 18 edges,
  32 of the 66. A reported `GraphCTLSchedule` run died this way: the user's `delEulerChi`
  returned `-1 + $Failed`, which flowed into `dH` as the Metropolis energy delta. It survived a
  spot check because that used one hand-picked pair whose link happened to be non-empty, and it
  survived the `DelHUsableQ` gate because that probes `dH[seed, 1, 2]` — one pair, on the seed
  graph. Neither can see a failure that depends on which pair and which state.

  Matched as a literal rather than by widening either `Depth` guard, so the three existing
  overloads are untouched. `{}` reads as either an empty facet list or an empty adjacency
  matrix; both denote the empty complex and both give 0, so one definition serves.
- **`EulerChi::argerr` named only the facet-list form** while the function accepts three. It now
  names all of them, and the usage message documents the `{}` convention.
- **The Metropolis acceptance test is compared in log space**, at all 17 sites: `Log[u] < -beta*dE`
  rather than `u < Exp[-beta*dE]`. `Log` is monotonic, so the two decide identically — verified
  over 34.5M plain comparisons and 220k covering the unlabeled path's automorphism-ratio factor,
  with zero disagreements, and end-to-end against the previous code at beta 7.5 and 20, where
  energies, magnetizations and edge counts came back identical.

  The exponential is not representable across much of this model's range. At beta 7.5 a single
  edge toggle on a 12-vertex graph reaches `beta*dE` of 2340 (mean 1570) against a threshold of
  about 709, so most proposals underflowed; a run at beta 20 emitted 55 underflow messages where
  the log form emits none. The cost is real but was never a correctness problem — `Exp[-775.]`
  returns exactly `0.`, so the comparison rejected, which is the right answer. What it cost was
  about **19x** on those calls: the penalty tracks the result leaving the normalized range, not
  the size of the argument (`Exp[-700.]` is as fast as `Exp[-5.]`; `Exp[-730.]` is 19x slower).

  Three things this had to get right, none of them a search-and-replace:
  - `selectionProb` in the unlabeled path is an exact `Integer`/`Rational` ratio of automorphism
    group orders. `Log` of an exact ratio stays symbolic and drags symbolic arithmetic through
    every proposal at about 1.7x, so it is numericized first. The **product** also underflows in
    its own right when that ratio is below 1, so clamping `Exp` alone would not have sufficed.
  - The acceptance `If` has **three** branches; the fourth argument fires when the comparison is
    undetermined and it means reject. That is what silently swallowed the `EulerChi` `$Failed`
    described above, and it is preserved.
  - `SimulatedAnnealing` looked immune because it wrapped `bta` in `SetPrecision`. It was not:
    `N[machineReal, precision]` cannot *add* precision — only `SetPrecision` fabricates digits —
    so `delE` stayed at machine precision and `arbitrary*machine` collapsed back to machine. It
    underflowed like the rest. Raising both operands would have worked at about 5x the cost of a
    machine `Exp`, where the log form costs about 1x.
- **`SGradDescent`'s softmax over edges is shifted by its maximum before exponentiating.** This
  one was a genuine failure rather than noise. The edge to flip is drawn with probability
  proportional to `Exp[-beta*deltaE]`, and when every move is uphill — which any model whose
  single-edge energies run to the hundreds reaches at high beta — every weight underflowed to
  `0.`, `Total` came out `0.`, and `w/Total[w]` made every entry `Indeterminate`. `RandomChoice`
  then failed outright. Measured on a complete graph at beta 200, the old code raised ten
  distinct message types (`General::munfl`, `Power::infy`, `RandomChoice::wghtv`, `Part::partw`
  and more) and returned through a cascade of errors; the new code raises none and returns a
  proper distribution.

  Note the sign runs the opposite way from the acceptance test: the exponent is `-beta*deltaE`,
  so it is the **downhill** edges — the ones softmax exists to favour — that made it large and
  positive and pushed `Exp` into arbitrary precision at the other end.

  Subtracting the maximum is exact, since the common factor cancels from numerator and
  denominator; over 20,000 well-conditioned tables the two forms agree to 7.8e-16, and the
  benign smoke-test runs are unchanged. Entries more than 700 below the maximum are assigned the
  `0.` they round to rather than exponentiated, which is what they already were to machine
  accuracy and keeps the table off the underflow path entirely. Where the old form gave three
  `Indeterminate`s, the new one recovers the distribution those weights actually have —
  `1 : 1.9*^-98 : 2.7*^-261`.

### Added
- **`SGradDescent-softmax-survives-all-uphill`** in `Tests/MCSims.wlt`, asserting the run is
  message-free rather than merely non-crashing — the old failure announced itself loudly and a
  shape check passed straight through it. It fails against the previous code.
- **`GraphSweepReplica-deep-reject-at-extreme-beta`** in `Tests/MCSims.wlt`, pinning the regime
  where the old form stopped being representable: from the complete graph at beta 400 every
  proposal is uphill by `beta*dE = 3200`, and the chain must sit still and return exact numbers,
  on both the labeled and unlabeled paths.
- **`EulerChi-empty-complex-is-0`** in `Tests/PureComplexes.wlt`, asserting both the value and
  the concrete route that reaches it — the link-of-an-edge identity `chi(link) - 1`, on a graph
  whose vertices 1 and 2 share no neighbour. It fails against the previous code. Separately
  verified across 40 random 8-vertex graphs: 1120 toggled pairs, 302 of them with an empty link,
  and the incremental identity agrees with recomputing `chi(new) - chi(old)` in every one.
  188 tests, all passing.

## [1.9.0] - 2026-08-14

`GraphComputeCorrelationTime` was well behind its complex-space twin,
`RandomPureSimplicialComplexMCMCCorrelationTime`, which had been rewritten in 1.7.0. This brings
the four graph overloads up to it: the lazy lag walk, and the three defects the complex side had
already had fixed. Relevant if you raise `$ECGravMaxCorrelationSweeps` — the measurement was
quadratic in the run length, so a bigger cap cost as its square.

Every correlation time that previously computed is unchanged. Verified two ways: 147 synthetic
AR(1) series (rho 0 to 0.99, corrT up to 82) agree with the tabulated form on `corrT`, on the
integration cut-off, and on the plotted curve, term for term, with a negative control that
separates the conventions in 72 of them; and 96 end-to-end configurations across all four
overloads give 400 of 400 exact matches on every value the old code was able to produce.

### Changed
- **The lag table is walked, not tabulated.** The integral stops at the first non-positive
  autocorrelation, but all `numsweeps-10` lags were computed first and the tail discarded — about
  990 of them at the default cap, to use a handful. Each costs O(run length), so the measurement
  was quadratic. A new private `LazyCorrelationTime` walks the lags, stops at the turnover, and
  then extends the curve — without summing — over exactly the window the plot displays. Measured
  on the correlation step alone: 8-54x at 1000 sweeps, 20-155x at 5000, 70-581x at 20000, the
  range running from strongly autocorrelated to uncorrelated observables. Both conventions of the
  tabulated form are preserved deliberately: the terminating non-positive term is part of the
  integral, and an autocorrelation that never turns over stops one lag short.

  One case does not gain: an observable that never turns over at all — a drifting or monotone
  one — still walks the whole range, and the walk carries about 20-25% more interpreter overhead
  per lag than the vectorised `Table` it replaced. That is a constant factor, flat in the run
  length rather than growing with it (0.75x at 1000 sweeps, 0.81x at 10000), and the curve is
  filled by index rather than by `AppendTo`, which would otherwise have restored the quadratic
  cost precisely in that case.

### Fixed
- **`corrT` came back as an unevaluated expression instead of a number on short runs.** All four
  overloads read the integration cut-off from `FirstPosition[...,numsweeps-10][[1]]`, whose
  default was a bare integer rather than a position list. Whenever the autocorrelation never went
  non-positive — which every run with `numsweeps < 10` does, its lag table being empty — that
  indexed an integer, and `corrT` came back as a symbolic `Max[...]`. Callers then multiplied
  their sweep counts by it. `eqlT <= 1` reaches this on any input. The lazy walk removes the
  construct rather than patching it.
- **The overload with an operator list and no delH had no stuck-observable handling at all.** A
  frozen observable — a constant operator, and `Length[g]` is one, since the vertex count does not
  change during a sweep — put a `Null` in the autocorrelation table, which indexed its way into
  `corrT` and returned it unevaluated. It now filters frozen observables, reports them through
  `::stuck` and `::alldefault`, emits `::shortrun`, initialises `corrTValues` with one slot per
  observable, and draws the autocorrelation plot, matching the overload that takes a delH. Two
  driver paths reach it: the no-delH external-field overloads of `GraphMultiHistogram` and
  `GraphParallelTempering`.
- **`::stuck` reported the wrong value** in the operator-list overload with delH: the branch named
  `operators[[i-2]]` but printed `observablesTable[[i-2,1]]`, so a stuck operator was reported
  with the value of the observable two rows above it — for the first operator, the energy. Same
  fix the complex side received.
- **The autocorrelation plot's legend ran ahead of its curves.** Only fluctuating observables get
  a curve, but the legend was built from all of them, so any stuck observable shifted every label
  onto the wrong line.

### Added
- **`"corrTMeasured"` in the result of every correlation-time routine**, on both the graph and
  complex sides, and carried through `RandomPureSimplicialComplexMCMC`. It is `True` when `corrT`
  came from an autocorrelation integral and `False` when it is the floor of 2 that every
  correlation time defaults to — which happens when every observable was frozen, and when the run
  was too short to have a lag range at all. Neither was distinguishable from a chain that
  genuinely decorrelates in two sweeps by looking at `corrT` alone, and the second is reported by
  no message, so the flag was its only possible signal. Same remedy, for the same reason, as the
  `"converged"` key on the equilibriators.

  Deliberately a flag rather than a larger default. `corrT` is the sweep count the drivers hand
  to `GraphSweepReplica`, so cost is linear in it: inflating it for stuck runs would multiply the
  cost of every run that froze, and at high beta a frozen chain usually means the ground state was
  found rather than that anything went wrong. A single stuck observable alongside a live one does
  **not** clear the flag — `corrT` is the maximum over the observables that do fluctuate. What to
  do about an unmeasured correlation time is left to the caller. Additive: no existing value
  changes, re-verified against the same 96 end-to-end configurations (400/400).
- Seven tests in `Tests/MCSims.wlt` and one in `Tests/PureComplexes.wlt`, taking the complex
  side's as the template — the graph side had two, which is why all of this survived. They pin the
  integration rule against the prelude's independent oracle on both operator-list overloads (with
  a monotone ramp forcing the never-turns-over branch, and a shuffled-series vacuity guard), the
  `::stuck` value, the all-frozen `::alldefault` path, that short runs return a number, that both
  operator-list overloads emit their plot, and each of the four cases `corrTMeasured` separates.
  All six of the rewrite tests fail against the previous code. 185 tests, all passing.

## [1.8.1] - 2026-08-14

The energy-callback check added in 1.7.0 rejected **correct** hamiltonians in every external-field
driver, returning `$Failed` with `ECGrav::badham` before doing any work. If you run
`GraphMultiHistogram`, `GraphCTLSchedule` or `GraphParallelTempering` over a table of external
field values, this affects you on 1.7.0, 1.7.1 and 1.8.0.

Bug fix only. No public function was removed, no argument order changed, no new setting, and no
result that previously computed changes value.

### Fixed
- **The energy-callback probe used the wrong arity in every external-field driver.** In an
  external-field run the callbacks' parameters are a row of the field table, not `hparams`: the
  driver does `Apply[hamiltonian, externalFieldTable[[i]]]` and `hparams` is empty, so the
  callback is ultimately called as `h[graph, field]`. The probe nonetheless evaluated `h[graph]`,
  which for a hamiltonian written the way those drivers require — `h[am_List, f_Real]`, passed as
  `h[]` — comes back unevaluated and is not a number, so the driver rejected it. On a reported
  8-vertex `GraphCTLSchedule` over five `NT0` values, `h[seed, -0.5]` is `310.` while `h[seed]` is
  unevaluated, and the run died at the gate. Eight probe sites now pass a real field row: the
  `externalFieldTable` overloads of `GraphMultiHistogram` (2) and `GraphCTLSchedule` (4), and the
  two external-field `GraphParallelTempering` overloads that take replicas from a previous run,
  where the row comes from `inputReplicas[[1,"externalField"]]`. The beta-tempering overloads are
  untouched — there `hparams` really is the parameter list and their probes were always correct.

### Added
- **`EnergyCallbackProbe-external-field-arity`** in `Tests/MCSims.wlt`. There were no tests on
  this check at all, which is why the bug shipped. It asserts that the field-row arity is
  accepted, that the old `hparams` probe rejects it — so the test cannot pass vacuously — and that
  an undefined head is still rejected. It also records one limit rather than hiding it: a
  slot-based pure `Function` ignores extra arguments, so it satisfies this probe and is caught one
  level in, by the inner driver's own probe on the applied form. 177 tests, all passing.

### Notes
- The README's test count had gone stale at 148 and is corrected to 177.

## [1.8.0] - 2026-08-13

Performance release. `RandomUniformFacetLabeledPureSimplicialComplex` draws 2x to 103x faster.
Same distribution, same output, **no public function added, removed or renamed, no argument
order changed**, and the Wolfram Language 15.0+ requirement is unchanged — the draw was simply
doing work proportional to the number of candidate facets where the number of distinct *weights*
is a small constant.

### Changed
- **The facet-labeled draw groups candidates by completion key.** `RandFLPCOne` weighed every
  admissible facet separately: at `{4,6,15}` that is 441 candidates for each of 6 facets, and it
  was essentially the whole cost of a sample. But a candidate's weight enters only through
  `nu(uncovered \ S)`, a multiplicity vector of block sizes summing to at most `p`, so there are
  at most `Sum[PartitionsP[i],{i,0,p}]` distinct weights per facet — 4 at `p=2`, 7 at `p=3`, 12
  at `p=4` — however many hundreds of candidates exist. At `{4,6,15}` those 441 candidates carry
  3 distinct keys. The draw now picks a key group with weight (group size) x (completions) and
  then a member uniformly, which is the same distribution because members of a group have equal
  weight. Keys ride as base-`(p+1)` integer codes updated in place when a block is first covered,
  and the per-cycle-type candidate tables are memoised in the new private `RandFLPCCandTable`
  (released by `NumPCClearCache` with the rest). Measured warm against the previous body: 2.1x at
  `{2,4,6}`, 17.6x at `{3,6,12}`, 51.5x at `{4,6,15}`, 76.8x at `{3,12,24}`, **102.9x at
  `{4,8,20}`** (52.8 ms to 0.51 ms per sample), 87.5x at `{3,15,30}`.
- **`$RandomFacetLabeledParallelThreshold` raised from 100 to 2000.** What the parallel branch
  pays is fixed rather than per-sample — `DistributedContexts` ships the whole of
  `ECGrav`Private`, memo tables included, on every call — so cutting the serial cost by 20-100x
  moved the crossover out with it. At 100 samples the parallel branch had become a 3x to 16x
  *pessimisation*. Remeasured on 11 subkernels, the crossover now runs from 200 samples at
  `{2,3,4}` to past 3200 at `{4,8,20}`; the constant sits at the top of that range because the
  penalties are lopsided, too low costing up to 16x and too high costing at most the parallel
  speedup, which is only 2-3x anywhere near the crossover.

### Added
- **`FacetLabeledSampler.md`**, a specification of the sampler alongside `Theory.md`,
  `UnlabeledCount.md` and `UnlabeledSampler.md`: the pair-sampling principle, the lemma that
  licenses the sequential facet draw, the `O(n^p)` state indexing, the grouped draw, the
  identities that establish uniformity, a line-numbered implementation map, cost tables, and
  references. Linked from the README, which now points at all three sampler specifications.
- **`RandFLPCCandTable-keys-track-uncovered`**, pinning the two invariants the grouped draw rests
  on — that the candidate table lists exactly the admissible facets, and that the incrementally
  updated codes still decode to `nu(uncovered \ S)` after any sequence of facets — over every
  cycle type with `p <= 4`, `n <= 11`. Get either wrong and the sampler still emits perfectly
  plausible complexes from the wrong distribution, so the identities are checked directly rather
  than through the output. The test was confirmed to fail under an injected wrong packing radix
  and under skipped code updates before being trusted. 176 tests, all passing.

### Notes
- A rejection-free dynamic program over row-type multiplicities carrying the column-equality
  partition as DP state was evaluated as an alternative and rejected. It is correct — its counts
  agree with `NumFacetLabeledPureComplexes` on all 234 triples with `p <= 4`, `M <= 6` — but its
  state space is `2^M (n+1) (p+1)^M Bell(M)`, and since `n <= pM` the only free axes are `p` and
  `M`, both exponential. Its table costs 118 s and 1.8 GB at `{3,7,14}` and did not build in
  180 s / 4 GB at `{3,8,16}`, where this sampler now draws in 0.27 ms. `FacetLabeledSampler.md`
  section 11 records the measurements.

## [1.7.1] - 2026-08-13

The ground-state search reported minima that were not minima. A sweep weighed a state's energy
against its running minimum only for states it *stepped into* by an accepted move, so the state
it started from — the one the chain occupies through every rejected move, and the one a replica
swap hands it — was never a candidate. Anything that searches for ground states is affected:
`GraphCTLSchedule`, `GraphMultiHistogram`, `GraphParallelTempering`,
`GraphComputeCorrelationTime`. The returned minimum could sit above energies the run had
visited and recorded, with an empty `minEstates` to go with it.

Bug fixes only. No public function was removed, no argument order changed, no new setting.
Reported energies are now the hamiltonian's own values rather than incremental accumulations,
so charted energies, specific heats and free energies can differ from 1.7.0 in their last bits.

### Fixed
- **The ground-state search reported minima that were not minima.** `GraphSweepReplica` (both
  overloads) tested a state's energy against its running minimum only inside
  `If[accept==1,...]`, so a state could enter the minimum only by being *stepped into* by an
  accepted move. The state a sweep *starts from* was therefore never weighed — and it is the
  state the chain occupies through every rejected move. A sweep that began on a low
  configuration and accepted nothing returned `minEToBeat` itself, an energy no visited state
  ever had, with an empty `minEstates`. Seeded on K4, the ground state of `h[-1.,0.]`, at
  beta 10 where nothing is accepted, the sweep sat at `-12.` for all 120 steps and reported
  `-50.` when asked to beat `-50.`. The sweep now weighs its seed, which is what its own
  documentation always promised ("the minimum energy visited throughout the sweep").
- **Consequently, per-replica minima under parallel tempering were wrong in the reported way:
  the last states visited could have lower energy than the minimum returned.** A replica swap
  moves a configuration into a replica without any accepted move, and under that replica's
  external field its energy is a different number no sweep had ever weighed. Fixing the sweep
  covers every such configuration, since each one seeds the receiving replica's next sweep.
  The configuration held after the *final* swap is not swept again, so the six drivers that
  track minima per replica — `GraphMultiHistogram` and `GraphParallelTempering`, external-field
  overloads, with and without `delH` — now collect it once the measurement loop ends,
  recomputing its energy under that replica's own field. This is what `GraphCTLSchedule`
  returns, via `GraphMultiHistogram`. The beta-swap drivers were not affected: they keep one
  global minimum under a single hamiltonian, so a swap moves no energy they had not seen.
- `GraphComputeCorrelationTime` (both overloads) inherited the same hole directly: it starts
  its running minimum *at* `minEToBeat` with no states and relies entirely on the sweep to
  improve it, so a frozen chain returned the threshold. `GraphEquilibriate` was never
  affected — it initialises from `hamiltonian[seedGraph]` rather than from the threshold.

- **Reported energies are now the hamiltonian's own values rather than accumulations.** The
  sweep builds its energy incrementally (`energy += delE`) and drifts a few ulps from
  `hamiltonian[state]`. On machine reals `==` and `<` share one *relative* tolerance and never
  both hold, so that drift is normally invisible — but a relative tolerance has nothing to
  scale at zero. For a hamiltonian whose ground state sits at exactly `0.`, `GraphSweepReplica`
  reported minima such as `-4.163336342344337*^-16`: a minimum below the least energy the model
  can produce, which then compares unequal to `0.` and drops the states tying it. Both sweep
  overloads now recompute the reported minimum, and the returned state's energy, from the
  hamiltonian once the sweep ends — two calls per sweep, off the hot path, where canonicalising
  inside the step loop measured +24% to +51% (the hamiltonian is a matrix product, 25–50× a
  `delH` step). `SGradDescent` and `SimulatedAnnealing` (both overloads each) canonicalise
  in-loop instead, which measured within noise there because each of their steps already scans
  the whole edge list; for `SimulatedAnnealing` this covers `excitedEnergy` and `excitedStates`
  as well. Note `data`'s energy is now exact, so charted energies, specific heats and free
  energies can differ from earlier versions in their last bits.

### Added
- Two deterministic regression tests, `GraphSweepReplica-counts-its-seed` and
  `GraphComputeCorrelationTime-counts-frozen-seed`, built on a provably frozen chain rather
  than on a pinned seed. Both `GraphMultiHistogram-invariants` and
  `GraphParallelTempering-invariants` now also assert that the reported minimum is at or below
  every energy the run recorded in its own chart, and that `minEstates` is non-empty —
  `groundStatesConsistentQ` passes vacuously when it is empty, which is exactly what the bug
  produced.
- `GraphSweepReplica-canonical-energies`, which asserts with `Order[a,b]==0` that the reported
  minimum is bit-for-bit the hamiltonian of the states reported with it. `==`, `===` and the
  suite's `floatEq` are all tolerant on machine reals and pass on precisely the values that
  test exists to reject.

## [1.7.0] - 2026-08-13

Equilibriation and correlation time, audited and rebuilt. The convergence test that every
Monte Carlo driver depends on was comparing a squared difference of window *means* against the
variance of *individual* samples -- quantities that differ by the effective sample size -- so it
was in effect testing whether the autocorrelation time was short rather than whether the
replicas agreed, and a systematic 1-sigma offset between half of them passed 99% of the time.
It now compares like with like, watches two observables instead of one, says so out loud when
neither varied, and runs once per window rather than once per sweep. The correlation-time stage
no longer takes its run length from the burn-in. Alongside that, the two user settings finally
reach parallel subkernels: `DistributeDefinitions` had been silently shipping nothing.

The complex-space MCMC family also gained real test coverage over this cycle -- all four
definitions of `RandomPureSimplicialComplexMCMC` and both stages beneath it, 172 tests in all.
No public function was removed and no argument order changed. One public setting is new,
`$ECGravMaxCorrelationSweeps`.

**Minimum Wolfram Language version is now 15.0**, up from 14.2. The paclet will not install on
14.x.

Note for anyone relying on reproducibility: equilibriation stops at a different sweep than it
did at 1.6.0, so `eqlT`, and everything downstream of it, differs for a given `SeedRandom`. Any
`eqlT` recorded from an earlier version is an underestimate. `GraphEquilibriate::noconv` gained
an argument slot, and `::nosignal` / `::shortrun` are new messages that can appear on runs that
used to be silent.

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

[Unreleased]: https://github.com/kassabetre/ECGrav/compare/v1.12.0...HEAD
[1.12.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.12.0
[1.11.1]: https://github.com/kassabetre/ECGrav/releases/tag/v1.11.1
[1.11.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.11.0
[1.10.1]: https://github.com/kassabetre/ECGrav/releases/tag/v1.10.1
[1.10.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.10.0
[1.9.1]: https://github.com/kassabetre/ECGrav/releases/tag/v1.9.1
[1.9.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.9.0
[1.8.1]: https://github.com/kassabetre/ECGrav/releases/tag/v1.8.1
[1.8.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.8.0
[1.7.1]: https://github.com/kassabetre/ECGrav/releases/tag/v1.7.1
[1.7.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.7.0
[1.6.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.6.0
[1.5.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.5.0
[1.4.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.4.0
[1.3.1]: https://github.com/kassabetre/ECGrav/releases/tag/v1.3.1
[1.3.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.3.0
[1.2.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.2.0
[1.1.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.1.0
[1.0.0]: https://github.com/kassabetre/ECGrav/releases/tag/v1.0.0
