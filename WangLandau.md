# Wang–Landau sampling for the manifold-emergence transition

> **What this is.** A specification of `GraphWangLandau`, flat-histogram sampling of the density of
> states $`g(E)`$ for `H2dCombManifold` with source terms. It exists to answer one question that
> parallel tempering answers badly: **is the manifold-emergence transition first or second order?**
>
> **Status: specified, implementation in progress.** `Kernel/Subpackages/MCSims.wl`, section
> *Wang–Landau*; usage and messages in `Kernel/ECGrav.wl`, chapter *Flat-Histogram MC*; tests in
> `Tests/MCSims.wlt`. Implements `ExpansionRoadmap.md` Track C (§C1–C5). Companion to
> `HomogeneousHamiltonian.md`, whose source operators supply the energy.

---

## 1. Why flat-histogram, and why now

Canonical sampling at a first-order transition is exactly where importance sampling struggles: the
two phases are separated by a free-energy barrier, and a canonical chain crosses it exponentially
rarely. Parallel tempering mitigates but does not remove this. Flat-histogram sampling attacks the
problem directly by building $`g(E)`$, from which every canonical quantity follows by reweighting:

```math
Z(\beta) \;=\; \sum_E g(E)\,e^{-\beta E}, \qquad
P_\beta(E) \;\propto\; g(E)\,e^{-\beta E}.
```

The order of the transition is then read off $`P_\beta(E)`$ at the transition temperature (§6).

**The ensemble is LABELED**, per `ExpansionRoadmap.md` §0.7. This is not a detail. It makes $`g(E)`$
a plain density of states rather than a $`\Gamma`$-weighted one, so every standard formula applies
unmodified, and it keeps `GroupOrder[GraphAutomorphismGroup[…]]` — whose own source comment flags a
memory leak — out of the inner loop entirely.

---

## 2. The energy spectrum decides the design

§C1 asks whether the spectrum is discrete enough to use exact levels, or whether bin width becomes a
parameter needing justification. **It is discrete, and sparse.** Measured over all $`2^{10} = 1024`$
labeled graphs at $`n = 5`$, with couplings $`(e_0, g_\chi, g_\partial) = (3.0, -0.5, 0.5)`$:

| $`n_{t0}`$ | states | distinct levels | min gap | span |
| --- | --- | --- | --- | --- |
| 0.0 | 1024 | 26 | 1/2 | 88.0 |
| 1.0 | 1024 | 25 | 1/2 | 94.5 |
| 0.7 | 1024 | 27 | 1/10 | 92.4 |
| −1.6 | 1024 | 27 | 3/10 | 78.4 |

Two consequences, and they point in opposite directions:

**Bin width is not a free parameter.** The energy is an integer combination of the source operators
weighted by the couplings, so the gaps are simple rationals fixed by the couplings' denominators
(at $`n_{t0} = 0.7`$, $`\mathrm{lcm}(1,2,2,10) = 10`$ gives gaps in multiples of 1/10). Levels are
exact and there is nothing to justify to a referee. **Keep the couplings rationally related** — an
irrational ratio degenerates the lattice and forces genuine binning.

**But the lattice is mostly empty, and that kills the obvious implementation.** At
$`n_{t0} = 0.7`$ the 27 realised levels sit on a 925-point lattice — **2.92 % occupancy**, with a
median gap of 19/10 against a minimum of 1/10. A dense array indexed by
$`\lfloor (E - E_{\min})/\delta \rfloor`$ would leave 898 slots permanently at $`H = 0`$, so the
flatness test $`\min H \ge c\,\overline{H}`$ could never be satisfied and the run would never
terminate. **`g` and `H` are therefore keyed by energy over realised levels only** (§4.1).

---

## 3. The algorithm

Standard Wang–Landau, with the 1/t modification-factor schedule.

**Acceptance.** A proposed move $`E \to E'`$ is accepted with probability
$`\min\bigl(1,\ g(E)/g(E')\bigr)`$, evaluated in log space:

```math
\Delta S \;=\; \ln g(E) - \ln g(E'), \qquad
\text{accept if } \Delta S \ge 0 \text{ or } \ln u < \Delta S,\ u \sim U[0,1).
```

The log form matches the package's log-space acceptance convention from 1.9.1 and cannot underflow
where $`\exp(\Delta S)`$ would.

**The update is applied at the current state whether or not the move was accepted.** This is the
classic Wang–Landau implementation error: updating only on acceptance biases $`g(E)`$ while leaving
the run looking healthy.

```math
\ln g(E_{\text{cur}}) \mathrel{+}= \ln f, \qquad H(E_{\text{cur}}) \mathrel{+}= 1.
```

**Schedule.** Two, selected by `"Schedule"`:

- `"Classic"` — $`f \to \sqrt{f}`$ on each flat histogram, resetting $`H`$. Simple, textbook, and the
  baseline §7 validates against.
- `"OneOverT"` (default) — run classic until $`\ln f \le N_{\text{lev}}/j`$ for $`j`$ steps, then set
  $`\ln f = N_{\text{lev}}/j`$ every step with no further flatness tests or histogram resets. This
  avoids the **error saturation** of the naive schedule, where the error stops decreasing below a
  floor set by the flatness criterion however long the run continues. §6's double-peak analysis is
  precisely what a saturated error obscures, so 1/t is the production path.

**`"FinalLogF"` does not mean the same thing to both schedules.** For `"Classic"` it is a genuine
convergence criterion. For `"OneOverT"` it is arithmetically a *step budget*: since
$`\ln f = N_{\text{lev}}/j`$, reaching a target $`\ln f`$ requires
$`j = N_{\text{lev}}/\ln f_{\text{final}}`$ steps — at 27 levels, $`2.7\times10^{7}`$ steps for
$`10^{-6}`$ and $`2.7\times10^{9}`$ for the default $`10^{-8}`$. Budget 1/t runs by `maxSteps`.

**Window.** Moves proposing $`E' \notin [E_{\min}, E_{\max}]`$ are rejected, with the update still
applied at the unchanged current state. The window is an argument in this version; automatic
bracketing from a short annealing run is a later overload.

---

## 4. Design decisions

### 4.1 `g` and `H` are keyed by energy, not indexed on a lattice

Forced by the 2.92 % occupancy of §2. The cost argument for arithmetic indexing does not survive
contact with the numbers: an incremental $`\Delta E`$ costs ~450 µs at $`n = 21`$ while an
Association lookup costs ~0.25 µs, so the indexing scheme is 0.06 % of a step. Correctness on a
sparse spectrum is worth vastly more than that.

### 4.2 Energy is accumulated incrementally, not recomputed

Textbook statements of the algorithm write $`E' = \mathrm{Energy}(X')`$ — a full evaluation per
proposal. Measured at $`n = 21`$, against the incremental `delHSourceTerms` + `delNumTriangles`:

| edge density | full $`H`$ | incremental $`\Delta E`$ | ratio |
| --- | --- | --- | --- |
| 0.2 | 2738 µs | 464 µs | 6× |
| 0.4 | 3473 µs | 609 µs | 6× |
| 0.6 | 4588 µs | 253 µs | 18× |
| 0.8 | 9002 µs | 450 µs | **20×** |

At $`10^7`$ steps that is ~25 hours against ~1.25. Full recomputation is not affordable.

### 4.3 New levels are re-anchored from the Hamiltonian

The price of §4.2 is float drift: $`10^7`$ accumulations of `+= dE` can move an energy off its true
level, and a drifted value would register as a **spurious new level** — silently inflating
$`N_{\text{lev}}`$ and corrupting both the flatness test and the 1/t schedule.

The fix exploits the fact that genuine new levels are rare (tens per run). **When a proposal lands on
an energy not already a key, recompute it from the Hamiltonian rather than trusting the running
sum.** If the recomputed value matches an existing key, the discrepancy was drift and is counted in
`levelMerges`; otherwise the level is genuinely new. Tens of full evaluations at 9 ms against
$`10^7`$ steps at 0.45 ms is negligible.

This is cheap *because* the spectrum is sparse. If a larger $`n`$ turns out to have thousands of
levels the assumption weakens and a cheaper discriminator is needed — `levelMerges` and the level
count are reported so this is visible rather than assumed.

### 4.4 A periodic assertion guards `delH`

§4.3 snaps away small discrepancies, which would equally mask a genuinely wrong `delH`. Every
`"CheckInterval"` steps the accumulated energy is asserted equal to `H[A]` recomputed, and the worst
discrepancy seen is reported as `energyDrift`. `GraphSweepReplica` carries the same check,
commented out; here it is live.

### 4.5 Normalisation is by the total state count, not the ground state

$`g`$ is determined only up to a constant. The usual advice — shift so that $`g(E_{\min}) = 1`$ —
is **wrong for this model**, because the labeled ground level is degenerate. Measured at $`n = 5`$,
$`n_{t0} = 0.7`$: $`g(E_{\min}) = 60`$, and only 2 of 27 levels have $`g = 1`$. Ground-state
degeneracy for a labeled graph system is $`N!/|\Gamma|`$ and grows fast.

For a labeled ensemble the total is known exactly, so the constant is fixed exactly:

```math
\sum_E g(E) \;=\; 2^{\binom{v}{2}}.
```

This doubles as a free correctness check — a run whose normalisation is badly off is reporting a
$`g(E)`$ that cannot be right.

### 4.6 The move set matches `GraphSweepReplica` exactly

Single edge toggle, with the same non-uniform-but-symmetric row/column draw
(`row = RandomInteger[{1,v-1}]`, `col = RandomInteger[{row+1,v}]`). Two reasons: the WL and
Metropolis chains then explore by the same moves and are directly comparable, and the later
repointing onto the state backend of `QuantumMCSimsPort.md` is a substitution rather than a rewrite.

State access is isolated into four local functions — `size`, `propose`, `applyMove`, `deltaE` —
matching that backend's primitives. **This is also the point:** those primitives were derived
entirely from Metropolis-flavoured drivers, and Wang–Landau is the first customer with a different
acceptance rule. It needs absolute $`E`$ and an energy grid, which no existing driver does, and
needs nothing from the automorphism ratio or the replica record. If the interface is wrong, it is far
cheaper to discover that here than after repointing 33 overloads.

---

## 5. Return contract

A four-element list, matching `GraphMultiHistogram` and `GraphParallelTempering`:

| slot | contents |
| --- | --- |
| 1 | `<\|"minEnergy", "minEstates"\|>` — same shape as the other drivers |
| 2 | **the density of states** — `<\|"logg", "levels", "window", "normalisation", "totalStates"\|>` |
| 3 | `<\|1 -> <\|"graph", "energy"\|>\|>` — final chain state, for restart |
| 4 | diagnostics — `<\|"hist", "logf", "steps", "phase", "schedule", "converged", "reductions", "energyDrift", "levelMerges"\|>` |

**Slot 2 is not a chart.** In the tempering drivers that slot holds a time series whose third column
is the energy, and `TemperingPlots.wls` reads it as `res[[2]][[All,All,3]]`. Wang–Landau's payload is
a function of energy, not a time series, so that access is meaningless here. Slot 2 is an Association
with named keys specifically so the access **errors** rather than silently returning something
plausible.

`minEstates` is deduplicated by **exact equality, not up to isomorphism**. The labeled ensemble is
what $`g(E)`$ counts, so isomorphism-dedup would report a different object from the one sampled.

---

## 6. Order of the transition

From $`g(E)`$, form $`P_\beta(E) \propto g(E) e^{-\beta E}`$ and locate $`\beta_c`$ by equal-weight
double peak.

- **First order** — the peak separation $`\Delta E`$ scales with system size (latent heat), and the
  depth of the histogram minimum between the peaks grows with the interface free energy.
- **Second order** — a single peak broadens, and the Binder energy cumulant minimum tends to 2/3.

**Report against several $`N`$.** One size cannot distinguish the two.

---

## 7. Validation

- **Exact, not statistical.** At $`n \le 5`$ every labeled graph is enumerable, so $`g(E)`$ is
  computed directly by tallying and compared against the sampled result. This is the standard the
  MCMC generators were held to.
- **Free check.** $`\sum_E g(E) = 2^{\binom{v}{2}}`$ exactly.
- **Independent route.** Thermodynamics reconstructed from $`g(E)`$ against the existing
  multi-histogram / MBAR path on a case where both are valid — two derivations of the same curve.
- **Schedule comparison — measured, not asserted.** At $`n = 5`$, both schedules given the same
  budget of 500 000 steps (~114 s each):

  | schedule | final $`\ln f`$ | max rel. error | mean rel. error |
  | --- | --- | --- | --- |
  | `"Classic"` | $`1.9\times10^{-9}`$ | **10.3 %** | 2.1 % |
  | `"OneOverT"` | $`5.4\times10^{-5}`$ | **4.8 %** | 1.8 % |

  The classic schedule drove $`\ln f`$ four orders of magnitude lower and is still **twice as
  wrong**. That is error saturation exhibited directly: a small $`\ln f`$ is not evidence of
  accuracy under the naive schedule, and any convergence claim resting on it is unfounded.
  Both runs recovered all 27 levels, the exact normalisation, and the exact ground-state
  degeneracy of 60.

---

## 8. Traps

1. **Update $`g`$ and $`H`$ at the current state whether or not the move was accepted.** Updating
   only on acceptance biases $`g`$ and the run still looks healthy.
2. **Do not index a dense lattice.** 2.92 % occupancy means the flatness criterion can never be met
   and the run never terminates (§2).
3. **Do not normalise by $`g(E_{\min}) = 1`$.** Measured $`g(E_{\min}) = 60`$ at $`n = 5`$ (§4.5).
4. **Every numeric option is gated on `NumericQ`, never `_Real`.** An integer-valued argument does
   not match `_Real`, the call returns *unevaluated*, and the acceptance test's deliberate fourth
   branch reads an undetermined comparison as **reject** — so the chain stops moving silently
   instead of erroring. Same trap as `HomogeneousHamiltonian.md` §4.2.
5. **The inner loop needs a step budget.** A flatness criterion that cannot be satisfied — for any
   of the reasons above — otherwise spins forever.

---

## 9. See also

`HomogeneousHamiltonian.md` for the source operators supplying $`E`$; `MBARWeights.md` for the
reweighting path §7 cross-checks against; `ExpansionRoadmap.md` §C for the track this implements;
`QuantumMCSimsPort.md` §2 for the state backend §4.6 anticipates.
