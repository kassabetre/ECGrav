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
- `"OneOverT"` (default) — Belardinelli–Pereyra. A first phase halves $`\ln f`$ each time *every
  level has been visited* since the last reset — weaker than the flatness test, because this phase
  only has to shape $`\ln g`$ roughly — and **after each halving** tests
  $`\ln f \le N_{\text{lev}}/j`$ at step $`j`$. Once that holds it sets
  $`\ln f = N_{\text{lev}}/j`$ every step with no further tests or histogram resets. This
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

Only an energy *inside* the window ever becomes a level (§4.3), and the window is tested on the
incremental energy *before* the key lookup, so an out-of-window proposal never pays for a full
Hamiltonian evaluation. A window the walk presses against also changes the normalisation (§4.5):
states then lie outside it and the total state count no longer applies.

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
count are reported so this is visible rather than assumed. **Measured at $`n = 10`$:** the lattice is
*not* sparse at low energy — 878 levels fill 97.6 % of the 0.1 lattice on $`[10.1, 100]`$ — yet a
300 000-step pilot still recorded zero drift and zero merges, so the recompute remains affordable.

**Two rules the full-window validation could not exercise,** since there no proposal ever leaves
the window:

- *Only an in-window energy registers.* Registering one outside the window plants a level the walk
  can never visit. Its histogram stays at zero, so the Classic flatness test can never pass; it
  inflates $`N_{\text{lev}}`$ in the 1/t rate; and it is returned in the density with
  $`\ln g = 0`$. Measured at $`n = 4`$ on the window $`[2.9, 16.2]`$ before the fix: 200 000 Classic
  steps left $`\ln f`$ at 1.0 with no reduction, 3 of the 8 levels returned lay outside the window,
  and $`g`$ was wrong by 99.9 %. After it, the run converges in about $`10^4`$ steps with the error of
  the full-window run.
- *A new level starts at its parent's $`\ln g`$, not at zero.* Under 1/t the established levels sit
  near $`\ln j`$ above zero, so a level found late at zero would trap the walk until it caught up —
  about $`\ln j / \ln f`$ visits, a share of the run that grows with $`j`$. The parent is one edge
  toggle away and the nearest estimate available; it also makes the first move into the new level
  certain, so every registered level is visited.

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

**Only while the window holds every state.** Once the walk has pressed against the window, states
lie outside it and the total describes more than the levels returned. The run then counts its window
rejections and normalises over the window instead — $`\sum_{E \in W} g(E) = 1`$, reported as
`"WindowRelative"` with `totalStates` Missing. Canonical averages and the §6.2 diagnostics are ratios
and do not notice; $`-\beta F`$ carries an unknown additive constant.

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
| 2 | **the density of states** — `<\|"logg", "levels", "window", "normalisation", "totalStates"\|>`; `"normalisation"` is `"TotalStates"` or, if the window rejected any move, `"WindowRelative"` (§4.5) |
| 3 | `<\|1 -> <\|"graph", "energy"\|>\|>` — final chain state, for restart |
| 4 | diagnostics — `<\|"hist", "logf", "steps", "phase", "schedule", "converged", "reductions", "energyDrift", "levelMerges", "windowRejections"\|>` |

**Slot 2 is not a chart.** In the tempering drivers that slot holds a time series whose third column
is the energy, and `TemperingPlots.wls` reads it as `res[[2]][[All,All,3]]`. Wang–Landau's payload is
a function of energy, not a time series, so that access is meaningless here. Slot 2 is an Association
with named keys specifically so the access **errors** rather than silently returning something
plausible.

`minEstates` is deduplicated by **exact equality, not up to isomorphism**. The labeled ensemble is
what $`g(E)`$ counts, so isomorphism-dedup would report a different object from the one sampled.

---

## 6. Order of the transition

### 6.1 One Wang-Landau run per coupling point — decided by measurement

$`g(E)`$ **cannot be reweighted across the source couplings.** Since
$`E(s) = O_0(s) + n_{t0} O_1(s)`$, changing $`n_{t0}`$ moves each state's energy by the
*state-dependent* amount $`O_1(s)`$, so states degenerate at one coupling split at another. This
follows from linearity; no reweighting trick recovers another coupling from a one-dimensional
$`g(E)`$.

That leaves two routes, and the choice was settled by measuring rather than by argument:

1. **One run per coupling point.** Each is independently validated by the invariants of §7 and the
   runs share no state, so a coupling grid parallelises trivially.
2. **A joint density $`g(O_0, O_1)`$.** Because $`H = c \cdot O`$ is linear in the couplings, one
   joint density would give $`Z(c)`$ for *every* coupling point at once — the flat-histogram
   analogue of what MBAR does for tempering.

Route 2 is attractive only if the joint support stays close to the one-dimensional level count. **At
$`n \le 6`$ it does, and that is a small-$`n`$ accident.** Measuring the conditional multiplicity —
for each observed $`O_0`$, how many distinct $`O_1`$ — over bins holding at least five samples:

| $`n`$ | mean multiplicity | max | samples/bin (median) |
| --- | --- | --- | --- |
| 5 (exact) | 1.04 | 2 | 265 |
| 6 (exact) | 1.20 | 2 | 72 |
| 7 | 1.92 | 6 | 28 |
| 9 | 3.97 | 12 | 10 |
| 12 | 5.01 | 15 | 9 |
| 15 | **5.66** | 14 | 10 |

Sampled and exact agree at $`n = 5, 6`$ (1.04 / 1.21 against 1.04 / 1.20), which is what makes the
larger sizes trustworthy. **The multiplicity grows steadily**, and at $`n = 15`$ it is a *lower*
bound — with ~10 samples per bin one cannot observe more than ten distinct $`O_1`$ values, and the
observed maximum of 14–15 sits at that ceiling. So the joint support is already ~6× the
one-dimensional level count at $`n = 15`$ and still climbing: **route 1**.

*Do not estimate this from raw support counts.* The naive ratio $`|(O_0,O_1)|/|O_0|`$ came out 1.78,
2.36, 1.73, 1.49 at $`n = 7, 9, 12, 15`$ — non-monotonic and *falling*, because unsampled support is
missing from both counts. The conditional multiplicity is the robust statistic.

### 6.2 The diagnostics

From $`g(E)`$, form $`P_\beta(E) \propto g(E) e^{-\beta E}`$ and locate $`\beta_c`$ by equal peak
heights.

- **First order** — the peak separation $`\Delta E`$ scales with system size (latent heat), and the
  depth of the histogram minimum between the peaks, $`\ln(P_{\max}/P_{\min})`$, grows with the
  interface free energy.
- **Second order** — a single peak broadens, and the minimum of the Binder energy cumulant
  $`V = 1 - \langle E^4\rangle / (3\langle E^2\rangle^2)`$ tends to 2/3.

**Report against several $`N`$.** One size cannot distinguish the two.

#### Telling a double peak from lattice granularity

The attainable energies are irregularly spaced on a sparse lattice (§2), so $`g`$ jumps between
adjacent levels and $`P_\beta(E)`$ is jagged. A bare local-maximum test is therefore useless: on the
**exact** $`n = 6`$ density it found **31 "maxima" among 91 levels** and called every $`\beta`$
bimodal. `TransitionOrderDiagnostics` filters candidates by topographic **prominence** — a peak is
kept only if it stands `"MinProminence"` above the higher of the two saddles flanking it, in log
units.

That threshold cannot be chosen once and for all. **Scan it and look for a plateau:**

| | peak count as `"MinProminence"` goes 0.5 → 4.0 |
| --- | --- |
| synthetic double peak | **2, 2, 2, 2** — a plateau |
| exact $`n = 6`$, $`n_{t0} = 0.7`$ | 28, 22, 15, 4 — no plateau |

A genuine bimodality is robust to the threshold; lattice granularity slides away as the threshold
rises. **The $`n = 6`$ result is a negative control, not a measurement**: six vertices cannot host a
transition, so there is nothing there to find, and the diagnostic correctly declines to find it once
granularity is filtered out.

The detector itself is validated against a **synthetic** density built so the answer is known —
$`\ln g = \ln(\text{two Gaussians}) + \beta_0 E`$, making $`P_{\beta_0}(E)`$ exactly that mixture.
It recovers the peak positions, the separation, the log depth and $`\beta_0`$ to float precision,
and returns `Missing["NoBimodalBeta"]` on a single Gaussian.

### 6.3 Where to point the telescope

Two parameters fell out of the $`n = 6`$ scan. They are finite-size numbers, so they say where to
look, not what is there.

**The $`\beta`$ grid must reach at least ~10.** At $`n_{t0} = 0.7`$ the global $`C_v`$ peak sits at
$`\beta = 7.0`$. A first pass capped at $`\beta = 3`$ reported the peak at 0.26 with a twentieth of
the height — a pure truncation artifact. **Report every local maximum, never a bare argmax.**

**The $`n_{t0}`$ grid should be non-uniform.** The $`C_v`$ peak moves $`0.64 \to 1.10 \to 1.70 \to
7.00`$ across $`n_{t0} = -1.6 \to 0.7`$, most of it between 0 and 0.7, and by $`n_{t0} = 1.2`$ the
high-$`\beta`$ structure is gone entirely ($`C_v(\beta = 12) = 0.008`$ against 0.936). Uniform
spacing wastes runs at the flat end and under-resolves the steep one.

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

  **The `"OneOverT"` row measured a defective schedule** — the one of §8 trap 6, whose first phase
  never ran; its final $`\ln f = 27/500\,000`$ exactly says so. At $`n = 4`$, same budget of
  200 000 steps, the defective and corrected schedules give worst errors of 2.5–9.3 % (3 seeds) and
  2.9–4.8 % (6 seeds): a small $`\ln g`$ range hides the defect, which is why the regression test
  counts reductions — 0 before, 8 after — rather than measuring error.

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
6. **Test the 1/t switch only after a reduction.** Tested on every step, $`\ln f \le N_{\text{lev}}/j`$
   holds at $`j = 1`$ — $`1 \le N_{\text{lev}}`$ — so the first phase never runs and the run is 1/t
   from its first step. That converges in principle but not in practice once $`\ln g`$ spans many
   units. If a whole band of $`m`$ levels is off by $`\Delta`$, the walk spends almost all its time
   in that band, so each of its levels gains $`\ln f \cdot N_{\text{lev}}/m`$ per step relative to
   the rest: with $`\ln f = N_{\text{lev}}/j`$ and $`m \approx N_{\text{lev}}/2`$ that is $`2/j`$,
   integrating to only $`2 \ln(j_2/j_1)`$, so clearing $`\Delta`$ takes a factor
   $`e^{\Delta/2}`$ in run length. At $`n = 10`$ the band spans $`\Delta = 5.9`$ across the window
   $`[0, 50]`$ — a factor of about 19 — and more across a wider window, since the span grows with it.
   The first phase clears it in a few traversals at large $`\ln f`$. A 300 000-step $`n = 10`$
   pilot on the defective schedule missed the canonical $`U`$ at $`\beta = 0.5`$ by
   $`z \approx 300`$ and called **every** $`\beta`$ bimodal — but that pilot was also short, 342
   visits per level on 878 levels, so it shows what an unconverged density does rather than
   isolating the schedule. The corrected schedule on $`[0, 50]`$, $`8\times10^6`$ steps: 11 halvings,
   switched to 1/t by step 729 000, histogram min/mean 0.894, and $`U`$ within $`|z| \le 2.3`$ of
   independent canonical chains at all ten $`\beta`$ checked.
7. **Never register an energy outside the window as a level** (§4.3).

---

## 9. See also

`HomogeneousHamiltonian.md` for the source operators supplying $`E`$; `MBARWeights.md` for the
reweighting path §7 cross-checks against; `ExpansionRoadmap.md` §C for the track this implements;
`QuantumMCSimsPort.md` §2 for the state backend §4.6 anticipates.
