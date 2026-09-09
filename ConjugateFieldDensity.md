# The density of the conjugate field

> **What this is.** A specification of `ConstrainedProbConjugateField`, which reweights a
> multi-state Monte Carlo run to a target point and returns the density of the conjugate variable
> there. It has two forms: the original one, for an ensemble of runs at different **external
> fields** at fixed $\beta$, and a beta-only form for an ensemble at different **inverse
> temperatures**, which gives $P_\beta(E)$.
>
> **Status: implemented.** `Kernel/Subpackages/MCSims.wl`; usage and messages in
> `Kernel/ECGrav.wl`; two tests in `Tests/MCSims.wlt`. Companion to `MBARWeights.md`, whose
> helpers the beta form is built on, and to `HomogeneousHamiltonian.md`, where the
> one-component identity §2 rests on is set out.

Notation: the run visited $K$ states, drawing $n_j$ samples at state $j$; $S=\sum_j n_j$ is the
total. $\Phi_j$ is `minusBetaF[j]`, that is $-\beta F$ at state $j$, so $\Phi_j = \log Z_j$ up to
the gauge constant the solver fixes. $O_s$ is the conjugate observable of sample $s$ and $E_s$ its
energy.

---

## 1. What it computes

The reweighted density of the conjugate variable at a target the run never visited. For the
external-field form the target is a field vector $x$ at the run's fixed $\beta$; for the beta form
it is an inverse temperature $\beta_T$, and the conjugate variable is the energy.

Both return a `SmoothKernelDistribution` built on MBAR-weighted samples, together with a plot
range. Nothing is binned: every sample from every state contributes to every target, carrying the
weight that says how likely it would have been there.

## 2. Beta tempering is the one-component homogeneous form

This is the whole content of the beta overload, and it is why there is no new estimator.

Write the hamiltonian in the homogeneous form $H = c\cdot O$. A beta-only run is the case

$$c=\{\beta\},\qquad O=\{E\},\qquad \beta_{\text{fixed}}=1,$$

so $c \cdot O = \beta E$ — exactly the Boltzmann exponent, and exactly the quantity the
external-field form reweights on. The rung betas take the place of the external fields, the energy
takes the place of the conjugate field, and the estimator is unchanged.

| homogeneous form | beta-only run |
| --- | --- |
| `betaFixed` | `1.0` |
| field key $h_j$ | rung inverse temperature $\beta_j$ |
| conjugate observable $O_s$ | energy $E_s$ |
| target field $x$ | target inverse temperature $\beta_T$ |
| $c\cdot O$ | $\beta E$ |

`BetaSurfaceData` in `TemperingPlots.wls` rests on the same identity.

## 3. The estimator

The MBAR weight of sample $s$ at target $\beta_T$ is

$$\log u_s(\beta_T) \;=\; -\beta_T E_s \;-\; D_s,\qquad
D_s \;=\; \log\sum_{j=1}^{K} \exp\!\big(\log n_j - \Phi_j - \beta_j E_s\big),$$

$D_s$ being the log of the MBAR mixture density at $E_s$. The density is then the weighted kernel
estimate

$$P_{\beta_T}(E) \;=\; \frac{\sum_s u_s\,K_h(E-E_s)}{\sum_s u_s},$$

and the effective sample size is Kish's

$$\mathrm{ESS}(\beta_T)\;=\;\frac{\big(\sum_s u_s\big)^2}{\sum_s u_s^2}.$$

Every quantity here is a ratio, so the weights need no normalisation and the free energy of the
target never has to be computed.

**$D_s$ carries no target.** That is the point of the list form. The sum over $K$ states — the only
part quadratic in $K$ — depends on the sample and not on $\beta_T$, so it is built once by
`MBARWeightBasis` and reused at every target, where the per-target step is one matrix–vector
product. Sweeping $\beta$ is what this overload exists for, and this is the reuse `MBARWeights.md`
was written for. The four-argument form still carries the pre-1.11.0 nested `Sum`, which repeats
that work at every call.

The helpers need no adaptation: `MBARWeightBasis` already coerces bare-real keys and scalar
measurements to matrices, so a beta chart can be handed to it as it comes off the run.

## 4. The plot range, and why it is computed here

The beta form does **not** reuse the four-argument form's range, and the reason is not stylistic.

That form widens by $0.8\,\min$ and $1.2\,\max$, then runs two sign-fixing passes in which the
second reads the minimum the first has already overwritten. Those factors widen an interval only
when the data is positive. **Energies are negative**, so they shrink it, and the sign fixes then
invert it. On a real run with energies in $[-60,-2]$:

| step | min | max |
| --- | --- | --- |
| after $0.8\min$ / $1.2\max$ | $-48.0$ | $-2.4$ |
| after the two sign fixes | $0.48$ | $-0.096$ |

leaving $\min > \max$ and neither bound near the data. The beta form pads by a fraction of the
**span** instead,

$$\text{lo}=\min E-0.2\,(\max E-\min E),\qquad \text{hi}=\max E+0.2\,(\max E-\min E),$$

which is sign-agnostic. The range does not depend on the target, so the list form returns one
range shared by every $\beta_T$ — which is what makes a stack of densities, or a `DensityPlot` of
$P(E;\beta)$, coherent.

## 5. The bandwidth, and the energy lattice

A kernel density estimate places a bump of width $h$ — the **bandwidth** — on every sample and
adds them up, each scaled by its weight. Whether the result reads as a curve or as a comb is
decided by $h$ against $\Delta$, the spacing between *attainable* energies:

$$h \gtrsim \Delta \;\Rightarrow\; \text{adjacent bumps overlap and merge};\qquad
h \ll \Delta \;\Rightarrow\; \text{one separated bump per occupied level}.$$

The energy of a graph Hamiltonian is lattice-valued, so $\Delta$ is finite and this is a live
concern rather than a technicality. `SmoothKernelDistribution` chooses $h$ by Silverman's rule
applied to the weighted sample, roughly $h \approx 1.06\,\sigma_w\,\mathrm{ESS}^{-1/5}$. That rule
was derived for continuous data and has no notion of a lattice — and $\sigma_w$ **collapses as the
target gets colder**, because the reweighted energy distribution narrows onto fewer levels.
Measured on a run with $\Delta = 1$:

| $\beta_T$ | ESS | levels carrying $>0.1\%$ | $\sigma_w/\Delta$ | $h$ | $h/\Delta$ |
| --- | --- | --- | --- | --- | --- |
| 0.05 | 151 | 71 | 19.8 | 3.57 | 3.6 |
| 0.15 | 271 | 72 | 23.1 | 4.24 | 4.2 |
| 0.35 | 1247 | 7 | 7.9 | 1.46 | 1.5 |
| 0.6 | 2254 | 2 | 1.83 | 0.50 | **0.5** |
| 1.0 | 2168 | 1 | 0.17 | 0.50 | **0.5** |

Counting local maxima over a fixed window shows the same thing directly: $h=0.5$ gives 12 bumps,
$h=1.0$ gives 5, $h=1.5$ gives 1, $h=2.0$ gives none.

**Note this is not fixed by rescaling the observable.** Silverman's rule is scale-invariant: divide
the energy by any constant and $h$, $\sigma_w$ and $\Delta$ all shrink together, leaving $h/\Delta$
exactly where it was. What matters is how many lattice levels the *reweighted* distribution spans,
which is a property of the target, not of the units.

Hence the two options and the warning:

- `"Bandwidth"` overrides the automatic choice. A value between $\Delta$ and $2\Delta$ is usually
  right: wide enough to bridge adjacent levels, narrow enough to keep genuine structure. Beyond a
  few $\Delta$ real features start being washed out.
- `ConstrainedProbConjugateField::lowbandwidth` fires when the bandwidth actually used falls below
  $\Delta$, reporting both numbers. It is raised **once per call**, naming the narrowest bandwidth
  in the sweep, because a sweep over sixty betas would otherwise emit sixty copies and then be cut
  off by `General::stop`.
- `"Form" -> "PMF"` sidesteps the question: for a lattice-valued energy the exact per-level
  probabilities are the true object, and the density is only a way of drawing them.

$\Delta$ is measured as the smallest positive gap between distinct sample values. For a genuinely
continuous observable the closest pair out of thousands of samples is minute, so $h/\Delta$ is
enormous and the warning cannot misfire. Floating-point noise between two values meant to be equal
shrinks $\Delta$ and *suppresses* the warning — a false negative, which is the safe direction.

## 6. What is refused

Each of these otherwise fails far downstream, as a symbolic expression or an error at plot time
rather than at the call.

| condition | message | why it cannot be allowed through |
| --- | --- | --- |
| `minusBetaF` and the measurements disagree on keys | `::keys` | `Lookup` yields `Missing`, which `Join` and `N` carry into a symbolic observable matrix; that reaches `Dot`, does not evaluate, and returns a blob rather than an error |
| measurements carry more than one component | `::notbeta` | the target $\{\beta_T\}$ is then the wrong width, `Dot` does not evaluate, and every density is a symbolic expression holding the whole sample matrix |
| every measurement list empty | `::nosamples` | there is nothing to estimate, and `MinMax[{}]` fails |
| every energy identical | `::degenerate` | a point mass, not a density: `SmoothKernelDistribution` returns a `DataDistribution` that looks fine and raises `InverseFourier::fftl` the first time `PDF` is asked for a value |
| `"Form"` is neither `"Density"` nor `"PMF"` | `::form` | an unrecognised form would otherwise be silently treated as a density |

**The component check must run before `MBARWeightBasis`, not after.** That helper forms
$O\cdot h^{\mathsf T}$ itself, so a wrong width raises `Dot::dotsh` inside it and hands back a
symbolic basis — at which point the width can no longer be read off it. The test asserts that
`Dot::dotsh` does not escape.

## 7. Reading the ESS

The list form returns $\mathrm{ESS}(\beta_T)$ per target. It is the number of independent samples
the reweighting has left at that $\beta$, and it falls away quickly outside the sampled ladder: a
density at a $\beta_T$ where the ESS is a handful of samples is extrapolation, whatever it looks
like. It is the same quantity `MBARSurfaces` masks its surfaces on.

Two things inflate it. Charts record every sweep, and consecutive sweeps are correlated — at a rung
whose `corrT` is a large fraction of the run length the ESS overstates the independent content by
roughly that factor. And the ESS says nothing about whether the states *overlap in the right
region*: a target between two well-sampled rungs can carry a healthy ESS while the tail of the
density it is asked about is supported by very few samples.

## 8. Getting the inputs from a run

Energy is chart column 3 in both the pre- and post-1.13.0 beta layouts — it is the only column the
package itself reads out of a beta chart — so for a run `res`:

```wl
en  = res[[2]][[All, All, 3]];
mbf = ComputeMinusBetaTimesFreeEnergy[en];
ConstrainedProbConjugateField[Subdivide[0.1, 5.0, 60], mbf, en]
```

This is what the beta drivers do internally. Run unpacking is deliberately not done inside the
function: `MCSims.wl` takes prepared associations throughout, and unpacking a driver return lives
in the `.wls` layer with `BetaSurfaceData` and `MixingData`.

## 9. Verification

| check | result |
| --- | --- |
| pointwise identical to the four-argument form on the same data expressed both ways | $<10^{-12}$; on one fixture the means agreed to `0.` |
| weights against the old nested-`Sum` formula, normalised | $2.6\times10^{-18}$ |
| $\langle E\rangle$ against $-\,\mathrm{d}(-\beta F)/\mathrm{d}\beta$, differenced from `NegativeBetaTimesFreeEnergy` | $4.4\times10^{-7}$, at $\varepsilon=10^{-4}$ |
| reweighting to a sampled rung reproduces its plain sample mean | $-59.04$ against $-58.92$, the gap being kernel smoothing |
| returned range is ordered and contains the data; PDF integrates over it | $1.0000000000000009$ |
| list form against mapping the scalar form | ESS and range identical, PDF difference `0.` |
| ESS inside the ladder exceeds ESS far outside it | asserted, so the ESS cannot be a constant |
| `::lowbandwidth` fires at a cold target and is **silent** at a hot one | both asserted, so a warning hardwired either way fails |
| `"Bandwidth"` reaches the estimator | asserted behaviourally: a wider kernel gives strictly fewer local maxima |
| PMF sums to 1, carries one atom per attained level, and reproduces the MBAR weighted mean | $<10^{-12}$, exact count, $<10^{-9}$ |

Three tests in all. Each behaviour was confirmed by reverting it separately: moving the component check after the
basis build, and restoring the four-argument range heuristic, each independently fails
`ConstrainedProbConjugateField-beta-form-range-and-guards`; and disabling the warning, forcing it
to fire unconditionally, ignoring `"Bandwidth"`, or leaving the PMF unnormalised each
independently fails `ConstrainedProbConjugateField-bandwidth-warning-and-PMF`. Forcing the warning
on also breaks the other two tests, which are silence checks by construction.

## 10. API

```wl
(* external-field form: the density of the conjugate field at a target field, fixed beta *)
ConstrainedProbConjugateField[betaFixed, targetExtField, minusBetaF, measurements]
    ==> {distribution(s), min, max}

(* beta-only form: the density of the ENERGY at a target inverse temperature *)
ConstrainedProbConjugateField[targetBetas_List, minusBetaF, energyMeasurements, opts]
    ==> {distributions, min, max, ess}
ConstrainedProbConjugateField[targetBeta_?NumericQ, minusBetaF, energyMeasurements, opts]
    ==> {distribution, min, max, ess}

    "Bandwidth" -> Automatic (default) | anything SmoothKernelDistribution accepts
    "Form"      -> "Density" (default) | "PMF"
```

Under `"Form" -> "PMF"` the first return value is a sorted `Association` `<|E -> p|>` per target
instead of a distribution — the exact reweighted probability of each attained energy. The other
three values are unchanged, and the scalar form unwraps as before. Plot it as atoms, not as a
curve:

```wl
{pmf, lo, hi, ess} = ConstrainedProbConjugateField[4.0, mbf, en, "Form" -> "PMF"];
ListPlot[Normal[pmf], Filling -> Axis, PlotRange -> All]
```

The external-field form's `min` and `max` are vectors, one entry per component of the conjugate
field. The beta form's are plain numbers, the conjugate variable being scalar by construction, so
`Plot[PDF[d, e], {e, min, max}]` takes them directly; `ess` is a list in the list form and a number
in the scalar one.

The beta form's keys are bare inverse temperatures and its values are flat energy lists — the shape
a beta chart already has, and the shape `ComputeMinusBetaTimesFreeEnergy` and
`NegativeBetaTimesFreeEnergy` already take. Its argument order matches the latter.

## 11. What this does not do

- **It does not check that the target is inside the sampled range.** It returns the ESS and leaves
  the judgement to the caller; there is no mask and no warning.
- **It does not decorrelate.** Samples are taken as given, so a chart written every sweep is used
  every sweep, and §7 applies.
- **It does not choose a bandwidth for you beyond what `SmoothKernelDistribution` does.** §5 is the
  whole of what it offers: a warning when the automatic choice falls below the lattice spacing, an
  override, and an exact discrete alternative. It does not silently substitute a lattice-aware
  rule, because that would change results for every existing caller.
- **The external-field form is untouched**, including its range heuristic. §4 is the reason the
  beta form does not delegate to it; whether that heuristic should itself be fixed is a separate
  question, and fixing it would change the plot ranges of existing external-field callers.
