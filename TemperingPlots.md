# Plotting a tempering run

`TemperingPlots.wls` turns the return of `GraphCTLSchedule` or `GraphParallelTempering` into
continuous surfaces over the physical $(\beta, n_{t0})$ plane, with the per-replica values overlaid
as the reference they must pass through.

It is the batched form of the recipes in `HomogeneousHamiltonian.md` §7.1–7.4. Those give one
surface at one point; this gives every surface over a whole grid in a single pass, because the MBAR
weights do not depend on the observable and do not have to be rebuilt per point.

```wl
Needs["ECGrav`"];
Get["/path/to/ECGrav/TemperingPlots.wls"];
```

## 1. Requires 1.11.0 or later, and it fails silently below that

The pipeline is built on the private MBAR helpers `MBARWeightBasis`, `MBARWeights`,
`MBARWeightedMean`, `MBARMomentsAt` and `MBARESSAt`, which were added in **1.11.0**. On an older
paclet `ECGrav`Private`MBARWeightBasis` has no definition at all, so

```wl
basis = ECGrav`Private`MBARWeightBasis[bt, mbf, conj]
```

returns *unevaluated* — the symbol with its arguments, wearing the head as data. Nothing errors.
Every later definition built on `basis` then inherits the unevaluated head, `ListPlot3D` is handed a
list of symbolic expressions, and the cell produces no output and no message. Check with

```wl
PacletObject["ECGrav"]["Version"]
```

and if it is below 1.11.0, install the current build and restart the kernel:

```wl
PacletInstall["/path/to/ECGrav/build/ECGrav-1.12.1.paclet", ForceVersionInstall -> True]
```

`CTLSurfaceData` checks for this and says so rather than returning a broken bundle.

## 2. The four steps

```wl
data = CTLSurfaceData[res,
         "ObservableNames" -> {"numEdges", "numTriangles", "EulerChi", "numBdryEdges",
                               "AvgKDim", "PureGraphQ", "DGraphQ", "ClosedDGraphQ"}];

data = AddObservable[data, "mag", #numEdges/78. &];      (* derived quantities *)

grid = SurfaceGrid[{0.2, 4.0, 40}, {-1.2, 1.0, 40}];     (* beta range, nt0 range *)

surf = MBARSurfaces[data, grid, "SiteCounts" -> <|"mag" -> 78|>];

plots = SurfacePlots[data, grid, surf, "SiteCounts" -> <|"mag" -> 78|>];
```

`res` is the `{schedule, hist}` pair `GraphCTLSchedule` returns, or `hist` itself.
`"ObservableNames"` names chart column 6 in the order the `obs` list was given to the driver; the
names become the surface names. Without it the columns are `"obs1"`, `"obs2"`, …

For a `GraphParallelTempering` return, the only change is the entry point:

```wl
data = TemperingSurfaceData[res, "ObservableNames" -> {...}];
```

Its chart has the same six columns but carries no free energies, so they are solved with
`ComputeMinusBetaTimesFreeEnergy` — sub-second since 1.12.0.

### 2.1 Beta-only runs

A beta ladder does not look like an external-field run — its chart column 4 is the *observable
list*, not the conjugate observables, and it has no residual column — so reading it with
`TemperingSurfaceData` is silently wrong. Use:

```wl
data = BetaSurfaceData[res, "ObservableNames" -> {...}];
```

Everything downstream is unchanged, because beta tempering **is** the one-component homogeneous
form: the field is $c = \{\beta\}$, the conjugate observable is the energy, $\beta_{\text{fixed}} = 1$,
and so $c\cdot O = \beta E$. `BetaSurfaceData` performs exactly that mapping, and the surfaces are
then functions of $\beta$ alone — pass targets as `{{beta}}`.

Three details it handles that are easy to get wrong by hand.

**Both row layouts.** Since 1.13.0 a beta row is `{sweep, beta, E, {obs…}, tag}`; before that the
observables were spliced in, giving `{sweep, beta, E, obs…}`. Both are read.

**Observables that cannot be averaged are dropped, with a message naming their positions.** An
`obs` entry returning a graph is the usual cause — it makes the chart enormous and there is nothing
to average. Judged on the first row of *every* replica, so a column that is a graph in one and a
number in another cannot slip through.

**The free energies are always re-solved**, never taken from the run, because a run made before
1.12.0 stored under-converged ones (§8). On a $K = 11$, $N = 2000$ run the whole call takes 0.54 s.

Verified on both layouts: on a 1.13.0 run the solved `minusBetaF` matched a separately computed one
to exactly `0.`, and $\langle E\rangle = -\partial\log Z/\partial\beta$ reproduces to $10^{-10}$.

## 3. What comes out

`MBARSurfaces` returns an association of names to value lists, one value per grid point, in the
grid's own order. `SurfaceNames[data]` lists them.

| name | quantity |
| --- | --- |
| `"ESS"` | MBAR effective sample size — how much data stands behind the point |
| `"logZ"` | $\log Z = -\beta F_{\text{phys}}$ |
| `"betaF"` | $\beta F_{\text{phys}}$ |
| `"FPhys"` | $F_{\text{phys}}$ |
| `"EPhys"` | $\langle H_{\text{phys}}\rangle$ |
| `"EHom"` | $\beta\langle H_{\text{phys}}\rangle = \langle c\cdot O\rangle$ |
| `"Cv"` | $\beta^2\,\mathrm{Var}[H_{\text{phys}}]$, i.e. $C/k_B$ |
| `"BinderEnergy"` | $1 - \langle H^4\rangle/(3\langle H^2\rangle^2)$ |
| `"O1"`, `"O2"`, … | the conjugate observables |
| *name* | each registered observable |
| `"Var"`*name* | its variance |
| `"chi"`*name* | $\beta\,N\,\mathrm{Var}$, with $N$ from `"SiteCounts"` (default 1) |
| `"Binder"`*name* | its Binder cumulant |

`BootstrapValues[data]` gives the same names evaluated directly at the $K$ bootstrap field points
from each replica's own samples, as `{beta, nt0, value}` triples. `SurfacePlots` overlays them.

## 4. Why it is fast

Measured on a $K=16$, $N=800$ run — $S = 12800$ samples — over a $12\times12$ grid:

| | time |
| --- | --- |
| `ExtrapolatedExpectationValue` per point, 13 surfaces (the old notebook, serial) | **16.3 s** |
| `MBARSurfaces`, all 52 surfaces | **0.085 s** |

That is 190x on the wall clock for four times as many surfaces — about 770x per surface-point —
and it needs no parallel kernels. Three things account for it.

**The basis is built once.** That is §7.1's advice and worth 10x on its own: the same 144 points for
one observable take 0.83 s through the public entry point and 0.078 s with the basis hoisted.

**The grid is one matrix product.** The MBAR weight of sample $s$ at target $c$ is
$\exp[-(\beta\, c\cdot O_s + \texttt{logBasis}_s)]$ up to a per-target constant that cancels out of
every ratio. For $T$ targets that whole exponent is one $S\times T$ product, not $T$ separate
$S$-vectors.

**Every surface rides the same weights.** MBAR reweighting is observable-agnostic, so all 52
surfaces are matrix products against one weight matrix. Asking for 52 costs what asking for one
costs: 0.078 s either way.

Scaling, all 52 surfaces: $12^2$ points 0.085 s, $30^2$ 0.24 s, $40^2$ 0.45 s, $60^2$ 0.87 s.
Peak memory stays near 240 MB — `"ChunkSize"` (default 512 targets) bounds the $S\times T$
intermediate.

## 5. Traps this handles for you

**The internal energy's observable moves with the plot coordinate.**
$H_{\text{phys}}(n_{t0}) = O_0 + n_{t0}O_1$, so there is no single energy list to reweight
(§7.3). The pipeline takes the outer powers of $O$ once and contracts them with
$v = \{1, n_{t0}\}$ per target, which gives $\langle H\rangle$, $\langle H^2\rangle$ and
$\langle H^4\rangle$ — and so the internal energy, the heat capacity and the energy Binder cumulant
— exactly, at every point.

**Units.** The stored $\beta$ is 1.0 and chart column 3 is $c\cdot O = \beta H_{\text{phys}}$.
`Mean[chartColumn3]` is therefore $\beta\langle H_{\text{phys}}\rangle$, and
`SpecificHeat[chartColumn3, 1, beta]` is $\beta^2$ times the real heat capacity. **Per-replica plots
built that way disagree with the reweighted surface by exactly those factors** — verified on the
$K=16$ fixture, where the ratio of the two heat capacities is $\beta^2$ to the last digit (up to
`Variance`'s $n-1$). Every surface and every bootstrap value here is in physical units.

**Masking.** Outside the sampled region the surfaces return `Indeterminate`, below an ESS floor of
half a replica's worth of samples (`"ESSFloor"`, `None` to disable). On a regular grid
`SurfacePlots` passes the values to `ListPlot3D` as an **array**, which leaves a genuine hole;
handed triples instead, `ListPlot3D` drops the masked points and triangulates straight across the
gap, drawing exactly the invented surface the mask exists to suppress.

**Key order.** Sample vectors are flattened with `Lookup[..., Key /@ Keys[mbf]]`, never `Values`.
The basis stacked its samples in that order and a mismatch is silent.

**Boolean observables.** `PureGraphQ` and friends are `True`/`False`; they are coerced with `Boole`.
A replica whose column is all `False` has $\langle x^2\rangle = 0$ and no Binder cumulant, which
comes back as `Indeterminate` rather than as `Power::infy` at every point.

**The Binder cumulants use raw moments**, $1 - \langle x^4\rangle/(3\langle x^2\rangle^2)$, which is
the form the old notebook used and is the standard one only for an observable whose mean vanishes by
symmetry. `EulerChi` and `numBdryEdges` both cross zero on this run, and where
$\langle x^2\rangle \to 0$ the ratio diverges — the spikes at the cold, low-$n_{t0}$ corner of
`"BinderEulerChi"` and `"BindernumBdryEdges"` are that, not a reweighting failure. Read
`"Binder"`*name* only for quantities where the form means something; `"Var"`*name* is always safe.

**Single-field runs** key their chart on bare reals rather than one-element lists. Both shapes are
accepted, and the physical coordinate is then just $\beta$ — but `SurfacePlots` builds `ListPlot3D`
and wants two axes, so plot a one-dimensional grid yourself from `MBARSurfaces`'s output.

**`"SiteCounts"` belongs on both calls.** `MBARSurfaces` and `SurfacePlots` each need it, or the
overlaid points are drawn with $N=1$ against a surface built with the real site count and sit an
order of magnitude below it — which reads as a reweighting error and is not one.

## 6. Verification

Against the old per-point method on the $K=16$, $N=800$ run, over the $12\times12$ grid, for all
thirteen surfaces the old notebook drew — free energy, internal energy, specific heat, both
conjugate observables and their susceptibilities, `numTriangles`, `EulerChi`, `mag` and their
susceptibilities:

> max relative difference **$0$ to $2.7\times10^{-11}$**, the largest being the variances, where the
> difference of two large numbers costs the digits.

Independently, $\langle E_{\text{phys}}\rangle = -\partial\log Z/\partial\beta$ at fixed $n_{t0}$
reproduces to $10^{-10}$ by central differencing, which checks the coordinate change, the units and
the extraction at once. And `numTriangles` is both `obs[[2]]` and the second conjugate observable,
so `"numTriangles"` and `"O2"` are computed by two different routes and agree to $7\times10^{-15}$.

## 7. Diagnostics

`SwapAcceptance[res]` reads `swapAccept/swapTry` per replica out of the run's histories. Constant
thermodynamic length is supposed to equalise it, so an uneven profile says the schedule did not do
its job. It is necessary and not sufficient: what matters is **round trips**, and nothing counts
those yet — the data is in `hist[[4]][[k, "history"]]`, a buffer of the field values replica index
$k$ has held, updated only on accepted swaps and already trimmed of its `-1.0` padding by the
driver.

## 8. A caveat on runs made before 1.12.0

`hist[[2, All, "minusBetaF"]]` is whatever the run stored. Runs made before 1.12.0 used the old
convergence test, which stopped roughly 12x short; on the fixture here the stored free energies are
**$5.3\times10^{-5}$** away from the converged values, and the surfaces move by up to
$1.6\times10^{-4}$. Re-solving costs 0.6 s and needs no new sampling:

```wl
mbf = ComputeMinusBetaTimesFreeEnergy[hist[[2, All, "data", All, 4]], 1.0];
```

`TemperingSurfaceData` always solves, so it is one way to get the current numbers out of an old run.

## 9. See also

- `HomogeneousHamiltonian.md` §6 (units), §7.1–7.4 (the per-point recipes), §8 (practical limits).
- `MBARWeights.md` — the shared weight vector and the derivation the batching rests on.
- `CHANGELOG.md` 1.10.1 — the four multi-field schedule fixes this recipe depends on.
