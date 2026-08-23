# Tempering in temperature and external field at once

> **What this is.** A recipe, and its justification, for running parallel tempering over
> **both** inverse temperature and an external field using only the shipped API — no new
> machinery. It works by rewriting the hamiltonian in *homogeneous* form, so that $\beta$ becomes
> one of the tempered couplings rather than a fixed scalar. Companion to `Theory.md`. Everything
> here is verified against `Kernel/Subpackages/MCSims.wl` and pinned by a test in
> `Tests/MCSims.wlt`.
>
> **Requires 1.10.1 or later.** Four defects in the multi-field schedule builder were fixed in
> that release, and every one of them is either the identity or unreachable with a single
> external field — which is to say, they are exactly the defects this recipe would have walked
> into. On 1.10.0 or earlier the construction below is sound but the schedule built from it is
> not. See §9.

Notation: $\beta$ inverse temperature; $H$ the physical hamiltonian; $O_i$ the observables
conjugate to the external fields $d_i$; $c$ the field vector the package tempers over; `bt` the
package's fixed inverse-temperature argument. Vectors are written without decoration and $\cdot$
is the dot product.

---

## 1. The problem

The package tempers along one axis at a time, and the two axes are mutually exclusive by
construction. The beta-tempering overloads of `GraphParallelTempering` require

```wl
Length[DeleteDuplicates[inputReplicas[[All,"beta"]]]] > 1 &&
Length[DeleteDuplicates[inputReplicas[[All,"externalField"]]]] == 1
```

(`MCSims.wl:7521`, `7696`) and the field-tempering overloads require the reverse
(`MCSims.wl:8478`, `8729`). So there is no ladder that varies both.

Multidimensional tempering *does* exist — the field vector may have any number of components,
with one conjugate observable each, and the swap topology is an arbitrary graph rather than a
1-D chain. What it cannot do is include $\beta$, because $\beta$ multiplies the whole hamiltonian
while the fields enter linearly.

This matters when a single fixed $\beta$ cannot serve. At high $\beta$ the chain freezes; at low
$\beta$ it never reaches the states of interest. Tempering in the field alone does not help,
because every replica shares the one $\beta$ and so every replica freezes together.

## 2. The construction

### 2.1 The idea

The requirement the multi-field path actually imposes is that the hamiltonian be **linear in the
tempered couplings**:

$$H(x;c) = \sum_i c_i\, O_i(x)$$

Nothing says one of those couplings cannot be $\beta$ itself. Split the physical hamiltonian into
the part $\beta$ multiplies and the field terms,

$$H_{\text{phys}}(x; d) = O_0(x) + \sum_i d_i\, O_i(x),$$

promote $O_0$ to a conjugate observable with its own coupling $c_0$, and set **`bt = 1`**:

$$H_{\text{hom}}(x;c) = c_0\,O_0(x) + \sum_i c_i\, O_i(x).$$

Then, identically,

$$e^{-1\cdot H_{\text{hom}}(x;c)} \;=\; e^{-\beta H_{\text{phys}}(x;d)}
\qquad\text{under}\qquad c_0 = \beta,\quad c_i = \beta\, d_i,$$

with the inverse map $\beta = c_0$, $d_i = c_i/c_0$. A $(\beta, d)$ grid is laid down in
$c$-space as $\{\beta,\ \beta d\}$, and the existing multi-field machinery then tempers in
temperature and field simultaneously.

Geometrically the reachable set in $c$-space is a **fan** through the origin — the rays are
lines of constant $d$ — not a box. That matters in §8.

### 2.2 The case in hand

For the hamiltonian used in the manifold-targeting runs,

$$H_{\text{phys}} = H_{\text{2dCombManifold}}(e_0) + n_{t0}\,T + g_N\,\chi + g_{N,\partial}\,B,$$

with $T$ the triangle count, $\chi$ the Euler characteristic and $B$ the boundary-edge count, the
split is

| | |
| --- | --- |
| $O_0$ | `H2dCombManifold[am,e0] + gN0*EulerChi[am] + gNBdry0*numBdryEdges[am]` |
| $O_1$ | `numTriangles[am]` |

so $c = (\beta,\ \beta\,n_{t0})$ and `bt = 1.0`. Note $g_N$ and $g_{N,\partial}$ stay at their
fixed physical values *inside* $O_0$ — they are not tempered. Promoting them to their own axes is
possible in principle but see §8 on why two components is the practical ceiling.

## 3. Why it is exact

Three separate weights have to transform correctly. All three do, and none of it is a
coincidence — each follows from the same identity $\beta H_{\text{phys}} \equiv 1\cdot
H_{\text{hom}}$.

### 3.1 The single-replica Metropolis weight

`GraphSweepReplica` reduces the entire acceptance decision to one line
(`MCSims.wl:2181`, `2347`):

```wl
delE = delH[data[[Key["graph"]]], delHparams, row, col];
logRatio = Log[N[selectionProb]] - delE*beta;
```

**`beta` appears nowhere else.** It enters only as the product `delE*beta`, so a hamiltonian
scaled by $\beta$ and evaluated at `bt = 1` gives a bit-identical decision. The code cannot tell
the difference.

This is pinned by `GraphSweepReplica-homogeneous-reparameterisation-is-identical` in
`Tests/MCSims.wlt:322`, so a future change that gave `beta` a second use would fail rather than
silently invalidate this recipe.

`selectionProb` is the automorphism ratio used for unlabeled runs and is hamiltonian-independent,
so it is untouched. The numeric range is unchanged too — $\beta\,\Delta H$ and
$1\cdot\Delta H_{\text{hom}}$ are the same magnitude — so the log-space acceptance introduced in
1.9.1 continues to do its job.

### 3.2 The replica-swap weight

The multi-field swap criterion is the vector form (`MCSims.wl:8016`, `8314`, `8569`, `8819`):

```wl
betaDelHDelConj = bt*(c_i - c_j) . (O(x_i) - O(x_j))
```

Substituting $c_0 = \beta$, $c_k = \beta d_k$ and `bt = 1`:

$$(\beta_i - \beta_j)\,\Delta O_0 \;+\; \sum_k (\beta_i d_{k,i} - \beta_j d_{k,j})\,\Delta O_k,$$

which is exactly the correct acceptance for swapping two replicas at different temperatures
*and* different fields. It reduces to $(\beta_i-\beta_j)\Delta H$ when the fields agree, and to
$\beta(d_i-d_j)\Delta O$ when the temperatures do — i.e. to each of the two shipped cases.

The swap topology is drawn from an arbitrary graph by
`ChooseRandomIndependentEdgeSet[replicaSwapEdgesLabels[[1]]]`, so a 2-D ladder needs no special
handling.

### 3.3 MBAR and the thermodynamic metric

MBAR reweighting assumes *precisely* this linear structure. It reweights through

```wl
Exp[betaFixed*(targetExtField - j) . conjugateExtFieldMeasurements]
```

on the conjugate-observable column of the chart, so with the whole hamiltonian written as
$c\cdot O$ the reweighting is exact rather than approximate.

Better still, with `bt = 1` the pair $(c, O)$ is the natural-parameter / sufficient-statistic
pair of an exponential family, $P \propto e^{-c\cdot O}$. The thermodynamic metric the schedule
builder constructs — the covariance $\Sigma$ of the conjugate observables, `bt^2*Cov(O)` — **is**
the Fisher information in that parameterisation. The homogeneous form is not a trick played on
the API; it is the parameterisation the machinery is written for.

## 4. The code

These are the definitions the §5 equivalence checks were run against. The helpers they call —
`numBdryEdges`, `delnumBdryEdges`, `numTriangles`, `delnumTriangles`, `delEulerChi` — are the ones
already used in the run notebooks and are unchanged by any of this.

### 4.1 The conjugate observables and their deltas

```wl
(* O0 -- everything beta multiplies, at the fixed physical couplings *)
obs0[am_List] := ECGrav`H2dCombManifold[am, e0]
    + gN0*ECGrav`EulerChi[am] + gNBdry0*numBdryEdges[am];
delObs0[am_List, i_Integer, j_Integer] := ECGrav`delH2dCombManifold[am, e0, i, j]
    + gN0*delEulerChi[am, i, j] + gNBdry0*delnumBdryEdges[am, i, j];

(* O1 -- the field-conjugate observable *)
obs1[am_List] := numTriangles[am];
delObs1[am_List, i_Integer, j_Integer] := delnumTriangles[am, i, j];
```

### 4.2 The homogeneous pair

```wl
hHom[am_List, c0_Real, c1_Real] := c0*obs0[am] + c1*obs1[am];

dHHom[am_List, c0_Real, c1_Real, i_Integer, j_Integer] :=
    c0*delObs0[am, i, j] + c1*delObs1[am, i, j];

conjObsHom = {obs0[#] &, obs1[#] &};
```

The driver is then called with `bt = 1.0`, `hHom[]`, `dHHom[]` and `conjObsHom`.

> **Pattern-matching trap.** `c0_Real` and `c1_Real` will not match an integer-valued field. A
> `dHHom` call that fails to match returns *unevaluated*, and the acceptance test's deliberate
> fourth branch treats an undetermined comparison as **reject** — so the chain silently stops
> moving instead of erroring. Wrap the field table in `N[...]`.

### 4.3 Coordinate helpers

```wl
cOf[beta_, nt0_]  := N[{beta, beta*nt0}];      (* (beta,nt0) -> c *)
physOf[c_List]    := {c[[1]], c[[2]]/c[[1]]};  (* c -> (beta,nt0) *)
cTable[betas_List, nt0s_List] := Flatten[Outer[cOf, N@betas, N@nt0s, 1], 1];
```

So `cTable[{1.0, 8.0}, {-1.6, 1.2}]` gives
`{{1., -1.6}, {1., 1.2}, {8., -12.8}, {8., 9.6}}`.

## 5. Equivalence, verified

Measured on 40 random 12-vertex graphs, $\beta \in \{0.5, 1, 3, 8\}$,
$n_{t0} \in \{-1.6, -0.5, 0, 0.7, 1.2\}$:

| check | result |
| --- | --- |
| $\beta\,h(am;n_{t0}) = h_{\text{hom}}(am;\beta,\beta n_{t0})$ | 800 comparisons, max abs. difference $2.3\times10^{-13}$ |
| $\beta\,dH(am;n_{t0},i,j) = dH_{\text{hom}}(\dots)$ | 800 comparisons, max abs. difference $2.8\times10^{-14}$ |
| **matched-seed chain**, 30 sweeps, 20 $(\beta,n_{t0})$ pairs | **final graph identical in every case** |
| reported energy | scales by exactly $\beta$ in every case |
| $\text{physOf}(\text{cOf}(\beta,n_{t0})) = (\beta,n_{t0})$ | exact to $10^{-12}$ |

The third row is the one that matters: it is not an agreement of formulas but of *trajectories*.
The two runs accept and reject the same moves in the same order and end on the same graph.

## 6. Units: what changes in the output

**This is the one thing that is genuinely different, and it is easy to miss.** Everything the
drivers record as an energy — `"energy"`, `minEnergy`, `minEstates`, `minEToBeat` — is now
$c\cdot O = \beta H_{\text{phys}}$, not $H_{\text{phys}}$.

Consequences:

1. **Within a replica it is consistent**, because $c$ is fixed there. Convergence tests,
   minimum-tracking and degeneracy collection all behave.
2. **Across replicas it is meaningless once $\beta$ varies.** A cold replica reports a
   numerically much lower energy for the *same physical state*. Do not compare `minEnergy`
   between rungs, and do not compare it against physical reference values such as the
   ground-state energies used for validating manifold runs.
3. **Chart column 5 becomes identically zero.** That column is
   `energy - field.conjObs`, the field-independent part, and in homogeneous form there is no
   field-independent part. The reconstruction `E(g) = col5 + g*triangles` therefore degenerates.

The replacement is better than what it replaces, because both components are recorded
explicitly. From the chart's conjugate-observable column (column 4):

$$E_{\text{phys}}(n_{t0}) \;=\; \texttt{conjObs[[1]]} \;+\; n_{t0}\cdot\texttt{conjObs[[2]]}$$

which is $\beta$-independent and valid at *any* field value, not just the one sampled.

## 7. Running it

1. **Bootstrap run + schedule.** Call the multi-field `GraphCTLSchedule`
   (`MCSims.wl:6321` with `delH`, `6562` without), which matches on
   `Length[externalFieldTable[[1]]] > 1`. A third overload, `MCSims.wl:6799`, builds a schedule
   from a previous run's histogram output instead of sampling afresh:

   ```wl
   GraphCTLSchedule[seedGraph, 1.0, hHom[], dHHom[],
                    cTable[betas, nt0s], conjObsHom, obs, NN, unlabeled, numReplicas]
   ```

   The bootstrap `cTable` must already span the $c$-region of interest — including varying
   $c_0$, or the schedule has no temperature information to work from.

2. **Read the schedule.** The return is `{{centers, edges, replicaLabels}, hist}`, where
   `centers` are the replica field vectors, `edges` the swap graph, and `hist` the
   multi-histogram output.

3. **Temper.** Feed the schedule to `GraphParallelTempering` — either the seed-graph form
   (`MCSims.wl:7880`) or the continue-from-replicas form (`MCSims.wl:8478`), both with
   `bt = 1.0`.

4. **Judge it by round trips, not swap acceptance.** Constant thermodynamic length equalises
   acceptance between neighbours, which is necessary and not sufficient. What matters is whether
   configurations actually travel from the hot end to the cold end and back. A healthy swap rate
   with no round trips means the cold replica is still frozen.

## 8. Practical limits

- **Two tempered components is the practical ceiling.** The schedule builder's diagnostic plot
  is `ContourPlot[PDF[dist,hVars], y]`, which takes exactly two ranges. $\beta$ plus one external
  field fits; adding a third axis requires editing those `Print` lines. The metric interpolation
  grid is $(1.4\,\texttt{numInterpolatingSamples}+1)^{\text{numVars}} = 29^{\text{numVars}}$ — the
  box is widened 20% each way at the unwidened step — which is a second reason not to reach for
  more.
- **The interpolation domain is the bounding box of the bootstrap table, while the reachable set
  is a fan.** The box corners correspond to low $\beta$ with large $|d|$ — parameters nothing was
  sampled at. Since 1.10.1 the replica density is masked by the MBAR **effective sample size**, so
  replicas are not placed there; but the interpolation grid is still spent on the whole box, so a
  wide $\beta$ range wastes most of it.
- **Narrow the $\beta$ window.** The fan's share of its bounding box shrinks as
  $\beta_{\max}/\beta_{\min}$ shrinks. Two or three overlapping windows also improve MBAR overlap,
  which is independently weakest at the cold end — precisely where the metric matters most.
- **Cost — no longer a limit.** 1.10.1 took the metric from 3 interpolant grids to 6 at
  `numVars = 2` (`numVars` first moments, `numVars(numVars+1)/2` second moments, one effective
  sample size), each an $O(K^2N)$ pass in the bootstrap table size $K$ — 7.2 hours at $K = 100$.
  All six are now filled in a single pass off one shared MBAR weight vector, in a fraction of a
  second. See **`MBARWeights.md`**. What remains is that a wide $\beta$ range still spends most of
  the grid on unsampled parameters, which is a quality argument, not a cost one.

## 9. Why 1.10.1 is a hard requirement

All four defects fixed in 1.10.1 live in the multi-field path and are invisible with one field —
so they had never been exercised before this recipe existed.

| | defect | effect on a $(\beta, d)$ schedule |
| --- | --- | --- |
| A | replica field-vector components could be permuted | $c_0 = \beta$ swapped with $c_1 = \beta d$, giving **negative $\beta$** on some replicas and not others |
| B | metric built from second moments, not covariances | with $O_0$ carrying the whole energy, the mean term dominates and placement stops tracking fluctuations |
| C | metric off-diagonals clipped to zero | the $O_0$–$O_1$ cross term is negative and was zeroed, discarding the coupling between the two axes |
| D | replicas placed on an extrapolated metric | replicas land in the box corners outside the fan |

B also affects single-field schedules, so **any schedule generated before 1.10.1 should be
regenerated** regardless of whether you use this recipe.

## 10. What this does not fix

Tempering in two dimensions is still tempering. If the transition of interest is first-order-like
— a bimodal energy histogram, hysteresis between hot and cold starts — then round trips die
however many axes the ladder has, because the tunnelling time between the two peaks grows
exponentially with system size. The constant-acceptance schedule does not address that: acceptance
between neighbouring rungs can look healthy while nothing traverses the ladder.

The instrument for that case is a density-of-states method — multicanonical or Wang–Landau
sampling, which walks in energy with weight $\propto 1/g(E)$, has no temperature and therefore no
freezing, and yields $\langle O\rangle_\beta$ at every $\beta$ from one run. **No such machinery
exists in ECGrav.** Before investing in it, run the cheap diagnostic: two chains at the target
$\beta$, one from the best known ground state and one from a random configuration. If their energy
averages meet, you are equilibrated and none of this is needed; if they do not, they bracket the
true value and the gap is a real bound.

## 11. See also

- `Theory.md` — the model, its observables and its energy functionals.
- `CHANGELOG.md` 1.10.1 — the four schedule-builder fixes, with their failure modes.
- `ReleaseNotes/v1.10.1.md` — the same, written for a reader upgrading.
- `Tests/MCSims.wlt:322` — `GraphSweepReplica-homogeneous-reparameterisation-is-identical`, the
  test this recipe rests on.
