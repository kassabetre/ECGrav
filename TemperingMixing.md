# Judging a tempering run by how it mixed

`TemperingMixing.wls` answers the question `HomogeneousHamiltonian.md` §7 step 4 poses and then
leaves open: not what the swap acceptance was, but whether configurations actually travelled from
the hot end of the ladder to the cold end and back.

```wl
Needs["ECGrav`"];
Get["/path/to/ECGrav/TemperingMixing.wls"];

mix = TemperingMixing`MixingData[res];        (* res = a GraphParallelTempering return *)
TemperingMixing`MixingReport[mix]
TemperingMixing`MixingPlots[mix]
```

It works on the return of `GraphParallelTempering` and on the four-element multi-histogram output
`GraphCTLSchedule` produces, which carries the same `histories` association.

---

## Part I — the physics

### 1. What parallel tempering is buying, and what can go wrong

A single Metropolis chain at large $\beta$ gets stuck: the moves that would carry it between
distant low-energy basins are rejected, so it samples one basin and reports its properties as if
they were the ensemble's. Parallel tempering runs a ladder of replicas at different $\beta$ and
periodically proposes exchanging their configurations. A configuration that is trapped at the cold
end can ride up the ladder to where the barriers are surmountable, cross, and come back down. The
cold replica then samples a genuinely different basin.

That is the whole mechanism, and it has a single point of failure: **the ladder has to actually
transport configurations.** Everything below is about measuring whether it did.

### 2. Why swap acceptance is not the answer

The obvious number to look at is the fraction of proposed exchanges that were accepted. It is
necessary — an acceptance of zero at some rung severs the ladder there — but it is emphatically
not sufficient, for a reason worth stating precisely.

Acceptance is a *local* quantity. Two neighbouring replicas at nearly the same $\beta$ exchange
almost every time, because their energy distributions overlap almost completely. A ladder of
sixteen replicas packed into a narrow $\beta$ window will show near-perfect acceptance everywhere
and transport nothing, because the whole ladder spans a temperature range too small to unfreeze
anything. Conversely a ladder with mediocre acceptance spread evenly across a wide range can move
configurations end to end very efficiently.

What acceptance *does* tell you is whether the schedule did its job. Constant thermodynamic length
places replicas so that acceptance is equal between every neighbouring pair; a large spread means
the metric the schedule was built on did not describe the region the replicas ended up in. That is
a diagnosis of the schedule, not of the sampling.

### 3. Round trips

The quantity that matters is the **round trip**: a configuration reaching the hottest rung and
returning to the coldest, or the reverse. It is the direct measurement of the mechanism in §1 —
one round trip is one opportunity for the cold replica to have swapped basins.

Two things follow.

**The round-trip rate sets the real sample size at the cold end.** Cold-replica measurements taken
between two round trips are correlated with each other in a way no autocorrelation time computed
*within* a replica can see: the replica's own chain may decorrelate in ten sweeps while every one
of those samples sits in the same basin. The number of independent cold configurations a run
produced is closer to its round-trip count than to `numsweeps/corrT`. A run with a healthy `corrT`
and zero round trips has one cold sample, however many measurements it wrote down.

**Zero round trips with healthy acceptance is the signature of a first-order-like transition.**
This is `HomogeneousHamiltonian.md` §10's warning: when the energy histogram is bimodal, the
tunnelling time between the peaks grows exponentially with system size, and adding rungs does not
help, because the bottleneck is not between two rungs — it is between two basins at the same
temperature. Neighbouring replicas still exchange happily; they just exchange configurations from
the same peak. Round trips are what distinguishes this case from a merely coarse ladder, and
nothing else in the output does.

### 4. Replica flow

Round trips give one number per configuration. The **replica flow** localises the problem along
the ladder.

Label each configuration with the end it visited most recently: *up* once it has touched the
coldest rung (it can only head toward hot next), *down* once it has touched the hottest. Then for
each rung $k$ count the visits made under each label and form

$$f(k) \;=\; \frac{n_{\uparrow}(k)}{n_{\uparrow}(k) + n_{\downarrow}(k)}.$$

By construction $f = 1$ at the coldest rung and $0$ at the hottest. In between, a ladder with no
bottleneck gives a smooth monotone descent — the two counter-propagating streams of configurations
are balanced everywhere and the walker performs an unbiased random walk in rung index. A **plateau**
in $f$ says configurations are arriving at that rung and turning back rather than passing through:
the diffusion is slow there and that is where extra rungs would buy the most.

This is the standard diagnostic from the replica-flow literature and it is more informative than the
round-trip count alone, because it says *where* the ladder is failing.

Properly $f$ is a histogram over *time* — a configuration that lingers at a rung should count for
as long as it lingered. That is what is computed when the run recorded its per-sweep occupancy
(§7); when the trajectories had to be replayed, the counts are per arrival instead, which
over-weights the high-degree rungs. On the four-replica cross-validation run of §7 the two agree to
within 0.013, but they are not the same quantity, and nothing guarantees they stay close.

**It assumes a linear ladder.** With one tempered parameter the rungs form a chain and $f$ is a
function of position along it. With two — $\beta$ and an external field — the rungs form a *graph*,
and ordering them by $\beta$ alone does not produce a path: two rungs at nearly the same $\beta$ can
sit at opposite ends of the field axis. A non-monotone $f$ on a two-dimensional schedule is
therefore **not by itself evidence of a bottleneck**. Read the round trips first, and treat $f$ as
suggestive until the ladder is a chain in the coordinate you plotted it against.

### 5. Occupancy

In a well-mixed run every configuration spends comparable time at every rung. `OccupancyMatrix`
gives the arrivals of each configuration at each rung; a run that mixed shows a broadly flat matrix,
and one that did not shows block structure — clusters of configurations confined to a range of
rungs. It is a coarser instrument than the flow, but it needs no labelling convention and it makes
a severed ladder obvious at a glance.

Whether that matrix is in **sweeps** or in **arrivals** depends on what the run recorded, and the
difference matters. Arrivals carry a structural bias: a rung with six swap-graph neighbours is
entered and left more often than one with four, whatever the sampling did. Sweeps do not.
`OccupancyKind[mix]` says which you have; §7 says when each is available.

### 6. How this relates to `corrT`

`corrT` is measured *within* a replica, from the autocorrelation of its observables, and it sets the
sweep interval between recorded measurements. It cannot see basin trapping, for the reason in §3.
The two diagnostics answer different questions and a run needs both to be healthy: `corrT` says
consecutive measurements at a rung are not the same configuration twice, round trips say the
configurations at the cold rung came from more than one basin.

`MixingReport` carries `corrT` and `corrTMeasured` through so the two can be read together.

> **`corrTMeasured` before the current build.** The overload of `GraphComputeCorrelationTime` that
> every external-field tempering run goes through never assigned the flag — it initialised it to
> `False` and left it there. Runs made before that was fixed report `corrTMeasured -> False` for
> every replica regardless of what was measured; read `corrTValues` instead, which was always
> correct. See `CHANGELOG.md`.

---

## Part II — the implementation

### 7. Two records, and which one you have

A run records the swap history twice over, and `MixingData` reports which it used as
`mix["source"]`.

**A tag column in the chart** — which configuration occupies each slot at each sweep, written
immediately after that sweep's `Swap[]`. This is the direct record: trajectories, dwell times and a
time-weighted flow all read straight off it. `MixingData` prefers it, reports
`"source" -> "chart"`, and `OccupancyKind` returns `"sweeps"`.

It sits at **column 7** in the external-field layout and **column 5** in the beta-schedule one, so
`MixingData` does not look at a fixed position. It tests the invariant: a tag column holds a key of
the chart in every row, and at every sweep the column read across slots is a **permutation** of the
keys, since a swap exchanges configurations and never creates or destroys one. Candidates are tried
from the last column backwards, so the detection survives further columns being added.

Testing the row *length* instead is not merely weaker, it is wrong, and the beta chart is the
counterexample: it splices nothing now, but before it was unflattened an eight-observable run had
eleven columns, and a `Length[row] >= 7` test read observable four as a tag and produced
`Missing[KeyAbsent, …]` where a slot index belonged.

**`histories`** — present in every run, and the only record in runs made before column 7 was added.
It requires the reconstruction of §§7.1–8. `"source" -> "replay"`, and occupancy falls back to
arrivals.

The two are cross-validated against each other on a run that has both: the walks agree element for
element, and so do the round-trip counts and the recovered swap graph. That is a check of two
independent implementations against two independent records in the same return, and it is worth
re-running after any change to either.

One asymmetry to know about. The chart's first row is written *after* the first sweep's `Swap[]`, so
the initial permutation is not in it; the convention that configuration $j$ starts at slot $j$
supplies it, which is exactly what the `histories` seeding at `MCSims.wl:7954` encodes. `MixingData`
prepends it to the walk — and to the walk only, since residence and the flow are counted over
recorded sweeps, where every sweep is a whole sweep. Without that prepend the chart walk is the
replay walk missing its first element, which is how this was found.

### 7.1 What `histories` actually holds

This is the thing to get right, and it is the opposite of what it looks like.

`histories` is keyed by replica index with `"externalField"`, `"swapAccept"`, `"swapTry"` and
`"history"`. The natural reading of `"history"` — a list of field values — is that it is the
trajectory of a replica through field space. It is not.

On an accepted swap between slots $k$ and $l$ the driver exchanges the two **states** and leaves
each slot's `"externalField"` where it is (`MCSims.wl:7999-8006`; the energies are recomputed
against the slot's own unchanged field). So **a slot is a fixed field value and what moves is the
configuration.** Then (`MCSims.wl:8011-8020`):

```wl
histories[[k, "history", 1]] = histories[[l, "history", -1]];
histories[[l, "history", 1]] = histories[[k, "history", -1]];
RotateLeft both
```

each slot's newest entry becomes the *other* slot's newest entry. Writing $x_k$ for the newest entry
at slot $k$, the update is a transposition of $x_k$ and $x_l$ — so the multiset $\{x_1,\dots,x_n\}$
is invariant, and since `history[[k,1]] === k` initially, the values are a permutation of the slots
at all times. Each value is a **tag** that moves exactly when a configuration moves. Therefore

> `histories[[k, "history"]]` is the ordered list of configuration **tags** that occupied slot $k$,
> each tag naming the slot its configuration started in.

It is the transpose of a trajectory, and the empirical check is decisive. A configuration moves only
along swap-graph edges, so under the trajectory reading *consecutive* entries would be graph-adjacent
with probability 1. Measured on the 16-replica run below: **33.1%** overall, and between 28.8% and
39.7% for every individual slot. The reconstruction of §7.2 is the positive control — its
trajectories are **100.0%** adjacent, as they must be.

The buffer is a rotating one of length `NN` initialised to `-1.0`, and the driver trims the unused
padding at exit (`MCSims.wl:8127-8129` for tempering, `4987-4991` for the histogram path), so a
returned history has length `swapAccept + 1` and contains no `-1.0`. Check anyway; never treat a
`-1.0` as a field value.

### 7.2 Replaying the run

Per-slot tag sequences are not directly usable: round trips are a property of a *configuration's*
path, and the buffers carry no timestamps and are not aligned with each other — only the two slots
involved in a swap advance.

They are nonetheless enough to reconstruct the run exactly.

**The admissibility rule.** An event at slots $k, l$ can be the next one only if each slot's next
unread entry is the other slot's current occupant:

$$\texttt{tags}[k][p_k + 1] = \texttt{at}[l] \quad\text{and}\quad \texttt{tags}[l][p_l + 1] = \texttt{at}[k].$$

**The disjointness lemma.** Any two admissible events are slot-disjoint. Suppose two shared slot
$k$. Both would require $\texttt{tags}[k][p_k+1]$ to equal the current occupant of their other slot;
but that entry names one tag, and occupants are distinct because the tags form a permutation. So the
two other slots would have to be the same slot, making the events identical. $\square$

Disjoint transpositions commute, so **the replay is order-independent**: whichever admissible event
is taken first, the per-configuration trajectories come out the same. Verified on the run below —
identical under first, last and random tie-breaking, with all 2402 swaps consumed and every buffer
read to its end.

The loop is therefore just:

```
pointers all at 1;  at[k] = k;  traj[k] = {k}
repeat
    find any admissible pair (k,l)          (* all admissible pairs commute *)
    swap at[k], at[l];  advance both pointers
    append k to traj of the tag now at k, and l to the tag now at l
until none admissible
```

`MixingData` reports `"complete"`, which is `True` only when every pointer reached the end of its
buffer. A `False` means the histories are not consistent with a pure exchange of states — a resumed
or edited run — and nothing downstream should be trusted.

**The swap graph is not an input.** The admissibility rule alone selects the events; the graph falls
out as the set of pairs that fired, returned as `"swapGraph"`. On the run below it recovered 41
edges, exactly the schedule's. This matters because the `GraphParallelTempering` return does not
carry the schedule, so requiring it would mean carrying the `GraphCTLSchedule` output around.

**Cost** is one pass over the admissible pairs per event: $O(E \cdot S)$ for $S$ swaps, 1.1 s for
2402 swaps over 16 replicas.

### 8. Rungs, ends and counting

Rungs are numbered by increasing $\beta$ — in homogeneous $c\cdot O$ form that is $c_0$, the
component the ladder is a ladder in — so rung 1 is the **hottest** and rung $n$ the coldest.
`mix["order"]` holds the permutation; the raw arrays stay in the chart's own key order.

A round trip is counted by walking each trajectory with a single label: touching the coldest rung
sets the label, touching the hottest increments the count if the label is set. `RoundTrips` reports
both directions (`"perConfiguration"` for cold→hot→cold and `"perConfigurationReverse"`); on a run
that mixed they agree to within one per configuration, since they count the same traversals offset
by half a period.

`"meanRoundTripSweeps"` divides the total sweep count by round trips per configuration. It is in
sweeps, not accepted swaps, because sweeps are what a production run is budgeted in.

### 9. Worked example

The 16-replica, 1000-sweep run in `PTSim…numEvents1000.mx`, on the schedule from the matching
`CTLSim…numEvents800.mx`:

| | |
| --- | --- |
| accepted swaps replayed | 2402 of 2402, replay complete |
| swap graph recovered | 41 edges, exactly the schedule's |
| round trips | **87** |
| configurations with at least one | **16 of 16** (2 to 8 each) |
| mean round-trip time | 184 sweeps, or 27.6 accepted swaps |
| rungs visited per configuration | 16 of 16, all of them |
| swap acceptance | 0.159 – 0.606, a **3.8×** spread |

The verdict is that the ladder transported configurations well: every configuration made the full
traverse several times, so the cold rung saw on the order of five independent basins rather than
one. The acceptance spread is the weak point and it is a statement about the **schedule** (§2): this
one was built over the full bounding box of the bootstrap table, so seven of its rungs sit outside
the sampled fan — see `HomogeneousHamiltonian.md` §8.

The flow $f$ came out non-monotone, and per §4 that is expected here rather than diagnostic: the
schedule is two-dimensional, with rung 4 at $(\beta, n_{t0}) = (0.857, 1.61)$ and rung 5 at
$(0.898, -2.22)$.

### 10. API

| | |
| --- | --- |
| `MixingData[res]` | replays the run; returns the bundle everything else takes |
| `RoundTrips[mix]` | per-configuration counts, totals, rate, visits to each end |
| `ReplicaFlow[mix]` | $f$ at each rung, in the chart's key order |
| `OccupancyMatrix[mix]` | `m[[rung, configuration]]`, arrivals |
| `MixingReport[mix]` | the summary table above |
| `MixingPlots[mix]` | `"RoundTrips"`, `"Acceptance"`, `"Flow"`, `"Occupancy"`, `"Trajectory"` |
| `TrajectoryPlot[mix, j]` | configuration `j`'s rung against swap number |

`MixingData` returns, among others: `"trajectories"` (per configuration, the rungs it occupied in
order), `"order"` (rungs sorted hot→cold), `"betas"`, `"nt0"`, `"swapGraph"`, `"acceptedSwaps"`,
`"sweeps"`, `"complete"`, and the `corrT` fields passed through from `replicas`.

### 11. What this does not do

- **Dwell times in a run made before column 7 existed.** The histories record the order of moves,
  not their duration, and nothing recovers it after the fact — the replay is exact about *what*
  happened and silent about *when*. Re-run to get residence-weighted occupancy and a time-weighted
  flow; there is no way to retrofit them.
- **Per-edge acceptance.** `swapAccept`/`swapTry` are per replica, so a single bad edge is visible
  only as depressed acceptance at both its ends.
- **A flow diagnostic for a genuinely two-dimensional ladder.** §4's $f$ wants a chain.

### 12. See also

- `HomogeneousHamiltonian.md` §7 step 4 (judge it by round trips), §8 (why this schedule's rungs
  landed where they did), §10 (first-order transitions, which no ladder fixes).
- `TemperingPlots.md` — the thermodynamic surfaces from the same run.
- `CHANGELOG.md` — the `corrTMeasured` fix referred to in §6.
