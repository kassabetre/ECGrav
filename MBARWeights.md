# One weight vector: MBAR reweighting for the CTL metric

> **What this is.** A specification of a **proposed** rewrite of the MBAR reweighting that the
> constant-thermodynamic-length schedule builder runs on. It records the algebra the rewrite
> rests on, the algorithm, what was measured, and what it does and does not change about the
> schedules produced. Companion to `HomogeneousHamiltonian.md`, which is the recipe this cost
> problem was found under.
>
> **Status: implemented.** Landed in `Kernel/Subpackages/MCSims.wl` as `MBARWeightBasis`,
> `MBARWeights`, `MBARWeightedMean`, `MBARMomentsAt` and `MBARESSAt`, with the three public MBAR
> entry points and all three multi-field `GraphCTLSchedule` copies rewritten on top of them.
> 193 tests green. §1 describes the nested-`Sum` form it replaced.
>
> **Why it matters.** The metric build is the dominant cost of a multi-field schedule, and it
> grows as the *square* of the bootstrap table. At one external field it is a few seconds and
> nobody noticed. At two it reaches hours.

Notation: $\beta$ is the package's fixed `betaFixed`/`bt`; $x \in \mathbb{R}^d$ the target field
at which the metric is wanted; $d$ = `numVars`, the number of tempered components. The bootstrap
run visited $K$ field points $h_1,\dots,h_K$ (the keys of `minusbetaF`), drawing $n_j$ samples at
$h_j$. Samples are indexed by one flat index $s = 1,\dots,S$ with $S = \sum_j n_j = KN$, running
over all field points concatenated in key order. $O_s \in \mathbb{R}^d$ is the conjugate-observable
vector of sample $s$ (chart column 4), $F_j$ is `minusBetaF[h_j]`, and $\odot$ is the elementwise
product.

---

## 1. The problem

`GraphCTLSchedule`'s multi-field overloads (`MCSims.wl:6244`, `6469`, `6690`) build the
thermodynamic metric by interpolating MBAR-extrapolated moments over a grid. Three things
multiply together, and each is larger than it looks.

**The grid is $29^d$, not $20^d$.** The `PrintTemporary` at `MCSims.wl:6309` announces
`numInterpolatingSamples^numVars` $= 20^d$ points. But the iterator widens the box by 20% on each
side while leaving the step at the *unwidened* range over 20 (`MCSims.wl:6320`):

```wl
{hMin - 0.2*range, hMax + 0.2*range, range/numInterpolatingSamples}
```

so each axis carries $1.4\times20+1 = 29$ points. At $d = 2$ that is **841**, not 400 — the
printed figure understates by $1.4^d$.

**There are six grids at $d = 2$, not one.** Three second moments (`MCSims.wl:6316`), two first
moments (`6332`), one effective sample size (`6346`) — $d(d{+}1)/2 + d + 1$ in general. So
$5\times841$ calls to `ExtrapolatedExpectationValue` and $841$ to `MBAREffectiveSampleSize`.

**Each call is $O(K^2N)$.** `NegativeBetaTimesFreeEnergy` (`MCSims.wl:485`) builds a $K\!\times\!N$
table of `LogSumExp`s over $K$ terms each; `ExtrapolatedExpectationValue` (`MCSims.wl:551`) then
calls it and repeats the same double sum. All of it is interpreted `Sum[]` over Association
`[[Key[j], s]]` lookups — nothing is a packed array. Measured:

| $K$ | $N$ | `ExtrapolatedExpectationValue` | `MBAREffectiveSampleSize` |
| --- | --- | --- | --- |
| 16 | 50 | 0.0740 s | 0.0415 s |
| 25 | 100 | 0.3471 s | 0.1930 s |
| 49 | 100 | 1.3367 s | 0.7239 s |
| 100 | 100 | 5.4246 s | 2.9471 s |
| 100 | 200 | 10.9331 s | 6.1405 s |

Clean $K^2N$. The full $d=2$ metric build is $5\cdot841\cdot t_{\text{EEV}} + 841\cdot t_{\text{ESS}}$:

| $K$, $N$ | serial | wall, ~10 kernels |
| --- | --- | --- |
| 16, 100 | 692 s | ~70 s |
| 25, 100 | 1 783 s | ~3 min |
| 49, 100 | 6 279 s | ~11 min |
| 100, 100 | 26 190 s | ~44 min |
| 100, 200 | 51 661 s | ~80 min |

For contrast, $d = 1$ is 29 points over 3 grids at $K\approx12$: about **5 seconds**.

One further trap while reading a stalled run: `ParallelTable` defaults to `Method -> Automatic`,
which is coarse-grained. The 29 outer rows are split into ~3-row chunks, one per kernel, and
nothing returns until a chunk completes — so the front end's
`"Distributing definitions, sent: 0, recvd: 0, ... 0 bytes each"` panel sits unchanged for the
whole build. It is stale, not deadlocked; `0 bytes` means the subkernels already had the
definitions.

## 2. The one quantity everything is built from

Define the **log MBAR denominator** of sample $s$ at target $x$:

$$D_s(x) \;=\; \log \sum_{j=1}^{K} n_j\, e^{-F_j}\, e^{\beta (x - h_j)\cdot O_s}$$

Every one of the three shipped functions is $D_s$ in disguise:

| function | where | in terms of $D_s$ |
| --- | --- | --- |
| `NegativeBetaTimesFreeEnergy` | `MCSims.wl:485` | $\log \sum_s e^{-D_s(x)}$ |
| `ExtrapolatedExpectationValue` | `MCSims.wl:551` | $e^{-\text{NBTFE}}\sum_s v_s\,e^{-D_s(x)}$ |
| `MBAREffectiveSampleSize` | `MCSims.wl:622` | its `logu` **is** $-D_s(x)$ |

In the second, `v` is `measuredObservableValues`; the inner
`Sum[nValues*Exp[-minusBetaF]*Exp[betaFixed*(target-j).conj], {j,...}]` is exactly
$e^{D_s(x)}$.

## 3. Three reductions

### 3.1 The normalisation cancels, so the free energy is never needed

Write $u_s = e^{-D_s(x)}$. Then $\text{NBTFE}(x) = \log\sum_s u_s$, and

$$\langle v\rangle_x \;=\; e^{-\text{NBTFE}(x)} \sum_s v_s\, e^{-D_s(x)} \;=\; \frac{\sum_s v_s u_s}{\sum_s u_s} \;=\; \frac{v\cdot u}{\text{Total}[u]}$$

An ordinary **weighted average**. The free energy was only ever the normaliser of $u$ and it
divides straight out. The shipped `ExtrapolatedExpectationValue` computes it anyway, at
$O(K^2N)$, and then rebuilds the same denominators a second time — a factor of two thrown away
before anything else.

The effective sample size is the same vector. Unwinding `MCSims.wl:645-648`, with
$A = \text{LogSumExp}(\log u)$ and $B = \log\sum_s u_s^2$, the returned `Exp[2A - B]` is the Kish
form

$$\text{ESS}(x) \;=\; \frac{\big(\sum_s u_s\big)^2}{\sum_s u_s^2}$$

### 3.2 One weight vector serves all six grids

$u_s$ depends on the target $x$ and on nothing else — it has no knowledge of which observable is
about to be averaged against it. But the six grids at `MCSims.wl:6316-6352` each rebuild the
identical denominators from scratch. Given $u$:

$$E[O_i] = \frac{O_{\cdot i}\cdot u}{\Sigma u}, \qquad
E[O_iO_j] = \frac{(O_{\cdot i} \odot O_{\cdot j})\cdot u}{\Sigma u}, \qquad
\text{ESS} = \frac{(\Sigma u)^2}{\Sigma u^2}$$

Six dot products against one vector. Factor of six.

### 3.3 The $K$-sum does not depend on the target

This is the reduction that carries the result. Expand the exponent inside $D_s$:

$$D_s(x) = \log \sum_j \exp\big(\log n_j - F_j + \beta\, x\cdot O_s - \beta\, h_j\cdot O_s\big)$$

The term $\beta\,x\cdot O_s$ carries **no $j$**, so it factors out of the sum:

$$\boxed{\;D_s(x) \;=\; \underbrace{\beta\, x\cdot O_s}_{b_s(x)} \;+\; \underbrace{\log \sum_j \exp\big(\log n_j - F_j - \beta\, h_j\cdot O_s\big)}_{L_s}\;}$$

$L_s$ contains no $x$. **The entire sum over the $K$ bootstrap states — the only quadratic part of
the whole computation — is evaluated once for the entire grid**, not once per grid point per
observable. What remains per grid point is $b_s(x) = \beta\,x\cdot O_s$: a single matrix–vector
product.

In matrix form, with $\mathbf{O}$ the $S\times d$ packed observable matrix and $\mathbf{H}$ the
$K\times d$ key matrix:

$$\mathbf{A} = \beta\,\mathbf{O}\mathbf{H}^{\mathsf{T}} \in \mathbb{R}^{S\times K},
\qquad L_s = \operatorname*{LogSumExp}_{j}\big(\log n_j - F_j - \mathbf{A}_{sj}\big),
\qquad b(x) = \beta\,\mathbf{O}x$$

## 4. Numerical stability, which improves for free

$D_s$ spans hundreds in the exponent once $x$ moves away from the sampled fields — which is
precisely the regime the interpolation grid is built for, since it covers the bounding box and the
reachable set is a fan inside it.

Because *every* output in §3 is a ratio of $u$, the vector may be rescaled freely: both
$(v\cdot u)/\Sigma u$ and $(\Sigma u)^2/\Sigma u^2$ are invariant under $u \to cu$. So shift
$\log u$ by its own maximum before exponentiating. Every entry then lies in $(0,1]$ and cannot
overflow, and `Clip[..., {-700., 0.}]` keeps the small end off the subnormal boundary — which is
where the cost cliff actually is, not at large arguments.

Two notes on the existing code. `ExtrapolatedExpectationValue` does none of this: it divides by
raw denominators, which is why it emits `General::munfl` under any benchmark that reaches outside
the sampled box. And the package's own `LogSumExp` (`MCSims.wl:45`) shifts by the **mean**, which
does not survive this spread — `MBAREffectiveSampleSize` already carries a local max-shift for
that reason, and the form below does the same.

## 5. The algorithm

```wl
(* ---------- once per schedule build ---------- *)
mbarPrecompute[bt_, minusBetaF_Association, conj_Association] :=
  Module[{keys = Keys[minusBetaF], Hm, Om, logn, F, A, L},
    Hm   = N[keys];                                       (* K x d *)
    Om   = N @ Join @@ Lookup[conj, Key /@ keys];         (* S x d, packed *)
    logn = Log[N[Length /@ Lookup[conj, Key /@ keys]]];   (* K *)
    F    = N @ Values[minusBetaF];                        (* K *)
    A    = bt * (Om . Transpose[Hm]);                     (* S x K, x-independent *)
    L    = Module[{t = ConstantArray[logn - F, Length[A]] - A, m},
             m = Max /@ t;
             m + Log[Total[Exp[Clip[t - m, {-700., 0.}]], {2}]]];
    <|"Om" -> Om, "L" -> L, "bt" -> bt|>];

(* ---------- once per grid point ---------- *)
mbarWeights[pc_, x_List] :=
  With[{lu = -(pc["bt"] * (pc["Om"] . N[x]) + pc["L"])},
    Exp[Clip[lu - Max[lu], {-700., 0.}]]];

(* ---------- all six quantities from that one vector ---------- *)
mbarMean[u_, v_]  := (v . u) / Total[u];
mbarESS[u_]       := Total[u]^2 / Total[u^2];
```

The observable vector `v` must be flattened in the **same key order** as `Om`, i.e.
`Join @@ Lookup[obsAssoc, Key /@ Keys[minusBetaF]]`. Second moments are
`mbarMean[u, Om[[All,i]] * Om[[All,j]]]`; first moments are `mbarMean[u, Om[[All,i]]]`.

## 6. Cost

Per grid point the shipped path is $O(K^2Nd) = O(SKd)$ and runs $d(d{+}1)/2 + d + 1$ times; the
vectorised path is $O(Sd)$ once, against a one-time $O(SK)$ setup. With $G$ grid points and
$G_{\text{grids}}$ grids the algorithmic ratio is

$$\frac{G_{\text{grids}}\cdot G\cdot SK}{SK + G\cdot S} \;=\; \frac{G_{\text{grids}}\, G\, K}{K + G}$$

which at $d=2$, $G=841$, $K=100$ is about **536×**. The remaining factor comes from replacing
interpreted `Sum` over Association lookups with `Dot` on packed reals. Measured, for the complete
841-point six-quantity build:

| $K$, $N$ | shipped (serial) | vectorised | ratio |
| --- | --- | --- | --- |
| 25, 100 | 1 783 s | 0.03 s | $5.9\times10^4$ |
| 49, 100 | 6 279 s | 0.05 s | $1.3\times10^5$ |
| 100, 100 | 26 190 s | 0.10 s | $2.6\times10^5$ |
| 100, 200 | 51 661 s | 0.15 s | $3.4\times10^5$ |
| 225, 200 | 257 582 s | 0.34 s | $7.6\times10^5$ |

$536 \times \approx 500 \approx 2.6\times10^5$ — the two factors account for the measurement.

**Measured after landing.** The same six quantities at 25 grid points, serial, through the public
API, against the pre-change file checked out beside it:

| $K$, $N$ | before | after | ratio |
| --- | --- | --- | --- |
| 16, 120 | 24.18 s | 0.072 s | 336× |
| 25, 100 | 47.94 s | 0.078 s | 614× |
| 49, 100 | 180.76 s | 0.269 s | 672× |
| 100, 100 | 767.50 s | 1.215 s | 632× |

That is §3.1 and the packed-array change alone — this route rebuilds the basis at every call, so
§3.3 contributes nothing to it. Scaled to the real 841-point grid the before column is 813 s at
$K=16$ and **25 819 s (7.2 h)** at $K=100$, against the 26 190 s the model in §1 predicts. The
schedule builder amortises the basis across the grid on top of this, which is where the remaining
factor of ~400 in the table above comes from.

**Memory.** $\mathbf{A}$ is $S\times K$ machine reals: 81 MB at $K=225$, $N=200$. It is transient —
only $L$, of length $S$, survives — so build it in row chunks if $S K$ grows beyond comfort. $L$
itself is negligible.

## 7. Verification

Against the shipped `ExtrapolatedExpectationValue` and `MBAREffectiveSampleSize`, over
$K \in \{4,9,25,49\}$ at $N=50$ and five target fields including two well outside the sampled box:

| $K$ | max rel. difference, EEV | max rel. difference, ESS |
| --- | --- | --- |
| 4 | $2.19\times10^{-14}$ | $8.78\times10^{-15}$ |
| 9 | $3.51\times10^{-14}$ | $1.32\times10^{-14}$ |
| 25 | $1.16\times10^{-14}$ | $2.01\times10^{-14}$ |
| 49 | $1.40\times10^{-14}$ | $1.11\times10^{-14}$ |

Machine epsilon. A second, independent check uses an exactly-solvable fixture: take the density of
states to be $O \sim \mathrm{MVN}(0,\Sigma)$ with reduced potential $c\cdot O$, so that samples at
field $c$ are $\mathrm{MVN}(-\Sigma c, \Sigma)$ and $\texttt{minusBetaF} = \tfrac12 c^{\mathsf T}\Sigma c$ in closed
form, with $E[O](x) = -\Sigma x$ exactly. This is the exponential family the CTL assumes, so it
tests the reweighting rather than only the refactor. MBAR reproduces the exact mean to 0.4% at
$K=20$, $N=400$.

The §5 listing was then run **verbatim** — copied out of this file unedited — against the shipped
functions on that fixture at $K=25$, $N=80$, covering all three quantity types and four targets:

| target $x$ | $E[O_1]$ | $E[O_1O_2]$ | ESS | max rel. difference |
| --- | --- | --- | --- | --- |
| $(2, -1)$ | $-6.46981$ | $-37.59571$ | 292.50 | $1.5\times10^{-15}$ |
| $(0.5, -1)$ | $-0.52868$ | $-3.14017$ | 95.87 | $1.9\times10^{-15}$ |
| $(3, 3)$ | $-16.04105$ | $511.00210$ | 90.16 | $5.6\times10^{-15}$ |
| $(-1, 6)$ | $-12.85510$ | $448.17483$ | **2.01** | $2.4\times10^{-15}$ |

The last row is the case that matters: an ESS of 2.01 means two samples carry essentially all the
weight, which is the deep-extrapolation corner the log-space shift exists for. The shipped form
agrees there — it simply pays `General::munfl` to get the same answer.

The fixture is four lines, and needs no MCMC run:

```wl
Sig    = {{4.0, 1.5}, {1.5, 9.0}};
fields = Flatten[Outer[{#1, #1*#2} &, Subdivide[0.5, 3., 4], Subdivide[-1., 1., 4], 1], 1];
conj   = Association @ Table[c -> RandomVariate[MultinormalDistribution[-Sig . c, Sig], 80], {c, fields}];
mbf    = Association @ Table[c -> 0.5*(c . Sig . c), {c, fields}];
```

There is no `.wlt` coverage yet, because nothing is landed. Anything landed from this document
should pin the identity $\langle v\rangle = (v\cdot u)/\Sigma u$ against the shipped
nested-`Sum` form as a `VerificationTest`, over this fixture, including a target outside the
bootstrap box.

## 8. Implementation map

| what | where | change |
| --- | --- | --- |
| `NegativeBetaTimesFreeEnergy`, list overload | `MCSims.wl:485` | keep the signature; body becomes `Log[Total[u]]` on unshifted weights, or leave it and stop calling it from `ExtrapolatedExpectationValue` |
| `ExtrapolatedExpectationValue`, list overload | `MCSims.wl:551` | body becomes $(v\cdot u)/\Sigma u$; drops the `NegativeBetaTimesFreeEnergy` call entirely |
| `MBAREffectiveSampleSize` | `MCSims.wl:622` | body becomes $(\Sigma u)^2/\Sigma u^2$ |
| the six grids | `MCSims.wl:6316-6352`, `6535-6571`, `6732-6767` | collapse to one `mbarPrecompute` plus a single `ParallelTable` over grid points returning all six values |
| `essInterpolationFn` | `MCSims.wl:6355`, `6574`, `6770` | delete — see §9 |

The scalar `bf_Real` overloads (`MCSims.wl:448`, `525`) are the beta-tempering path and take the
same treatment with $O_s$ the energies and $h_j$ the betas, but they were never the bottleneck and
can be left alone.

All three schedule-builder copies are near-identical; whatever is done to one must be done to all
three, and 193 tests must stay green.

## 9. Consequence: the ESS interpolant goes

`essInterpolationFn` existed only because an exact ESS call used to cost seconds. It costs
**0.107 ms** once §5 is in place, and `rho` is evaluated on the $21^d$ `densityPoints` grid —
**441 exact calls, about 0.05 s for the whole mask**. So it is now removed, and `rho` calls
`MBARESSAt` directly.

**This is not a bug fix, and an earlier draft of this document wrongly implied it was.** The
interpolant was harmless in the shipped code, for a reason worth recording: `densityPoints` steps
by the same $(h_{\max}-h_{\min})/\texttt{numInterpolatingSamples}$ as the interpolation grid and
merely starts four steps inside it, so its points coincide with interpolation **nodes** to
$9\times10^{-16}$, and `Interpolation` reproduces its own node values exactly. Measured across
$K \in \{16,25,49,100\}$: **0 of 441 mask verdicts differ** between the interpolant and exact ESS.
Schedules do not move.

What makes removal worth doing is that this correctness rested entirely on a coincidence between
two grid definitions written in different places, with nothing recording the dependency. Evaluated
anywhere *off* node — which any change to either grid would cause — the interpolant is poor.
Against exact ESS at 4000 off-grid points, $K=96$, $N=100$:

| ESS surrogate | negative values | mask errors / 4000 | false masks |
| --- | --- | --- | --- |
| `Interpolation`, order 3 (what shipped) | 43 | 246 | 140 |
| order 1 | 0 | 229 | 131 |
| interpolate $\log$ ESS | 0 | 232 | 141 |
| $\log$ ESS, order 1 | 0 | 223 | 152 |
| **exact** | **0** | **0** | **0** |

The ESS surface runs from 1 to 604 across a $29\times29$ grid, so this is a **resolution** limit,
not an interpolation-order one — lowering the order barely moves it. A *false mask* is a point
with genuine MBAR support where `rho` would be zeroed anyway. Order 3 additionally returns
**negative** ESS, meaningless for a Kish count bounded below by 1.

## 10. What this does not change

- **The schedules themselves**, up to floating-point noise. §3 is an algebraic identity, verified
  to $2\times10^{-14}$; the metric, the density $\rho$, the K-means clustering and the swap graph
  are untouched. §9 does not change them either, as measured there — it removes a hidden coupling
  rather than a wrong answer.
- **The $K^2$ in the setup.** $\mathbf{A}$ is still $S\times K$, so the *one-time* cost still grows
  with the bootstrap table. It is one matmul rather than $6G$ interpreted double sums, which is the
  whole point, but a very large bootstrap table is still a large bootstrap table.
- **`essFloor` itself** (`MCSims.wl:6361`). It is `essMinFraction*Total[Length/@chartConjugateField]`
  $= 0.01KN$, and its meaning drifts with $K$: at $K=20$ it keeps 9.8% of the box where a
  $K$-independent $0.5N$ keeps 5.0%, while at $K=400$ it keeps 23.5% against 36.6%. Refining the
  bootstrap grid silently redefines "supported" by about 1.5×. Worth revisiting, but it is a
  separate question, and measurement does **not** support the stronger claim that the floor goes
  degenerate at large $K$: the kept fraction is stable at 22–24% from $K=48$ to $K=400$, because
  ESS grows with $K$ ($\text{ESS}_{\max}\approx 5.5K$) as neighbouring replicas start to overlap.
  As a sanity check, the fan occupies $88.2/307.3 = 28.7\%$ of its bounding box analytically, and
  the mask keeps 22–24% — it is doing roughly the right thing.
- **The grid size.** $29^d$ rather than the advertised $20^d$ is unaffected by any of this. Fixing
  the `PrintTemporary` at `MCSims.wl:6309` to report real grid dimensions is a one-line change
  worth making at the same time, and `numInterpolatingSamples` is a hard-coded `Module` local
  (`MCSims.wl:6277`, `6496`, `6709`) that could reasonably become an option.
- **Sampling quality.** Nothing here helps a frozen chain at high $\beta$. If the bootstrap run
  reports `GraphEquilibriate::noconv` or a `GraphComputeCorrelationTime::shortrun` with a
  correlation time comparable to the run length, the metric is being built from data that is not
  from the target distribution, and a faster wrong answer is still wrong.

## 11. See also

- `HomogeneousHamiltonian.md` — the two-component tempering recipe this was found under, and §8
  there for the practical limits of $d = 2$.
- `Theory.md` — the model, observables, and the thermodynamic-length argument.
- `CHANGELOG.md`, 1.10.1 — the four multi-field CTL correctness fixes. Two of them
  (`25f91a8`, `140031a`) are what took the metric from three grids to six; they are correct and
  should stay, and this document is how to afford them.
