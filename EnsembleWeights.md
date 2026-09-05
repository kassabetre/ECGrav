# Ensemble Weights for Random Pure Complexes at Fixed Facet Count

How the vertex-labeled, facet-labeled and unlabeled ensembles of pure simplicial complexes are
related to one another, and what that implies for the distributions of topological and combinatorial
statistics. Written for P1 (see `ExpansionRoadmap.md`, Track A).

The headline: **the three ensembles are related exactly, by automorphism group orders.** The
convergence of the facet-labeled and unlabeled ensembles as the facet count grows is not merely an
empirical observation to be tested — it is a theorem with a rate that the package can compute from
two counter calls, with no Monte Carlo at all.

_Written 2026-09-04. Every numeric claim below was re-run against the shipped code before being
recorded._

---

## 1. The three ensembles

Fix a **purity** `p` (facets are `p`-subsets) and a **facet order** `M`. Throughout, the vertex
count `n` is **free**: it is determined by the complex rather than imposed, and ranges over
`p <= n <= pM` with `C(n,p) >= M`.

A pure complex is a set of `M` pairwise distinct `p`-subsets whose union is the whole vertex set
(the covering condition is what ties `n` to the complex).

| | ensemble | objects | package |
|---|---|---|---|
| `S_V(p,M)` | vertex-labeled | `M`-element sets of `p`-subsets of `[n]` covering `[n]`, over all `n` | `NumVertexLabeledPureComplexes[p, M]` |
| `S_F(p,M)` | facet-labeled | `M` **labeled**, pairwise distinct `p`-subsets covering the vertex set, up to vertex permutation | `NumFacetLabeledPureComplexes[p, M]` |
| `S_U(p,M)` | unlabeled | isomorphism classes | `NumUnlabeledPureComplexes[p, M]` |

Each is sampled uniformly by the corresponding exact, rejection-free sampler
(`RandomVertexLabeledPureSimplicialComplex[{p,M}, k]` and the two `RandomUniform…` samplers), with
the vertex count left free by omitting it.

## 2. Two automorphism groups, not one

For a complex `C` on vertex set `V(C)` with `n(C) = |V(C)|`:

- **`Γ(C) <= Sym(V(C))`** — the **vertex** automorphism group: permutations of the vertices carrying
  the facet set to itself. Package: `PureComplexAutomorphismGroupOrder`.
- The action of `Γ(C)` on the `M` facets is a homomorphism `ρ : Γ(C) → S_M`. Its image
  **`Γ_F(C) = ρ(Γ(C))`** is the **facet** automorphism group. Package:
  `PureComplexFacetAutomorphismGroupOrder`.

`ker ρ` is the subgroup fixing every facet *setwise*, and it is **not** generally trivial:

$$|Γ(C)| \;=\; |\ker ρ(C)| \cdot |Γ_F(C)|.$$

The single triangle `{{1,2,3}}` is the cleanest example: `|Γ| = 6` (all of `S_3`) while `|Γ_F| = 1`,
since there is only one facet to permute. Two disjoint triangles give `|Γ| = 72` and `|Γ_F| = 2`,
so `|ker ρ| = 36 = 3! · 3!`. Both are in the test suite.

**The distinction matters, because the two labelings see different groups.** This is the whole
content of what follows.

## 3. Orbit counts

**Lemma.** Let `C` be an isomorphism class of pure complexes with purity `p`, facet order `M`, and
`n = n(C)` vertices. Then the class contains

- **`n! / |Γ(C)|`** distinct vertex-labeled complexes, and
- **`M! / |Γ_F(C)|`** distinct facet-labeled complexes.

*Proof.* For the first, `S_n` acts on vertex-labeled complexes with vertex set `[n]`; the stabilizer
of `C` is exactly `Γ(C)`, so orbit–stabilizer gives `n!/|Γ(C)|`.

For the second, a facet-labeled complex is a bijection `ℓ : [M] → F(C)` taken up to vertex
permutation. There are `M!` such bijections for a fixed representative `C`, and `Γ(C)` acts on them
through `ρ`. A labeling is fixed by `σ` exactly when `ρ(σ) = id`, i.e. when `σ ∈ ker ρ`, so every
orbit has size `|Γ(C)|/|\ker ρ(C)| = |Γ_F(C)|`, and the number of orbits is `M!/|Γ_F(C)|`. ∎

Summing over classes:

$$|S_V| = \sum_C \frac{n(C)!}{|Γ(C)|}, \qquad |S_F| = \sum_C \frac{M!}{|Γ_F(C)|}, \qquad |S_U| = \sum_C 1 .$$

## 4. The induced weights — the central table

Push each uniform ensemble forward onto isomorphism classes. Up to normalization, class `C` receives:

| ensemble | weight on class `C` | behaviour |
|---|---|---|
| `S_U` | `1` | uniform by definition |
| `S_F` | `1 / |Γ_F(C)|` | **bounded in `[1/M!, 1]`**, and `→ 1` as complexes lose facet-moving symmetry |
| `S_V` | `n(C)! / |Γ(C)|` | **factorial in a varying `n`** — unbounded across the class set |

Everything else in this document follows from those three rows.

Note what the middle row does *not* depend on: a complex with an enormous vertex automorphism group
but a trivial facet action — the single facet, again — is weighted **identically** in `S_F` and
`S_U`. The facet-labeled/unlabeled discrepancy is governed by `Γ_F` alone, never by `Γ`.

## 5. Exact identities

Dividing the sums of §3:

> $$\boxed{\;\frac{|S_F(p,M)|}{M!\,|S_U(p,M)|} \;=\; \mathbb{E}_U\!\left[\frac{1}{|Γ_F|}\right]\;}$$

> $$\frac{|S_V(p,M)|}{|S_U(p,M)|} \;=\; \mathbb{E}_U\!\left[\frac{n!}{|Γ|}\right]$$

The first is the important one. Since `|Γ_F| >= 1`, it lies in `(0, 1]`, and **it equals 1 exactly
when every class is facet-asymmetric — which is exactly the condition for `S_F` and `S_U` to
coincide.** It is computable from two counter calls. No sampling, no estimator, no error bar.

The Radon–Nikodym derivative is

$$\frac{dF}{dU}(C) \;=\; \frac{1/|Γ_F(C)|}{\mathbb{E}_U[1/|Γ_F|]},$$

so the total variation distance between the two ensembles, over isomorphism classes, is

$$d_{TV}(F, U) \;=\; \tfrac12\, \mathbb{E}_U\!\left[\left|\frac{1}{|Γ_F|\;\mathbb{E}_U[1/|Γ_F|]} - 1\right|\right].$$

Because the normalizer is known **exactly**, this is an importance-sampling estimate against a known
constant — far better conditioned than any two-sample test, and it measures the whole distribution
rather than one marginal. **This is the quantity P1 should lead with.**

## 6. Consequences

**Facet-labeled → unlabeled.** The weight ratio is bounded and tends to 1 as facet-moving
automorphisms die out, so the two ensembles converge, and every statistic computed on them converges
with them. The rate is `E_U[1/|Γ_F|]`.

**Vertex-labeled stays apart.** The weight `n!/|Γ|` grows factorially in a vertex count that varies
across classes, so `S_V` is pushed hard toward the top of the vertex-count range and cannot converge
to `S_U`. **The vertex count is therefore the sharpest discriminator, and also the cheapest
statistic to measure** — a good first figure. Observed at `p = 3, M = 5` (maximum possible
`n = pM = 15`), in 8 draws from each:

| ensemble | vertex counts drawn |
|---|---|
| `S_V` | 15, 12, 12, 13, 13, 13, 13, 14 |
| `S_F` | 11, 10, 6, 10, 6, 9, 9, 7 |
| `S_U` | 8, 8, 8, 8, 7, 8, 6, 9 |

## 7. Measured

`E_U[1/|Γ_F|] = |S_F| / (M! |S_U|)`, exact:

| p (d = p−1) | M=2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|
| 2 (graphs) | .500 | .300 | .265 | .254 | **.246** | .261 | .277 | .297 |
| 3 | .500 | **.403** | .420 | .531 | .628 | .721 | .788 | .835 |
| 4 | .500 | **.453** | .546 | .707 | .830 | .903 | .940 | .960 |

Three things to read off:

1. **It is not monotone.** It starts at 1 (`M = 1`, trivially one class), falls as symmetric
   configurations proliferate, bottoms out (bold), then climbs toward 1 as complexes become
   generically asymmetric.
2. **The rate depends sharply on purity.** By `M = 9` it is 0.96 at `p = 4`, 0.84 at `p = 3`, and
   only 0.30 at `p = 2`. **Convergence is dramatically slower for graphs than for higher-dimensional
   complexes**, and the `p = 2` row has barely turned the corner by `M = 9`.
3. **`M = 1` and `M = 2` are exactly 1 and 1/2, for every purity.** Provable by hand — see §8.

And `E_U[n!/|Γ|] = |S_V|/|S_U|`, showing the scale of the vertex-labeled distortion:

| p | M=2 | M=4 | M=6 | M=8 | for scale: `log10((pM)!)` at M=8 |
|---|---|---|---|---|---|
| 2 | 3.00 | 81.8 | 5.66e3 | 6.53e5 | 13.3 |
| 3 | 10.3 | 6.68e3 | 1.16e7 | 2.55e10 | 23.8 |
| 4 | 40.0 | 8.32e5 | 3.74e10 | 1.89e15 | 35.4 |

The last column matters: `E_U[n!/|Γ|]` sits far below `(pM)!`, so the mass is **not** concentrated
at exactly `n = pM`. Most classes have far fewer vertices, and the factorial weight pushes `S_V`
toward the high end of the achievable range rather than onto its endpoint — consistent with the
12–15 spread observed above.

## 8. Verification

- **Sampling check of the identity.** 2000 draws from `S_U` at `p = 3, M = 5`, averaging
  `1/|Γ_F|` via `PureComplexFacetAutomorphismGroupOrder`: **0.5236 sampled** against
  **0.5308 exact**. Agreement to ~1.4%, consistent with the sampling error on 2000 draws.
- **`M = 1`.** One facet, so `Γ_F <= S_1` is trivial and the ratio is exactly 1. Confirmed for
  `p = 2..5` in the test suite (`ensemble-ratio-M1-is-one`).
- **`M = 2`.** For any two distinct `p`-subsets `A ≠ B`, fix `A ∩ B` pointwise and swap `A \ B` with
  `B \ A` — they have equal size `p − |A ∩ B|`, so such a permutation always exists and swaps the
  two facets. Hence `|Γ_F| = 2` for *every* two-facet complex and the ratio is exactly 1/2.
  Confirmed for `p = 2..5` (`ensemble-ratio-M2-is-one-half`).
- The automorphism functions themselves are pinned against known groups: `S_4` on the tetrahedron
  boundary (24), the octahedral group (48), `S_3` on a single triangle (6), 72 on two disjoint
  triangles.

## 9. What P1 should compute

Statistics, all with `n` free: **number of vertices**, **number of connected components**,
**Betti numbers**, **automorphism group order**.

- [ ] The exact `E_U[1/|Γ_F|]` curve per purity, extended in `M` until each crosses ~0.99. This is
      the convergence result and it costs two counter calls per point.
- [ ] `d_TV(F, U)` by the §5 formula — sampling from `S_U` only, with the exact normalizer.
- [ ] Per-statistic empirical distributions for all three ensembles, and their pairwise distances.
- [ ] The vertex-count distribution as the `S_V` contrast (§6).

### Statistical methodology — one trap

**Every statistic here is discrete**, and `|Aut|` takes divisor values, so it is violently lumpy.
The classical two-sample Kolmogorov–Smirnov test is **conservative under ties**: p-values are
inflated and the test under-rejects. Since the conclusion being argued is *"F and U agree"* — that
is, a *failure to reject* — **a conservative test is biased toward the desired answer.** This will
not survive refereeing.

- Use **permutation-based p-values** for the KS (or Cramér–von Mises, or Anderson–Darling)
  statistic. They are exact under the null regardless of discreteness, and cost about ten lines.
- **Report the effect size, not the p-value.** Failure to reject at each `M` says little, since
  power is a function of sample size and can be manufactured. The evidence for convergence is
  `D_M → 0`: plot the statistic against `M` with bootstrap intervals, and overlay the exact
  `E_U[1/|Γ_F|]` curve as the theoretical prediction. Agreement between an exactly computed rate
  and an empirical one is a far stronger result than any hypothesis test.
- For `|Aut|`, compare **`log |Aut|`**: group orders are multiplicative and heavy-tailed.

## 10. Open

- Does `|Γ_F| → 1` in probability under `S_U` admit a clean proof for pure complexes at fixed `M`,
  in the spirit of the classical result that random graphs are asymmetric? The measured table says
  it happens; the rate's purity dependence is unexplained.
- Why is the dip at `M = 6` for `p = 2` but `M = 3` for `p = 3, 4`?
- The `p = 2` row is the slowest and is also the case that connects to *Dynamical Quantum
  Multigraphs*, where the labeled/unlabeled distinction produced qualitatively different
  thermodynamics. Whether the slow convergence here is the same phenomenon is worth a paragraph.

---

## Related

- `FacetLabeledCount.md`, `UnlabeledCount.md` — the counting formulas behind `|S_F|` and `|S_U|`.
- `FacetLabeledSampler.md`, `UnlabeledSampler.md` — the exact samplers.
- `ExpansionRoadmap.md` Track A — where this sits in P1.
