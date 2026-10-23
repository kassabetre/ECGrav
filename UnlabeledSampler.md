# Sampling fully unlabeled pure complexes, exactly

> **What this is.** A specification of `RandomUniformUnlabeledPureSimplicialComplex`, matching
> the implementation in `Kernel/Subpackages/PureComplexes.wl`. Companion to
> `UnlabeledCount.md`, which specifies the counter this sits on. It records the sampling
> principle, the five steps, the lemma that had to be corrected, why covering is enforced during
> the draw rather than by rejection, and how uniformity is established — by exact identities
> rather than by goodness-of-fit.

Notation as in `UnlabeledCount.md`: $`p`$ purity, $`M`$ facet order, $`n`$ vertex count,
$`\lambda \vdash n`$ an integer partition used as a cycle type, $`m_k`$ the number of parts of
$`\lambda`$ equal to $`k`$, and $`z_\lambda = \prod_k k^{m_k} m_k!`$, so $`n!/z_\lambda`$ permutations
have that type. $`U(p,M,n)`$ is `NumUnlabeledPureComplexes`, $`\mathrm{VL}(p,M,n)`$ is
`NumVertexLabeledPureComplexes`, $`A(p,M,n)`$ is the count with covering not required
(`UnlabeledCount.md` §2.2), and $`n_d`$ is the number of $`\langle\sigma\rangle`$-orbits of size $`d`$
on $`p`$-subsets (§2.4 there).

---

## 1. What is being sampled, and what "uniform" means

The objects are the **isomorphism classes** of pure complexes: sets $`S`$ of $`M`$ pairwise distinct
$`p`$-subsets of $`[n]`$ covering $`[n]`$, up to relabelling the vertices — the $`S_n`$-orbits, counted
by $`U(p,M,n)`$. Uniform means each of the $`U(p,M,n)`$ classes has probability exactly $`1/U(p,M,n)`$,
**not** that each labelled complex is equally likely. Those differ sharply: at $`(p,M,n)=(2,3,4)`$
there are 16 labelled complexes but two classes — the path $`P_4`$ with 12 labellings and the star
$`K_{1,3}`$ with 4 — so a uniform labelled draw splits them $`75/25`$.

Output is a facet list on vertex set $`[n]`$; it is a representative of its class, and the
representative's labelling is itself uniform.

---

## 2. Why the two obvious approaches fail

**Rejection from the vertex-labeled generator.** A class with automorphism group $`\mathrm{Aut}`$
has $`n!/|\mathrm{Aut}|`$ labellings, so accepting a uniform labelled draw with probability
$`|\mathrm{Aut}|/n!`$ cancels the over-representation exactly. Correct, and this is what the
package did before. But the acceptance rate collapses factorially: measured against that
implementation, **146 ms** per sample at $`\{3,4,6\}`$, **3.31 s** at $`\{3,4,8\}`$, **11.4 s** at
$`\{3,4,10\}`$.

**Porting the facet-labeled sampler.** That one draws a cycle type, cuts a uniform permutation
into blocks, and picks facets sequentially. Its state, and the exactness lemma licensing it, do
not survive here — see §7. The failure is silent: the naive port produces perfectly plausible
complexes from a subtly wrong distribution.

---

## 3. The principle: pair sampling

Let

```math
\Phi \;=\; \{(\sigma, S) \;:\; \sigma \in S_n,\ S \text{ a } \sigma\text{-invariant covering } M\text{-set}\}.
```

For any $`S_n`$-orbit $`O`$ (i.e. any isomorphism class),

```math
\#\{(\sigma,S) \in \Phi : S \in O\} \;=\; \sum_{S \in O} |\mathrm{Stab}(S)| \;=\; |O| \cdot \frac{n!}{|O|} \;=\; n!,
```

**the same number for every orbit, whatever its size.** So drawing $`(\sigma,S)`$ uniformly from
$`\Phi`$ and discarding $`\sigma`$ leaves the class exactly uniform. No automorphism weighting, no
rejection. ($`|\Phi| = \sum_\sigma \lvert\mathrm{Fix_{cov}}(\sigma)\rvert = n!\,U(p,M,n)`$, which is the identity
§10 checks.)

This is the same argument that made `RandomUniformFacetLabeledPureSimplicialComplex` exact; only
the inner object changes, from an ordered tuple to a $`\sigma`$-invariant set.

---

## 4. The five steps

```
1.  λ  ~  ∝ (n!/z_λ) · |Fix_cov(λ)|                      cycle type of σ
2.  work on the canonical σ_λ (cycles = consecutive blocks)
3.  a  ~  ∝ C_cov(λ, a)   over profiles with Σ_d d·a_d = M       size profile
4.  choose a_d distinct ⟨σ_λ⟩-orbits of size d, for each d, covering [n]
5.  relabel by a uniform random π ∈ S_n
```

A $`\sigma`$-invariant set is a union of $`\langle\sigma\rangle`$-**orbits** of $`p`$-subsets, which is
why steps 3 and 4 are phrased in orbits rather than facets.

---

## 5. Step 1 — the cycle type

$`\lvert\mathrm{Fix_{cov}}(\lambda)\rvert`$ is the number of $`\sigma`$-invariant **covering** $`M`$-sets for
$`\sigma`$ of type $`\lambda`$. The counter never forms this — it reaches covering by differencing
$`U(n) = A(n) - A(n-1)`$ over the vertex count, which says nothing per cycle type — so the sampler
supplies it, by inclusion–exclusion over **deleted cycles**:

```math
\lvert\mathrm{Fix_{cov}}(\lambda)\rvert \;=\; \sum_{\rho \subseteq \lambda} (-1)^{\ell(\lambda)-\ell(\rho)} \prod_k \binom{m_k(\lambda)}{m_k(\rho)} \cdot \lvert\mathrm{Fix}(\rho)\rvert,
```

$`\rho`$ ranging over sub-multisets, $`\ell`$ the number of parts, and $`\lvert\mathrm{Fix}\rvert`$ the counter's own
non-covering `NumULPCFixSets`. The identity behind it: the sets avoiding a cycle set $`T`$ are
exactly the sets of the permutation with those cycles removed, so they are counted at the reduced
cycle type. This one idea supplies every count in the sampler.

The weight for $`\lambda`$ is $`(n!/z_\lambda)\cdot\lvert\mathrm{Fix_{cov}}(\lambda)\rvert`$. Cost of the whole
table is about $`3\times`$ a single call to the counter.

---

## 6. Step 2 and 5 — canonical $`\sigma`$, relabel last

Instead of realising a uniform $`\sigma`$ of type $`\lambda`$, the draw works on the **canonical**
$`\sigma_\lambda`$ whose cycles are consecutive blocks $`\{1..k_1\}, \{k_1{+}1..k_1{+}k_2\}, \dots`$,
and applies a uniform $`\pi \in S_n`$ at the very end.

This is legitimate because conjugation by $`\pi`$ is a bijection from $`\sigma_\lambda`$-invariant
covering sets to $`\pi\sigma_\lambda\pi^{-1}`$-invariant ones, and $`\pi\sigma_\lambda\pi^{-1}`$
ranges uniformly over the permutations of type $`\lambda`$. The pair $`(\sigma,S)`$ therefore comes
out uniform over its slice of $`\Phi`$.

It matters for cost, not elegance: **every orbit computation becomes a function of $`\lambda`$
alone**, so the orbit list, orbit-size counts, profiles and completion counts are memoised across
samples instead of rebuilt per sample.

---

## 7. Step 3 — the size profile, and why it is drawn up front

### 7.1 The exactness lemma does not carry over

The facet-labeled sampler's sequential draw rests on: *completions depend on the history only
through (how many objects are used, which blocks are uncovered), never on which were used.* **That
is false here.** It held in only 11 of 16 tested cycle types.

Minimal counterexample. Take $`\sigma = (1\,2)(3\,4)`$, $`p = 2`$, $`M = 3`$. The
$`\langle\sigma\rangle`$-orbits on the 2-subsets of $`[4]`$ are

| orbit | size | support | members |
| --- | --- | --- | --- |
| #1 | 1 | cycle 1 | $`\{1,2\}`$ |
| #2 | 2 | cycles 1,2 | $`\{1,3\},\{2,4\}`$ |
| #3 | 2 | cycles 1,2 | $`\{1,4\},\{2,3\}`$ |
| #4 | 1 | cycle 2 | $`\{3,4\}`$ |

Three partial selections share the state *(2 units used, nothing uncovered)*:

| selection | size profile | completions |
| --- | --- | --- |
| $`\{\#2\}`$ | one of size 2 | **2** |
| $`\{\#3\}`$ | one of size 2 | **2** |
| $`\{\#1,\#4\}`$ | two of size 1 | **0** |

$`\{\#1,\#4\}`$ spent its budget in units of 1 and needs one more unit, but every remaining orbit
has size 2 — it has no completions at all. Same "amount used", different futures.

The reason the facet-labeled case escapes this: there every chosen object costs exactly one unit
of the budget, so the profile *is* the total. Here orbit sizes exceed 1 and $`M`$ is a sum of sizes,
so two histories can spend the same budget in different denominations.

The correct state is **(size profile $`\alpha_d`$, uncovered cycles)**, verified 16/16, and with the
uncovered set compressed to its **length multiset**, still 16/16.

### 7.2 Sampling the profile dissolves the problem

Rather than carry $`\alpha`$ through the draw, sample the whole profile first. That this is even
available comes from reading the count the other way round.

The isomorphism classes are also the $`S_M`$-orbits of **separating incidence tableaux** — transpose
the incidence matrix and a facet-labeled complex becomes a multiset of $`n`$ nonempty subsets of
$`[M]`$, one row per vertex, every label in exactly $`p`$ rows, columns pairwise distinct. Acting on
those by $`\tau \in S_M`$ and averaging gives the same count. And a $`\tau`$-cycle of length $`d`$ pins
exactly one $`\sigma`$-orbit of size $`d`$, so:

> $`\tau`$'s cycle type **is** the size profile.

The profile is therefore a first-class quantity with its own marginal, not hidden state. That
marginal is §5's inclusion–exclusion refined by profile:

```math
C_{\mathrm{cov}}(\lambda, a) \;=\; \sum_{T \subseteq \text{cycles}} (-1)^{|T|} \prod_d \binom{n_d(T)}{a_d},
```

with $`n_d(T)`$ the orbits of size $`d`$ avoiding $`T`$. Summing over profiles returns
$`\lvert\mathrm{Fix_{cov}}(\lambda)\rvert`$, which §10 checks.

*(As a counting algorithm the $`S_M`$ route is not better — the outer loop shrinks from
$`\mathrm{PartitionsP}(n)`$ to $`\mathrm{PartitionsP}(M)`$, but the inner sum over $`\sigma`$ no longer
collapses well. Its value is entirely this structural one.)*

---

## 8. Step 4 — drawing the orbits

With the profile fixed, the task is: pick $`a_d`$ distinct orbits of size $`d`$ for each $`d`$, such
that their supports cover every cycle.

**Completion counts.** For a partial selection with $`\alpha_d`$ orbits of size $`d`$ used and $`\mathcal{U}`$
still uncovered,

```math
K(\alpha, \mathcal{U}) \;=\; \sum_{T \subseteq \mathcal{U}} (-1)^{|T|} \prod_d \binom{n_d(T) - \alpha_d}{\,a_d - \alpha_d\,}.
```

The subtraction is exactly right because an orbit already chosen has its support inside the
*covered* cycles, so it avoids every $`T \subseteq \mathcal{U}`$ and survives into every term — only its count
per size enters. Verified against brute force on 333 partial selections, 288 of them nonzero,
zero mismatches.

**The draw.** Orbits are taken one at a time, each weighted by $`K`$ of the state it leads to. This
is exact by the chain rule, and it produces an *ordered sequence*; the order is then forgotten.
Every valid set has the same number of orderings, so the underlying set is uniform, and the
ordering factor cancels out of the weights. The invariant that makes this work is

```math
\sum_{O \text{ candidate}} K(\text{state after } O) \;=\; t \cdot K(\text{state}),
```

$`t`$ the number still to pick — each completing set counted once per choice of its next element.
Asserted at **every step of every draw** in testing: 0 violations in ~4,600 steps.

### A worked draw at $`\{p,M,n\} = \{2,4,5\}`$

Four edges covering five vertices. It is small enough to follow by hand, but every step does real
work: several cycle types carry weight, the inclusion–exclusion for covering has nonzero
correction terms, one cycle type has two admissible profiles, candidate orbits share buckets, and
one bucket is shut out by covering. Every number below comes from the shipped helpers
(`RandULPCTypeWeights`, `RandULPCOrbits`, `RandULPCProfiles`, `RandULPCCCov`, `RandULPCComps`), and the
draw is a seeded run of the same logic as `RandULPCDraw`.

**The target.** $`U(2,4,5) = 4`$, and the four classes have very different automorphism groups:

| class | $`\lvert\mathrm{Aut}\rvert`$ | labelled complexes | uniform labelled draw | target |
| --- | --- | --- | --- | --- |
| path $`P_5`$ | 2 | 60 | $`4/9`$ | $`1/4`$ |
| chair (a vertex of degree 3 with one leg extended) | 2 | 60 | $`4/9`$ | $`1/4`$ |
| triangle plus an edge, $`K_3 \cup K_2`$ | 12 | 10 | $`2/27`$ | $`1/4`$ |
| star $`K_{1,4}`$ | 24 | 5 | $`1/27`$ | $`1/4`$ |

The labelled complexes total $`\mathrm{VL}(2,4,5) = 135`$, which is why a uniform labelled draw is so
far from the target.

**Step 1: the cycle type.**

| $`\lambda`$ | $`n!/z_\lambda`$ | $`\lvert\mathrm{Fix_{cov}}(\lambda)\rvert`$ | weight |
| --- | --- | --- | --- |
| $`(5)`$ | 24 | 0 | 0 |
| $`(4,1)`$ | 30 | 1 | 30 |
| $`(3,2)`$ | 20 | 1 | 20 |
| $`(3,1,1)`$ | 20 | 3 | 60 |
| $`(2,2,1)`$ | 15 | 7 | 105 |
| $`(2,1,1,1)`$ | 10 | 13 | 130 |
| $`(1^5)`$ | 1 | 135 | 135 |

The weights sum to $`480 = 5! \cdot 4`$, which is I1. For $`(2,1,1,1)`$ the inclusion–exclusion of §5
has one nonzero correction. There are $`\lvert\mathrm{Fix}(2,1,1,1)\rvert = 22`$ invariant 4-sets, of
which $`3\,\lvert\mathrm{Fix}(2,1,1)\rvert = 3 \cdot 3 = 9`$ miss one of the three fixed points, so
$`\lvert\mathrm{Fix_{cov}}\rvert = 22 - 9 = 13`$. Every other term vanishes, because four edges cannot
fit on three vertices. The draw takes $`(2,1,1,1)`$ with probability $`130/480 = 13/48`$.

**Step 2: the canonical $`\sigma`$.** For $`(2,1,1,1)`$ this is $`\sigma_\lambda = (1\,2)`$, with cycles
$`c_1 = \{1,2\}`$, $`c_2 = \{3\}`$, $`c_3 = \{4\}`$, $`c_4 = \{5\}`$. Its $`\langle\sigma\rangle`$-orbits on
the ten edges are

| orbit | size | support | members |
| --- | --- | --- | --- |
| $`O_1`$ | 1 | $`c_1`$ | $`\{1,2\}`$ |
| $`O_2`$ | 2 | $`c_1, c_2`$ | $`\{1,3\},\{2,3\}`$ |
| $`O_3`$ | 2 | $`c_1, c_3`$ | $`\{1,4\},\{2,4\}`$ |
| $`O_4`$ | 2 | $`c_1, c_4`$ | $`\{1,5\},\{2,5\}`$ |
| $`O_5`$ | 1 | $`c_2, c_3`$ | $`\{3,4\}`$ |
| $`O_6`$ | 1 | $`c_2, c_4`$ | $`\{3,5\}`$ |
| $`O_7`$ | 1 | $`c_3, c_4`$ | $`\{4,5\}`$ |

so $`n_1 = 4`$ and $`n_2 = 3`$.

**Step 3: the size profile.** A profile needs $`a_1 + 2a_2 = 4`$:

| profile $`(a_1, a_2)`$ | full term | minus one fixed point deleted | $`C_{\mathrm{cov}}(\lambda,a)`$ | probability |
| --- | --- | --- | --- | --- |
| $`(4,0)`$ | $`\binom44 = 1`$ | $`3\binom24 = 0`$ | 1 | $`1/13`$ |
| $`(2,1)`$ | $`\binom42\binom31 = 18`$ | $`3\binom22\binom21 = 6`$ | 12 | $`12/13`$ |
| $`(0,2)`$ | $`\binom32 = 3`$ | $`3\binom22 = 3`$ | 0 | 0 |

The deleted-point counts use the orbit counts of the reduced type $`(2,1,1)`$, which has
$`n_1 = 2`$ and $`n_2 = 2`$. Profile $`(0,2)`$ is admissible in size but cannot cover: two pair orbits
reach only two of the three fixed points. The profile column sums to $`13 = \lvert\mathrm{Fix_{cov}}\rvert`$,
which is I2. The draw takes $`(2,1)`$: two orbits of size 1 and one of size 2.

**Step 4: the orbits.** Candidates are bucketed by (size, support $`\cap\,\mathcal{U}`$). A bucket's
weight is its size times $`K`$ of the state it leads to, and within a bucket the orbit is uniform.

*First orbit.* Nothing is covered, so $`\mathcal{U} = \{c_1,c_2,c_3,c_4\}`$, $`K = 12`$ and $`t = 3`$.
Every orbit is its own bucket:

| bucket | orbits | $`K`$ after | weight | probability |
| --- | --- | --- | --- | --- |
| $`(1, \{c_1\})`$ | $`O_1`$ | 3 | 3 | $`1/12`$ |
| $`(2, \{c_1,c_2\})`$, $`(2, \{c_1,c_3\})`$, $`(2, \{c_1,c_4\})`$ | $`O_2`$, $`O_3`$, $`O_4`$ | 4 each | 4 each | $`1/9`$ each |
| $`(1, \{c_2,c_3\})`$, $`(1, \{c_2,c_4\})`$, $`(1, \{c_3,c_4\})`$ | $`O_5`$, $`O_6`$, $`O_7`$ | 7 each | 7 each | $`7/36`$ each |

The weights sum to $`36 = t\,K`$, which is I3. The draw takes $`O_6 = \{3,5\}`$, covering $`c_2`$ and $`c_4`$.

*Second orbit.* Now $`\mathcal{U} = \{c_1, c_3\}`$, $`K = 7`$ and $`t = 2`$, and buckets start to merge:

| bucket | orbits | $`K`$ after | weight | probability |
| --- | --- | --- | --- | --- |
| $`(1, \{c_1\})`$ | $`O_1`$ | 1 | 1 | $`1/14`$ |
| $`(2, \{c_1\})`$ | $`O_2`$, $`O_4`$ | 2 | 4 | $`2/7`$ |
| $`(2, \{c_1,c_3\})`$ | $`O_3`$ | 3 | 3 | $`3/14`$ |
| $`(1, \{c_3\})`$ | $`O_5`$, $`O_7`$ | 3 | 6 | $`3/7`$ |

The weights sum to $`14 = t\,K`$. The draw takes the bucket $`(2,\{c_1\})`$ and, uniformly within it,
$`O_2 = \{\{1,3\},\{2,3\}\}`$. This covers $`c_1`$ and uses up the one orbit of size 2.

*Third orbit.* Now $`\mathcal{U} = \{c_3\}`$, $`K = 2`$ and $`t = 1`$:

| bucket | orbits | $`K`$ after | weight | probability |
| --- | --- | --- | --- | --- |
| $`(1, \varnothing)`$ | $`O_1`$ | 0 | 0 | 0 |
| $`(1, \{c_3\})`$ | $`O_5`$, $`O_7`$ | 1 | 2 | 1 |

$`O_1`$ is still available by size but would leave vertex 4 uncovered, so its completion count is 0.
This is covering being enforced during the draw rather than by rejection (§9). The draw takes
$`O_5 = \{3,4\}`$.

This ordered path had probability $`\tfrac{7}{36} \cdot \tfrac27 \cdot \tfrac12 \cdot \tfrac12 = \tfrac1{72}`$.
The same three orbits arrive in $`3! = 6`$ orders, each with probability $`1/72`$, so the set has
probability $`1/12 = 1/C_{\mathrm{cov}}(\lambda,a)`$, as every set with this profile does.

**Step 5: relabel.** The set is $`S = \{\{3,5\},\{1,3\},\{2,3\},\{3,4\}\}`$, the star $`K_{1,4}`$ centred at
vertex 3, which $`\sigma = (1\,2)`$ visibly fixes. The uniform relabelling drawn,
$`\pi = (1\,2)`$, happens to be an automorphism of $`S`$, so the output is the same facet list.

**Why every class comes out at exactly $`1/4`$.** For each class, count the $`\sigma_\lambda`$-fixed
covering sets of that class, for each cycle type, and weight by $`n!/z_\lambda`$:

| class | $`(4,1)`$ ×30 | $`(3,2)`$ ×20 | $`(3,1,1)`$ ×20 | $`(2,2,1)`$ ×15 | $`(2,1,1,1)`$ ×10 | $`(1^5)`$ ×1 | weighted total |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $`P_5`$ | 0 | 0 | 0 | 4 | 0 | 60 | 120 |
| chair | 0 | 0 | 0 | 0 | 6 | 60 | 120 |
| $`K_3 \cup K_2`$ | 0 | 1 | 1 | 2 | 4 | 10 | 120 |
| $`K_{1,4}`$ | 1 | 0 | 2 | 1 | 3 | 5 | 120 |
| column total | 1 | 1 | 3 | 7 | 13 | 135 | 480 |

Every row is $`\sum_{S \in c} \lvert\mathrm{Stab}(S)\rvert = n! = 120`$, the orbit–stabiliser count of §3, so
every class gets $`120/480 = 1/4`$; the column totals are the $`\lvert\mathrm{Fix_{cov}}(\lambda)\rvert`$ of
step 1. The identity column alone is the uniform labelled draw. The star collects only 5 of its
120 there; the other 115 come from the nontrivial cycle types, in proportion to its 24
automorphisms. Following every branch of the shipped draw exactly, rather than counting fixed
sets, gives the same $`1/4`$ for each class, with I3 holding at every step.

---

## 9. Covering during the draw, not after

Enforcing covering by rejecting non-covering draws is also exact and much simpler. Its expected
trial count has a closed form. Writing $`q(\lambda,a)`$ for the probability of drawing that pair and
$`C_{\mathrm{all}}(\lambda,a) = \prod_d \binom{n_d}{a_d}`$ for the unconstrained count, the
acceptance is $`C_{\mathrm{cov}}/C_{\mathrm{all}}`$ and

```math
\mathbb{E}[\text{trials}] \;=\; \sum_{\lambda,a} q(\lambda,a)\,\frac{C_{\mathrm{all}}}{C_{\mathrm{cov}}} \;=\; \frac{1}{n!\,U(p,M,n)}\sum_\lambda \frac{n!}{z_\lambda}\,\lvert\mathrm{Fix}(\lambda)\rvert \;=\; \frac{A(p,M,n)}{U(p,M,n)}.
```

Two things follow. First, $`C_{\mathrm{cov}}`$ **cancels**: weighting $`\lambda`$ and $`a`$ toward
covering buys nothing, since it raises $`q`$ exactly where acceptance is already high. Second, the
cost is the ratio of the non-covering to the covering count, and that is ruinous at the top of the
vertex range, where the facets must nearly partition $`[n]`$:

| $`\{p,M\}`$ | median over $`n`$ | worst | at $`n`$ |
| --- | --- | --- | --- |
| $`\{2,4\}`$ | 3.0 | 11 | 8 |
| $`\{3,4\}`$ | 3.7 | 66 | 12 |
| $`\{3,6\}`$ | 12.4 | 4,279 | 18 |
| $`\{4,7\}`$ | 56.3 | **11,662,671** | 28 |

The worst case is always exactly $`n = pM`$. Hence covering is built into the completion counts.

---

## 10. Why it is uniform

Uniformity is not inferred from goodness-of-fit. Three identities force it.

| | identity | what it pins |
| --- | --- | --- |
| **I1** | $`\sum_\lambda (n!/z_\lambda)\,\lvert\mathrm{Fix_{cov}}(\lambda)\rvert = n!\,U(p,M,n)`$ | step 1 marginal |
| **I2** | $`\sum_a C_{\mathrm{cov}}(\lambda,a) = \lvert\mathrm{Fix_{cov}}(\lambda)\rvert`$ | step 3 marginal |
| **I3** | $`\sum_O K(\text{after } O) = t\,K(\text{state})`$ | step 4 draw is uniform |

Given I3 the conditional draw is uniform, so $`\Pr[S \mid \lambda,a] = 1/C_{\mathrm{cov}}(\lambda,a)`$,
and the class probability telescopes:

```math
\Pr[c] \;=\; \sum_\lambda \frac{(n!/z_\lambda)\lvert\mathrm{Fix_{cov}}(\lambda)\rvert}{n!\,U(p,M,n)} \sum_a \frac{C_{\mathrm{cov}}(\lambda,a)}{\lvert\mathrm{Fix_{cov}}(\lambda)\rvert} \cdot \frac{\#\{S \in c \text{ with profile } a\}}{C_{\mathrm{cov}}(\lambda,a)} \;=\; \frac{1}{n!\,U(p,M,n)}\sum_\sigma \#\{S \in c : \sigma S = S\},
```

and the inner sum is $`\sum_{S \in c} |\mathrm{Stab}(S)| = n!`$. So $`\Pr[c] = 1/U(p,M,n)`$ **for every class**.

That was also computed directly rather than sampled — enumerating the classes, summing the
per-cycle-type fixed counts, and comparing — giving exactly $`1/U(p,M,n)`$ for every class on seven
parameter sets, including the 15 classes of $`\{3,4,6\}`$. A third identity worth noting:
$`\lvert\mathrm{Fix_{cov}}(1^n)\rvert = \mathrm{VL}(p,M,n)`$, since the identity
permutation's invariant covering sets *are* the vertex-labeled complexes — it pins the covering
count against a function computed a wholly different way.

---

## 11. Implementation map

All in `Kernel/Subpackages/PureComplexes.wl`; the derivation is the header comment at line 3682.
Private helpers are named `RandULPC*`.

| line | symbol | role |
| --- | --- | --- |
| 3735, 3738 | `$RandomUnlabeledParallelThreshold`, `RandULPCGoParallel` | sample-count gate, also gated on `$KernelCount` |
| 3742 | `RandULPCOrbits[λ,p]` | the orbits as `{size, support, members}` — the only place the $`\binom{n}{p}`$ subsets are touched |
| 3762 | `RandULPCNd[λ,p]` | $`n_d`$ by Möbius inversion off the counter's $`f(e)`$ (`UnlabeledCount.md` §2.5), so reduced cycle types cost nothing |
| 3779 | `RandULPCSubMulti[λ]` | sub-multisets with inclusion–exclusion signs and binomial weights |
| 3792 | `RandULPCFixCov[λ,p,M]` | $`\mathrm{Fix_{cov}}`$ (§5) |
| 3799 | `RandULPCTypeWeights[p,M,n]` | step-1 weights |
| 3810 | `RandULPCProfiles[λ,p,M]` | admissible size profiles |
| 3824 | `RandULPCCCov[λ,p,M,a]` | step-3 weights (§7.2) |
| 3834 | `RandULPCComps[...]` | completion counts (§8) |
| 3855 | `RandULPCDraw[λ,p,M,a]` | step 4 |
| 3900 | `RandULPCOne[p,M,n]` | one sample, steps 1–5 |
| 3919 | primary pattern | `[{p,M,n}, numSamples]` |
| 3943 | overload | `[{p,M}, numSamples]` |
| 3979, 3988 | overloads | single-sample forms |

All memo tables are released by the shared ``ECGrav`Private`NumPCClearCache[]``.

With the vertex count free, $`n`$ is drawn from `NumUnlabeledPureComplexes` directly. The previous
implementation drew it from the **vertex-labeled** counts and let rejection repair the difference,
which is not the distribution that was wanted.

---

## 12. Optimisations

Three, each exact rather than heuristic:

1. **Completion counts indexed by the uncovered *length multiset***, not the uncovered cycle set.
   Sound because $`n_d(T)`$ depends only on the shape of the surviving cycles, so states of equal
   shape share a memo entry. Checked separately: at size boundaries the compression holds 35/35;
   mid-class it fails 47/59, which is why per-size decisions are atomic.
2. **Inclusion–exclusion over sub-multisets** with binomial weights: writing $`u_k`$ for the number
   of uncovered cycles of length $`k`$, a sub-multiset is a choice of how many of each length to
   take, so there are $`\prod_k(u_k{+}1)`$ terms instead of $`2^{\sum_k u_k}`$ subsets.
3. **Candidates bucketed by $`(\text{size}, \mathrm{supp} \cap \mathcal{U})`$** — orbits sharing that pair lead
   to the same state, so one completion count serves a whole bucket.

Together these took $`\{3,5,12\}`$ from 455 ms to 2.1 ms per sample.

---

## 13. Cost

Measured warm, against the previous rejection implementation extracted from the prior commit:

| $`\{p,M,n\}`$ | before | after | speedup |
| --- | --- | --- | --- |
| $`\{2,3,4\}`$ | 3.7 ms | 0.12 ms | 31× |
| $`\{3,4,6\}`$ | 146 ms | 0.32 ms | 458× |
| $`\{2,6,6\}`$ | 188 ms | 0.38 ms | 496× |
| $`\{3,4,8\}`$ | 3.31 s | 0.68 ms | 4,800× |
| $`\{3,4,10\}`$ | 11.4 s | 0.82 ms | **14,000×** |
| $`\{2,6,10\}`$ | — | 0.81 ms | |
| $`\{3,5,12\}`$ | — | 2.14 ms | |
| $`\{3,6,14\}`$ | — | 4.58 ms | |

The first call at a given $`(p,M,n)`$ builds the memo tables; the figures above are warm, which is
the regime any real use is in.

**`$RandomUnlabeledParallelThreshold = 500`**, measured on 11 subkernels: serial/parallel runs
0.24–0.53 at 100 samples, 0.31–0.81 at 200, 0.78–1.48 at 400, 1.47–2.71 at 800, 2.7–4.2 at 3200.
It sits far above the facet-labeled generator's 100 because nearly all the work here is memo
tables keyed on the cycle type and **every subkernel rebuilds its own** — parallel is genuinely
slower below a few hundred samples.

---

## 14. Tests

In `Tests/PureComplexes.wlt`:

| TestID | what it pins |
| --- | --- |
| `-weight-identity` | I1, over 120 parameter sets |
| `-profile-marginals` | I2, every cycle type of every set |
| `-fixcov-at-identity` | $`\lvert\mathrm{Fix_{cov}}(1^n)\rvert`$ = vertex-labeled count |
| `-valid-samples` | $`M`$ distinct $`p`$-subsets covering $`[n]`$, across 9 sets including $`n = pM`$ |
| `-uniform-15-classes` | chi-square over the 15 classes of $`\{3,4,6\}`$ |
| `-vertex-count-marginal` | `{p,M}` draws $`n`$ from `NumUnlabeledPureComplexes` |
| `-uniform-over-classes` | the pre-existing $`P_4`$-vs-$`K_{1,3}`$ check |

I1–I3 are the ones that matter. A wrong cycle-type weight still emits plausible complexes, so
identity checks catch what output inspection cannot.

---

## 15. Notes

- **Chi-square is the weak test here.** Repeated runs at $`\{3,4,6\}`$ produced p-values as low as
  0.048 purely by chance, which cost real time before the exact route was taken. Prefer I1–I3 and
  the direct class-probability computation.
- **A doc page's example output is algorithm-dependent.** The reference page's `SeedRandom[1]`
  example had to be regenerated: the prose was still true, but the value was produced by the old
  sampler.
- **Guard against checks that cannot fail.** Two verification scripts silently passed everything
  before being fixed — a `GroupBy` whose third argument was not a `Function` (so the grouped
  result always had `Length` 1), and `Sort` on an `Association`, which orders by value and so
  never canonicalises key order for `===`. Any state-sufficiency check should include a control
  state known to be insufficient, and confirm the checker reports it as such.
