# Sampling facet-labeled pure complexes, exactly

> **What this is.** A specification of `RandomUniformFacetLabeledPureSimplicialComplex`, matching
> the implementation in `Kernel/Subpackages/PureComplexes.wl`. Companion to `Theory.md`, and to
> `UnlabeledCount.md` / `UnlabeledSampler.md`, which do the same for the fully unlabeled case.
> It records what is being sampled, the pair-sampling principle that makes the draw exact, the
> lemma that licenses the sequential facet draw, the two indexing tricks that make it fast, and
> how correctness is established — by exact identities rather than by goodness-of-fit.

Notation: $`p`$ purity (every facet has $`p`$ vertices), $`M`$ facet order (number of facets), $`n`$
vertex count, as in `FacetLabeledCount.md`, whose notation this document follows. $`\lambda \vdash n`$
is an integer partition, used as a cycle type, $`m_k`$ the number of parts of $`\lambda`$ equal to
$`k`$ (so $`\lambda`$ is also written as its multiplicity vector $`(m_1,\dots,m_p)`$), and
$`z_\lambda = \prod_k k^{m_k} m_k!`$, so $`n!/z_\lambda`$ permutations of $`[n]`$ have that type.
$`s_F(p,M,n)`$ is `NumFacetLabeledPureComplexes[p,M,n]`, $`x^{(r)} = x(x-1)\cdots(x-r+1)`$ is the
falling factorial, and $`P(k)`$ is the number of integer partitions of $`k`$. $`[n]`$ is
$`\{1,\dots,n\}`$; separately, $`[x^k]f`$ denotes the coefficient of $`x^k`$ in a series $`f`$.

---

## 1. What is being sampled

A **facet-labeled pure complex** is an ordered $`M`$-tuple of pairwise distinct $`p`$-subsets of a
vertex set, covering it, **up to relabelling the vertices**. The facets carry labels
$`1,\dots,M`$; the vertices do not. Formally, let

```math
X \;=\; \{(F_1,\dots,F_M) : F_i \subseteq [n],\ |F_i| = p,\ F_i \neq F_j,\ \textstyle\bigcup_i F_i = [n]\},
```

on which $`S_n`$ acts by relabelling. The objects are the **orbits**, counted by $`s_F(p,M,n)`$.

Uniform means each of the $`s_F(p,M,n)`$ orbits has probability exactly $`1/s_F(p,M,n)`$ — **not**
that each element of $`X`$ is equally likely. The two differ whenever orbits have different sizes,
i.e. whenever some complex has a nontrivial automorphism.

The output is a facet list on vertex set $`[n]`$: a representative of its orbit, whose labelling is
itself uniform among the $`n!/|\mathrm{Aut}|`$ labellings of that class.

### The incidence-tableau picture

Transposing the incidence matrix gives the picture the counter uses, and the one the DP
alternative in §11 is stated in. A complex becomes a **multiset of $`n`$ nonempty subsets of
$`[M]`$** — one row per vertex, listing the facets it lies in — in which every label occurs in
exactly $`p`$ rows. The rows form a multiset because the vertices are unlabelled; "the $`M`$ facets
are pairwise distinct" becomes "the $`M`$ columns are pairwise distinct", the **separating**
condition. This is the form in which $`s_F`$ is derived (see the header comment at
`PureComplexes.wl:1801`), and it is a complete invariant: two facet-labeled complexes are
isomorphic iff their row multisets agree. That is what the test suite canonicalises with.

---

## 2. Why the shipped sampler is not a port of the vertex-labeled one

`RandomVertexLabeledPureSimplicialComplex` reads its weights straight off the counting recursion,
because with labelled vertices "how many complexes put $`k`$ new vertices in the next facet" is a
plain product of binomials, and the count decomposes facet by facet.

Nothing like that was available when this sampler was written. With the vertices unlabelled, $`s_F`$
is a **Burnside average over cycle types**, and the pair-sampling scheme of §3 is what replaced the
missing recursion.

> **Earlier revisions said no such recursion exists and that looking for one was wasted effort.
> That was too strong, and it is now false.** `FacetLabeledCount.md` §9.3 gives a genuine
> facet-by-facet recursion: it carries a **partition** of the vertex count — the multiset of
> type multiplicities — and never forms a cycle type. The obstruction the old wording described is
> real only for a recursion in $`(p,M,n)`$ alone; carrying that partition removes it, exactly as
> `FacetLabeledCount.md` §9 now says for the counting side.
>
> **And it does drive an exact sampler.** Built and verified 2026-09-14: uniform over isomorphism
> classes on eight parameter sets, every class hit, $`\chi^2`$ p-values 0.09–0.53, class sets
> identical to this sampler's. The decomposition is positive because the subtraction in (9.3)
> removes $`M-j`$ *concrete* forbidden branches that all land on the same child state, so branch
> probabilities exist; the DP table stays indexed by the partition while the concrete complex rides
> alongside to say which all-or-nothing branches are the already-placed facets.
>
> **It is not currently a replacement**, on speed: measured 0.94× at $`(2,3,4)`$ but 8× at
> $`(3,6,12)`$, 19× at $`(4,6,15)`$ and 47× at $`(3,8,16)`$, widening with $`p`$ and $`M`$. What it offers is
> **independence**. `FacetLabeledCount.md` §12 notes that this sampler's strongest leg uses the
> shipped counter as its right-hand side, so an error there would be masked rather than caught; a
> sampler that is Burnside-free on both sides is the first that could disagree.

The older `RandomFacetLabeledPureSimplicialComplex` (no `Uniform` in the name) does exactly that
anyway, by iterative vertex addition off a shrinking pool. It is **not uniform**, which its
source comment says and its usage message does not. It is a different function and is not
specified here.

---

## 3. The principle: pair sampling

Let

```math
\Phi \;=\; \{(\sigma, x) : \sigma \in S_n,\ x \in X,\ \sigma \cdot x = x\}
```

be the set of **fixed pairs**. For any orbit $`O`$, the orbit–stabiliser theorem gives

```math
\#\{(\sigma,x) \in \Phi : x \in O\} \;=\; \sum_{x \in O} |\mathrm{Stab}(x)| \;=\; |O| \cdot \frac{n!}{|O|} \;=\; n!,
```

**the same number for every orbit, whatever its size.** So drawing $`(\sigma,x)`$ uniformly from
$`\Phi`$ and discarding $`\sigma`$ leaves the orbit exactly uniform. No rejection, and no
automorphism-order weighting to get wrong.

This is the double-counting that proves the Cauchy–Frobenius (Burnside) lemma [1,2], read as a
sampling recipe rather than a counting one; summing the display over orbits gives
$`|\Phi| = \sum_\sigma |\mathrm{Fix}(\sigma)| = n!\,s_F(p,M,n)`$, which is the identity §9 checks. The
same observation underlies Jerrum's Burnside process [3,4] — the difference is that the Burnside
process is a Markov chain that converges to the orbit-uniform distribution, whereas here $`\Phi`$ is
simple enough to sample from **directly and exactly**, in one shot.

The identical argument makes `RandomUniformUnlabeledPureSimplicialComplex` exact; only the inner
object changes, from an ordered tuple to a $`\sigma`$-invariant set (`UnlabeledSampler.md`).

---

## 4. Which $`\sigma`$ contribute

$`\sigma`$ fixes a tuple $`x`$ exactly when every facet $`F_i`$ is $`\sigma`$-invariant, i.e. is a union
of cycles of $`\sigma`$. Calling the cycles **blocks**, a facet is a set of blocks whose sizes sum
to $`p`$.

Two consequences:

- **Only cycle types with all parts $`\le p`$ occur.** A cycle longer than $`p`$ fits inside no facet,
  so it could never be covered, and $`|\mathrm{Fix}(\sigma)| = 0`$. This is why the type list is
  `IntegerPartitions[n, All, Range[p]]` and not all $`\lambda \vdash n`$.
- **Only the cycle *supports* matter, never the cyclic order.** Given the block partition, the
  number of $`\sigma`$ realising it is $`\prod_k ((k-1)!)^{m_k}`$, a constant depending only on
  $`\lambda`$. So sampling $`\sigma`$ uniformly among its type is the same as sampling the block
  partition uniformly, and the implementation only ever materialises blocks.

Writing $`N(\lambda,p)`$ for the number of admissible facets available to a block collection of type
$`\lambda`$, with $`m_k`$ blocks of size $`k`$,

```math
N(\lambda,p) \;=\; [x^p] \prod_{k \ge 1} (1 + x^k)^{m_k},
```

counting the sub-collections of **blocks** whose sizes total $`p`$ — a truncated cycle-index-style
generating function [5,6]. Blocks of equal size are distinguishable, which is what the exponents
$`m_k`$ supply: expanded, $`N(\lambda,p) = \sum_{\nu} \prod_k \binom{m_k}{\nu_k}`$ over
$`\nu`$ with $`\sum_k k\,\nu_k = p`$. It is *not* the number of sub-multisets of the block
*sizes*: for $`\lambda = \{2,1,1,1\}`$ — three blocks of size 1 and one of size 2 — only two size-multisets
total $`p = 2`$, namely $`\{2\}`$ and $`\{1,1\}`$, while $`N = 4`$, because the three unit blocks offer
$`\binom{3}{2} = 3`$ distinct pairs. `FacetLabeledCount.md` §3.3 states the same count for the same
reason. `RandFLPCWeightCounts` builds the whole coefficient list out to $`x^p`$
by truncated convolution, so the degree-$`n`$ polynomial is never formed.

---

## 5. The four steps

```
1.  λ  ~  ∝ (n!/z_λ) · #{tuples λ fixes}        cycle type of σ
2.  cut a uniform permutation of [n] into blocks of sizes λ
3.  draw the M distinct weight-p block-subsets one at a time,
    each weighted by how many ways the rest can still be finished
4.  facet i = union of the blocks in subset i
```

Step 2 is a Fisher–Yates shuffle [7,8] (`RandomSample`) followed by `TakeList`. Cutting a uniform
permutation word into consecutive blocks of the prescribed sizes realises each block partition
exactly $`z_\lambda`$ times, hence uniformly.

Step 3 is the recursive method for random generation [9,10]: at each stage choose a branch with
probability proportional to the number of completions below it, so that every complete object is
reached along exactly one path with probability $`1/(\text{total count})`$.

---

## 6. Completion counts

The quantity step 3 needs is: having spent $`j`$ subsets and with a set $`U`$ of blocks still
uncovered, how many ways are there to choose $`M - j`$ further weight-$`p`$ subsets, pairwise
distinct, distinct from the $`j`$ already spent, covering all of $`U`$?

Inclusion–exclusion over which uncovered blocks are *left* uncovered gives

```math
C(\lambda, M, j, \kappa) \;=\; \sum_{\omega \le \kappa} (-1)^{|\omega|} \left(\prod_k \binom{\kappa_k}{\omega_k}\right) \Big(N(\lambda - \omega, p) - j\Big)^{(M-j)},
```

where $`\kappa`$ is the multiplicity vector of the sizes of the blocks in $`U`$, $`\omega`$ runs
componentwise over $`0 \le \omega_k \le \kappa_k`$, $`|\omega| = \sum_k \omega_k`$, and
$`\lambda - \omega`$ is the type left after removing $`\omega_k`$ blocks of size $`k`$
(`RandFLPCCompletions`). The falling factorial counts ordered choices of distinct subsets, which
is what "the facets are labelled and pairwise distinct" asks for.

### 6.1 The lemma that makes step 3 exact

The $`-j`$ in $`\big(N(\lambda-\omega,p) - j\big)`$ is the whole subtlety, and it is exact for a reason
worth stating:

> **The number of completions depends on the draws so far only through $`(j, U)`$ — how *many*
> subsets have been used and which blocks remain uncovered — never on *which* subsets were used.**

Every used subset lies inside the already-covered blocks, so it automatically avoids any
$`\omega \subseteq U`$, and therefore drops out of each inclusion–exclusion term by itself; only its
*count* enters. Without this the state would have to remember the used set and the draw would not
decompose. It is the same observation that lets the vertex-labeled generator subtract $`(j-1)`$ for
a repeated facet.

Note this lemma is specific to the facet-labeled case. Its analogue **fails** in the unlabeled
sampler, where a naive port produced plausible complexes from a subtly wrong distribution
(`UnlabeledSampler.md` §7).

### 6.2 The $`n^p`$ indexing

The state is indexed by the **multiplicity vector of uncovered block sizes**, not by the uncovered
block *set*. Blocks of equal size are interchangeable at this point, so nothing is lost, and a
$`2^{\#\text{blocks}}`$ inclusion–exclusion becomes $`\prod_k(\kappa_k + 1) = O(n^p)`$ terms — the same
polynomial scaling, and the same near-independence of $`M`$, that
`NumFacetLabeledPureComplexes` has.

### 6.3 Step-1 weights

```math
w(\lambda) \;=\; \frac{n!}{z_\lambda}\, C(\lambda, M, 0, \lambda),
```

the number of permutations of type $`\lambda`$ times the number of tuples each fixes
(`RandFLPCTypeWeights`). By §3 these must satisfy

```math
\sum_{\lambda} w(\lambda) \;=\; n!\; s_F(p,M,n). \tag{I1}
```

---

## 7. The grouped draw

Steps 3's naive form weighs every candidate subset separately: at $`\{p,M,n\} = \{4,6,15\}`$ that is
441 candidates $`\times`$ 6 facets, and it was ~100% of the cost of a draw.

But a candidate's weight enters **only** through $`\kappa\big(U \setminus S\big)`$, the uncovered
blocks it would leave, and with $`U`$ fixed at this stage that is
$`\kappa(U) - \kappa(S \cap U)`$. So the weight is determined by the **key** $`\kappa(S \cap U)`$, the
multiplicity vector of the still-uncovered blocks that $`S`$ itself contains, which has
$`\sum_k k\,\kappa_k(S \cap U) \le p`$ because $`|S| = p`$. The number of distinct such keys is at most

```math
\sum_{i=0}^{p} P(i) \;=\; 4,\ 7,\ 12 \quad\text{for } p = 2,3,4,
```

**however many hundreds of candidates there are** ($`P`$ the partition function). At $`\{4,6,15\}`$
the 441 candidates carry 3 distinct keys.

So the implementation groups candidates by that key and draws

```math
\Pr[\text{group } g] \;\propto\; |g| \cdot C\big(\lambda, M, j{+}1, \kappa(U) - \text{key}(g)\big),
```

then a member of $`g`$ uniformly. Members of a group have equal weight, so this is the same
distribution — but with a handful of memo lookups per facet instead of one per candidate.

Two supporting details:

- **The key is carried as an integer.** Each entry of $`\kappa(S \cap U)`$ is at most $`p`$, so base
  $`p+1`$ packs the vector without carrying: $`\mathrm{code} = \sum_k \kappa_k(S \cap U)\,(p+1)^{k-1}`$, and `Tally` on a packed integer array does
  the grouping in one pass.
- **Codes are updated in place, not recomputed.** When a facet covers block $`b`$ for the first
  time, exactly the candidates containing $`b`$ lose it from their intersection with $`U`$, so
  `codes[[containing[[b]]]] -= radix[[sizes[[b]]]]`. The blocks-to-candidates incidence, the
  candidate list, and the initial codes depend only on the cycle type, so `RandFLPCCandTable`
  memoises them; there are only $`P_{\le p}(n)`$ cycle types.
- Spent candidates are marked with the code $`-1`$, which no real key can take.

A wrong key here would still emit perfectly plausible complexes from the wrong distribution, so
the test suite checks the two invariants directly (§10).

### A worked draw at $`\{p,M,n\} = \{3,3,5\}`$

Small enough to follow by hand, but every mechanism above fires. A nontrivial cycle type carries
real weight, the grouped draw offers more than one group at two of its three stages, and the
classes have different stabilisers, which is exactly where "uniform over tuples" and "uniform over
classes" part company. It is also the parameter set of the chi-square test in §12. Every number
below comes from the shipped helpers (`RandFLPCTypeWeights`, `RandFLPCCandTable`,
`RandFLPCCompletions`), and the draw is a seeded run of the same logic as `RandFLPCOne`.

**The target.** $`s_F(3,3,5) = 7`$. There are $`3! \cdot 100 = 600`$ tuples in $`X`$, since 100 is the
vertex-labeled count. Four of the seven classes have a stabiliser of order 2, and so 60 tuples each;
the other three have a trivial stabiliser and 120 tuples each. A uniform tuple would therefore put
$`1/10`$ on each of the first four classes and $`1/5`$ on each of the last three. The sampler must put
$`1/7`$ on every one.

**Step 1: the cycle type.** Only types with all parts $`\le 3`$ occur (§4):

| $`\lambda`$ | $`n!/z_\lambda`$ | $`N(\lambda,3)`$ | $`C(\lambda,3,0,\lambda)`$ | $`w(\lambda)`$ |
| --- | --- | --- | --- | --- |
| $`(3,2)`$ | 20 | 1 | 0 | 0 |
| $`(3,1,1)`$ | 20 | 1 | 0 | 0 |
| $`(2,2,1)`$ | 15 | 2 | 0 | 0 |
| $`(2,1,1,1)`$ | 10 | 4 | 24 | 240 |
| $`(1^5)`$ | 1 | 10 | 600 | 600 |

The first three types offer fewer than $`M = 3`$ admissible facets, so every falling factorial in
their inclusion–exclusion vanishes. The weights sum to $`840 = 5! \cdot 7`$, which is (I1). The draw
takes $`(2,1,1,1)`$ with probability $`240/840 = 2/7`$ and the identity type with probability $`5/7`$.
Suppose it takes $`(2,1,1,1)`$.

**Step 2: the blocks.** The shuffle gives $`(1,4,3,2,5)`$, and cutting it to the sizes
$`(2,1,1,1)`$ gives the blocks $`b_1 = \{1,4\}`$, $`b_2 = \{3\}`$, $`b_3 = \{2\}`$, $`b_4 = \{5\}`$. This is
$`\sigma = (1\,4)`$; the cyclic order inside a block never matters (§4).

**The candidates.** A candidate is a set of blocks of total size 3, so there are $`N(\lambda,3) = 4`$
of them. Their keys are packed in base $`p + 1 = 4`$, as
$`\mathrm{code} = \kappa_1 + 4\kappa_2 + 16\kappa_3`$:

| candidate | blocks | facet | key $`\kappa(S \cap U)`$ | code |
| --- | --- | --- | --- | --- |
| $`c_1`$ | $`b_1, b_2`$ | $`\{1,3,4\}`$ | $`(1,1,0)`$ | 5 |
| $`c_2`$ | $`b_1, b_3`$ | $`\{1,2,4\}`$ | $`(1,1,0)`$ | 5 |
| $`c_3`$ | $`b_1, b_4`$ | $`\{1,4,5\}`$ | $`(1,1,0)`$ | 5 |
| $`c_4`$ | $`b_2, b_3, b_4`$ | $`\{2,3,5\}`$ | $`(3,0,0)`$ | 3 |

Any three distinct candidates cover all four blocks, so all $`4 \cdot 3 \cdot 2 = 24`$ ordered choices
are $`\sigma`$-fixed tuples, which is the $`C(\lambda,3,0,\lambda) = 24`$ of step 1.

**Step 3, first facet ($`j = 0`$).** Nothing is covered yet, so $`\kappa(U) = (3,1,0)`$:

| group (code) | members | $`\lvert g \rvert`$ | $`\kappa(U) - \text{key}`$ | $`C(\lambda,3,1,\cdot)`$ | weight | probability |
| --- | --- | --- | --- | --- | --- | --- |
| 5 | $`c_1, c_2, c_3`$ | 3 | $`(2,0,0)`$ | 6 | 18 | $`3/4`$ |
| 3 | $`c_4`$ | 1 | $`(0,1,0)`$ | 6 | 6 | $`1/4`$ |

The weights sum to $`24 = C(\lambda,3,0,\lambda)`$, as the recursive method requires. The draw takes
group 5 and then, uniformly within it, $`c_2`$: so $`F_1 = \{1,2,4\}`$. This covers $`b_1`$ and $`b_3`$ for
the first time. The candidates containing $`b_1`$ ($`c_1, c_2, c_3`$) lose the place value 4, those
containing $`b_3`$ ($`c_2, c_4`$) lose 1, and $`c_2`$ is marked spent. The codes are now
$`(1, \text{spent}, 1, 2)`$.

**Step 3, second facet ($`j = 1`$).** $`b_2`$ and $`b_4`$ are uncovered, so $`\kappa(U) = (2,0,0)`$:

| group (code) | members | $`\lvert g \rvert`$ | $`\kappa(U) - \text{key}`$ | $`C(\lambda,3,2,\cdot)`$ | weight | probability |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | $`c_1, c_3`$ | 2 | $`(1,0,0)`$ | 2 | 4 | $`2/3`$ |
| 2 | $`c_4`$ | 1 | $`(0,0,0)`$ | 2 | 2 | $`1/3`$ |

The weights sum to $`6 = C(\lambda,3,1,(2,0,0))`$. The draw takes group 2, so $`F_2 = \{2,3,5\}`$. This
covers $`b_2`$ and $`b_4`$, and every remaining code drops to 0.

**Step 3, third facet ($`j = 2`$).** Everything is covered, and the one group holds $`c_1`$ and $`c_3`$,
each with $`C(\lambda,3,3,0) = 1`$. The draw takes $`c_1`$ uniformly, so $`F_3 = \{1,3,4\}`$.

Given the cycle type and the blocks, this path had probability
$`\tfrac34 \cdot \tfrac13 \cdot \tfrac13 \cdot \tfrac12 = \tfrac1{24} = 1/C(\lambda,3,0,\lambda)`$. So does
every other path: each of the 24 $`\sigma`$-fixed tuples is equally likely.

**Step 4: the output** is $`(F_1,F_2,F_3) = (\{1,2,4\}, \{2,3,5\}, \{1,3,4\})`$: $`F_1`$ and $`F_3`$ share
the edge $`\{1,4\}`$, and $`F_2`$ meets each of them in one vertex. Every facet is a union of blocks, so
$`\sigma = (1\,4)`$ fixes the tuple. The output therefore lies in a class with a stabiliser of order 2,
as every draw from a nontrivial cycle type must.

**Why every class comes out at exactly $`1/7`$.** Enumerating every shuffle and every branch of
step 3, for both cycle types, gives each tuple's probability exactly:

| $`\lvert\mathrm{Stab}\rvert`$ | via $`(1^5)`$ | via $`(2,1,1,1)`$ | per tuple | tuples per class | per class |
| --- | --- | --- | --- | --- | --- |
| 1 | $`\tfrac57 \cdot \tfrac1{600} = \tfrac1{840}`$ | 0 | $`\tfrac1{840}`$ | 120 | $`\tfrac17`$ |
| 2 | $`\tfrac1{840}`$ | $`\tfrac27 \cdot \tfrac1{240} = \tfrac1{840}`$ | $`\tfrac1{420}`$ | 60 | $`\tfrac17`$ |

The identity branch alone is a uniform tuple, which would give the $`1/10`$ and $`1/5`$ above. The
$`(2,1,1,1)`$ branch adds exactly one more pair per tuple on the classes that have an automorphism:
the 240 pairs of that slice of $`\Phi`$ are the $`4 \times 60`$ tuples with a stabiliser of order 2,
each with its one non-identity stabiliser element. This is the orbit–stabiliser count of §3, in
numbers: each tuple has probability $`|\mathrm{Stab}(x)|/(n!\,s_F) = |\mathrm{Stab}(x)|/840`$.

---

## 8. Vertex count free

`RandomUniformFacetLabeledPureSimplicialComplex[{p,M}, k]` draws $`n`$ from the facet-labeled counts
themselves, $`\Pr[n] \propto s_F(p,M,n)`$ over $`p \le n \le pM`$, then samples at that $`n`$. The
previous implementation drew $`n`$ from the **vertex-labeled** counts and let a rejection step
repair the difference; weighting by $`s_F`$ is the distribution that was wanted all along, and needs
no repair.

The emptiness guard tests $`s_F(p,M,n) \neq 0`$, not the vertex-labeled count. The two vanish
together on $`p \le n \le pM`$, but $`s_F`$ also needs $`\binom{n}{p} \ge M`$ distinct facets, so it is
the correct test.

---

## 9. Why it is uniform

Three independent legs, in decreasing order of how much they would catch:

1. **(I1)** $`\sum_\lambda w(\lambda) = n!\,s_F(p,M,n)`$, exactly, as integers. A wrong cycle-type
   weight still produces plausible-looking complexes, so this identity — not the outputs — is the
   thing to check.
2. **The completion-count invariants** of §7: the candidate table lists exactly the admissible
   facets, and the incrementally maintained codes still decode to $`\kappa(S \cap U)`$, so that
   $`\kappa(U)`$ minus the decoded key is $`\kappa(U \setminus S)`$, after any sequence of facets has
   been laid down.
3. **Goodness-of-fit** over exhaustively enumerated classes, as a backstop. At 60,000–100,000
   draws on $`\{2,4,5\}`$, $`\{3,3,5\}`$ and $`\{4,3,8\}`$ the draw is uniform and statistically
   indistinguishable from an independently derived exact sampler.

Chi-square is the weak leg. It is insensitive to exactly the errors that matter, and repeated runs
throw low p-values by chance; prefer legs 1 and 2.

---

## 10. Implementation map

All in `Kernel/Subpackages/PureComplexes.wl`; the derivation is the header comment at line 4002.
Private helpers are named `RandFLPC*`.

| line | symbol | role |
| --- | --- | --- |
| 4051 | `RandFLPCWeightCounts[m]` | coefficients of $`\prod_k (1+x^k)^{m_k}`$ out to $`x^p`$; last entry is $`N(\lambda,p)`$ (§4) |
| 4067 | `RandFLPCCompletions[λ,M,j,κ]` | $`C(\lambda,M,j,\kappa)`$, the inclusion–exclusion of §6 |
| 4082 | `RandFLPCTypeWeights[p,M,n]` | step-1 weights $`w(\lambda)`$, with the types and their multiplicity vectors |
| 4111, 4114 | `$RandomFacetLabeledParallelThreshold`, `RandFLPCGoParallel` | sample-count gate, also gated on `$KernelCount` |
| 4118 | `RandFLPCCandTable[sizes,p]` | per-cycle-type candidate list, packed keys, place values, block→candidate incidence (§7) |
| 4140 | `RandFLPCOne[p,M,n]` | one sample, steps 1–4 |
| 4191 | primary pattern | `[{p,M,n}, numSamples]` |
| 4235 | overload | `[{p,M}, numSamples]` (§8) |
| 4278, 4287 | overloads | single-sample forms |

The counter it sits on is `NumFacetLabeledPureComplexes` at line 2044, with its own derivation at
line 1861 and helpers `NumFLPC*`.

All memo tables are released by the shared ``ECGrav`Private`NumPCClearCache[]`` — which must list
every memoised helper, `RandFLPCCandTable` included.

---

## 11. Cost

Measured warm in a fresh kernel, against the previous body extracted from `39dc9d3`. "Before" is
the per-candidate weighing; "after" is the grouped draw of §7.

| $`\{p,M,n\}`$ | before | after | speedup |
| --- | --- | --- | --- |
| $`\{2,3,4\}`$ | 0.094 ms | 0.066 ms | 1.4× |
| $`\{2,4,6\}`$ | 0.189 ms | 0.089 ms | 2.1× |
| $`\{3,4,8\}`$ | 0.553 ms | 0.099 ms | 5.6× |
| $`\{2,8,12\}`$ | 1.202 ms | 0.163 ms | 7.4× |
| $`\{3,6,12\}`$ | 2.688 ms | 0.153 ms | 17.6× |
| $`\{3,7,14\}`$ | 4.819 ms | 0.190 ms | 25.4× |
| $`\{4,6,15\}`$ | 14.42 ms | 0.280 ms | 51.5× |
| $`\{3,10,20\}`$ | 19.72 ms | 0.328 ms | 60.1× |
| $`\{2,20,30\}`$ | 18.97 ms | 0.516 ms | 36.8× |
| $`\{3,12,24\}`$ | 46.53 ms | 0.606 ms | 76.8× |
| $`\{4,8,20\}`$ | 52.80 ms | 0.513 ms | **102.9×** |
| $`\{3,15,30\}`$ | 108.5 ms | 1.240 ms | 87.5× |

The first call at a given $`(p,M,n)`$ builds the memo tables; the figures are warm, which is the
regime any real use is in.

**`$RandomFacetLabeledParallelThreshold = 2000`.** It was 100, which suited a draw costing 3–60 ms.
What the parallel branch pays is fixed rather than per-sample — `DistributedContexts` ships the
whole of `` ECGrav`Private` ``, memo tables included, on every call — so cutting the per-draw cost by
20–100× moved the crossover out with it, and at 100 samples parallel had become a 3–16×
*pessimisation*. Remeasured on 11 subkernels with kernels up and tables warm, the crossover runs
200 at $`\{2,3,4\}`$, 300 at $`\{2,4,6\}`$, 800 at $`\{3,6,12\}`$, 1700 at $`\{4,6,15\}`$, 2800 at
$`\{3,10,20\}`$ and past 3200 at $`\{4,8,20\}`$. No single constant is right everywhere and the
penalties are lopsided — too low costs up to 16×, too high costs at most the parallel speedup,
only 2–3× anywhere near the crossover — so the constant sits at the top of the range. At 2000 the
serial/parallel ratios run 0.87–2.34; by 5000, 1.49–3.48.

### The alternative that was measured and rejected

A rejection-free **dynamic program over row-type multiplicities** carrying the current
column-equality partition $`\pi`$ as part of its state — conditioning on separation instead of
testing for it — is correct (its counts agree with $`s_F`$ on all 234 triples with $`p \le 4`$,
$`M \le 6`$) but is exponential where this sampler is polynomial. Its state space is
$`2^M (n+1) (p+1)^M \,\mathrm{Bell}(M)`$ with $`\mathrm{Bell}(M)`$ the Bell number, and since $`n \le pM`$ the only free axes are $`p`$
and $`M`$, both exponential.
Building the table costs 9.7 s / 164 MB at $`\{3,6,12\}`$, 118 s / 1.8 GB at $`\{3,7,14\}`$, and did
not finish in 180 s / 4 GB at $`\{3,8,16\}`$ — where this sampler draws in 0.27 ms. The reason it
loses as a *counter* is structural: the $`\big(N\big)^{(M)}`$ falling factorial in
`NumFacetLabeledPureComplexes` already performs the Möbius inversion over the partition lattice
[11,12] in closed form, which is exactly what the $`\pi`$ state re-derives by brute force.

---

## 12. Tests

In `Tests/PureComplexes.wlt`:

| TestID | what it pins |
| --- | --- |
| `RandomUniformFacetLabeled-burnside-weights` | **(I1)**, exactly, on 10 parameter sets |
| `RandFLPCCandTable-keys-track-uncovered` | the §7 invariants, over every cycle type with $`p \le 4`$, $`n \le 11`$ |
| `RandomUniformFacetLabeled-uniform-over-classes` | chi-square over the 7 classes of $`\{3,3,5\}`$, pinned to the serial branch |
| `RandomUniformFacetLabeled-vertex-count-marginal` | the `{p,M}` form draws $`n`$ from $`s_F`$ (§8) |
| `RandomUniformSamplers-structure`, `-parallel-branch`, `-empty-space` | shape, the parallel path, and the guards of §8 |

The counter has its own seven, including a brute-force enumeration of the defining description and
the Stirling identity $`B(p,M,n) = \sum_k S(M,k)\, s_F(p,k,n)`$, which exercises the derivation rather
than its endpoints.

**Non-vacuity.** The §7 invariant test was confirmed to fail under two injected defects — a wrong
packing radix, and codes never updated — before being trusted. A differential test that cannot
fail is worse than no test; this repo has been bitten by that twice
(`UnlabeledSampler.md` §15).

---

## 13. References

The algorithm is assembled from standard results; these are the ones it actually rests on.

1. W. Burnside, *Theory of Groups of Finite Order*, 2nd ed., Cambridge University Press, 1911,
   §145. — The orbit-counting lemma. Due to Cauchy and Frobenius, not Burnside; see
   P. M. Neumann, "A lemma that is not Burnside's", *Math. Sci.* **4** (1979), 133–141.
2. R. P. Stanley, *Enumerative Combinatorics*, Vol. 1, 2nd ed., Cambridge University Press, 2012,
   Ch. 1–2. — The falling factorial and the Stirling identity $`x^n = \sum_k S(n,k)\,x^{(k)}`$ that
   imposes "pairwise distinct facets" (§6), and inclusion–exclusion.
3. M. Jerrum, "Uniform sampling modulo a group of symmetries using Markov chain simulation",
   *DIMACS Ser. Discrete Math. Theoret. Comput. Sci.* **10** (1993), 37–47. — The Burnside
   process: the same fixed-pair set $`\Phi`$, sampled by Markov chain where §3 samples it directly.
4. L. A. Goldberg and M. Jerrum, "The 'Burnside process' converges slowly", *Combin. Probab.
   Comput.* **11** (2002), 21–34. — Why the chain is not a shortcut, and direct sampling from $`\Phi`$
   is worth the effort when $`\Phi`$ is tractable.
5. G. Pólya, "Kombinatorische Anzahlbestimmungen für Gruppen, Graphen und chemische
   Verbindungen", *Acta Math.* **68** (1937), 145–254; English translation in G. Pólya and
   R. C. Read, *Combinatorial Enumeration of Groups, Graphs, and Chemical Compounds*,
   Springer, 1987. See also Stanley, *Enumerative Combinatorics*, Vol. 2, §7.24, for cycle
   indicators in modern notation. — Cycle-index enumeration; the generating function
   $`\prod_k (1+x^k)^{m_k}`$ of §4 is a truncated instance.
6. J. H. Redfield, "The theory of group-reduced distributions", *Amer. J. Math.* **49** (1927),
   433–455. — The earlier, independent statement of the same theory.
7. R. A. Fisher and F. Yates, *Statistical Tables for Biological, Agricultural and Medical
   Research*, Oliver & Boyd, 1938, Example 12. — The shuffle behind step 2.
8. D. E. Knuth, *The Art of Computer Programming*, Vol. 2, 3rd ed., Addison-Wesley, 1997,
   §3.4.2, Algorithm P. — The in-place form, and the correctness argument for cutting a uniform
   permutation into blocks.
9. H. S. Wilf, "A unified setting for sequencing, ranking, and selection algorithms for
    combinatorial objects", *Adv. Math.* **24** (1977), 281–291; and A. Nijenhuis and H. S. Wilf,
    *Combinatorial Algorithms*, 2nd ed., Academic Press, 1978. — The recursive method: sequential
    choice with probability proportional to completion counts, which is step 3.
10. P. Flajolet, P. Zimmermann and B. Van Cutsem, "A calculus for the random generation of
    labelled combinatorial structures", *Theoret. Comput. Sci.* **132** (1994), 1–35. — The
    systematic form of the same method, and the cost model for it.
11. G.-C. Rota, "On the foundations of combinatorial theory I: Theory of Möbius functions",
    *Z. Wahrscheinlichkeitstheorie* **2** (1964), 340–368. — Möbius inversion over the partition
    lattice, with $`\mu(\hat 0, \pi) = \prod_B (-1)^{|B|-1}(|B|-1)!`$; the closed form that the
    falling-factorial weight encodes (§11).
12. R. P. Stanley, *Enumerative Combinatorics*, Vol. 1, 2nd ed., Ch. 3. — The partition lattice
    $`\Pi_M`$, its Möbius function, and the Bell-number growth that sinks the DP alternative.

Companion specifications in this repository: `Theory.md` (the model),
`UnlabeledCount.md` and `UnlabeledSampler.md` (the fully unlabeled case, where the §6.1 lemma
fails and a different scheme is needed).
