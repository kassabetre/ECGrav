# Counting fully unlabeled pure complexes

> **What this is.** A specification of `NumUnlabeledPureComplexes`, matching the
> implementation in `Kernel/Subpackages/PureComplexes.wl`. It records the derivation, the
> algorithm actually shipped, the conventions, what it was verified against, and what was
> ruled out. It is written to be the starting point for the uniform sampler, which needs the
> same decomposition for its weights.

Notation: $p$ is the **purity** (every facet has $p$ vertices), $M$ the **facet order**
(number of facets), $n$ the **vertex count**. The argument order is $(p, M, n)$ throughout,
matching `NumVertexLabeledPureComplexes` and `NumFacetLabeledPureComplexes`. $[n]$ is
$\{1,\dots,n\}$, $S_n$ the symmetric group, $\lambda \vdash n$ an integer partition, $m_k$
the number of parts of $\lambda$ equal to $k$, and $z_\lambda = \prod_k k^{m_k} m_k!$, so
that $n!/z_\lambda$ is the number of permutations of cycle type $\lambda$.

---

## 1. The object being counted

A **pure complex** here is a set $S$ of $M$ pairwise **distinct** $p$-subsets of $[n]$ whose
union is all of $[n]$. The facets are distinct (it is a set, not a multiset) and they cover
(no isolated vertices).

Three counts sit on the same objects, differing only in what carries a label:

| function | vertices | facets | counts |
| --- | --- | --- | --- |
| `NumVertexLabeledPureComplexes` | labelled | unlabelled | the sets $S$ themselves |
| `NumFacetLabeledPureComplexes` | unlabelled | labelled | $S_n$-orbits of ordered $M$-tuples |
| `NumUnlabeledPureComplexes` | unlabelled | unlabelled | $S_n$-orbits of the sets $S$ |

The third is the number of **isomorphism classes**, and it is the count the other two refine:
every class carries at least one vertex-labeled and at least one facet-labeled complex, and
distinct classes carry disjoint sets of them, so

$$U(p,M,n) \;\le\; \min\bigl(\mathrm{VL}(p,M,n),\ \mathrm{FL}(p,M,n)\bigr).$$

Note that $\mathrm{FL}$ can *exceed* $\mathrm{VL}$: at $(2,6,4)$ the only complex is all six
edges of $K_4$, giving $\mathrm{VL}=1$, $\mathrm{FL}=30$, $U=1$. Per class the facet-labeled
contribution is $M!/|H|$ with $H$ the **image** of $\mathrm{Aut}(S)$ in $\mathrm{Sym}(S)$ —
not $M!/|\mathrm{Aut}(S)|$, which is only correct when that action is faithful. For $p \ge 3$
a non-identity automorphism can fix every facet setwise, e.g. $(2\,3)$ acting on
$\{1,2,3\},\{2,3,4\}$. The vertex-labeled contribution is $n!/|\mathrm{Aut}(S)|$, which *is*
the plain stabiliser.

---

## 2. Derivation

### 2.1 The group is $S_n$, not $S_n \times S_M$

Acting on the ordered $M$-tuples that the facet-labeled count uses, the natural group is the
product $G = S_n \times S_M$, and $(\sigma,\tau)$ fixes $(F_1,\dots,F_M)$ exactly when
$\sigma(F_i) = F_{\tau(i)}$. Carrying that through would require tracking, for each
$\tau$-cycle of length $\ell$, a single facet $F$ with $\sigma^\ell(F) = F$, hence the cycle
structure of $\sigma^\ell$ and gcd bookkeeping.

**None of that is needed.** Forgetting the facet order maps tuples onto *sets* of facets;
that map commutes with $S_n$; and because the facets in a tuple are pairwise **distinct**,
$S_M$ acts **simply transitively** on each fibre. Hence the $G$-orbits on tuples are exactly
the $S_n$-orbits on sets, and a single Burnside average over $S_n$ suffices, with

$$|\mathrm{Fix}(\sigma)| \;=\; \#\{\sigma\text{-invariant } M\text{-element sets of distinct } p\text{-subsets covering } [n]\}.$$

The distinctness of the facets is what makes the fibres free; it is not a technicality.

### 2.2 Covering, removed and restored

Covering is not a per-orbit condition, so it is taken off before Burnside and put back
afterwards. Let $A(p,M,n)$ be the same count with covering **not** required. An unlabeled
object whose facets span $j$ vertices is an unlabeled *covering* object on $j$ vertices
together with $n-j$ **indistinguishable** isolated vertices — and indistinguishability is
exactly what makes the correspondence one-to-one, which it would not be with labels. So

$$A(p,M,n) = \sum_{j \le n} U(p,M,j), \qquad\text{hence}\qquad U(p,M,n) = A(p,M,n) - A(p,M,n-1).$$

This is the same padding trick `NumFacetLabeledPureComplexes` uses.

### 2.3 Burnside over cycle types

$$A(p,M,n) \;=\; \frac{1}{n!}\sum_{\sigma \in S_n} |\mathrm{Fix}(\sigma)| \;=\; \sum_{\lambda \vdash n} \frac{1}{z_\lambda}\,|\mathrm{Fix}(\lambda)|.$$

The sum over **integer partitions of $n$** is the outer loop of the algorithm, and the reason
$n$ is the limiting argument.

### 2.4 $|\mathrm{Fix}(\sigma)|$ is a coefficient

A $\sigma$-invariant set of $p$-subsets is precisely a union of $\langle\sigma\rangle$-**orbits**
of $p$-subsets. Choosing one is choosing which orbits to include, subject to the sizes
totalling $M$. If $\sigma$ has $n_d$ orbits of size exactly $d$ on the $p$-subsets, then

$$|\mathrm{Fix}(\sigma)| \;=\; [z^M] \prod_{d \ge 1} (1 + z^d)^{n_d}.$$

Only $d \le M$ can matter: a factor $(1+z^d)^{n_d}$ with $d > M$ contributes its constant
term and nothing else below $z^{M+1}$.

### 2.5 Fixed subsets of $\sigma^e$

Let $f(e)$ be the number of $p$-subsets fixed by $\sigma^e$. Since $\sigma^e$ splits each
$k$-cycle of $\sigma$ into $\gcd(k,e)$ cycles of length $k/\gcd(k,e)$, and a fixed subset is
a union of those cycles,

$$f(e) \;=\; [x^p] \prod_k \bigl(1 + x^{\,k/\gcd(k,e)}\bigr)^{\gcd(k,e)\,m_k}.$$

Cycles of $\sigma^e$ longer than $p$ can never be used, so those factors are $1$ up to degree
$p$ and are skipped.

### 2.6 From $f$ to the coefficient, by Newton's identity

The orbit counts $n_d$ *could* be recovered by Möbius inversion over the divisor lattice —
$f(e)$ counts the subsets of period dividing $e$, so $n_d = \frac{1}{d}\sum_{e \mid d} \mu(d/e) f(e)$
— but they are **never formed**. Taking the logarithmic derivative of the product in §2.4
turns it into a recurrence whose coefficients are the $f$ directly. Writing
$c(m) = [z^m]\prod_d (1+z^d)^{n_d}$:

$$m\,c(m) \;=\; \sum_{k=1}^{m} \ell(k)\, c(m-k), \qquad c(0) = 1,$$

$$\ell(k) \;=\; \begin{cases} f(k) & k \text{ odd} \\ f(k) - 2 f(k/2) & k \text{ even.}\end{cases}$$

*Proof sketch.* $\log \prod_d (1+z^d)^{n_d} = \sum_d n_d \sum_{j\ge1} (-1)^{j-1} z^{dj}/j$;
collecting $m = dj$ gives $\sum_m \frac{z^m}{m} \sum_{d \mid m} d\,n_d\,(-1)^{m/d-1}$. Now
$\sum_{d \mid m} d\,n_d = f(m)$, because the orbits whose size divides $m$ consist exactly of
the points fixed by $\sigma^m$. Splitting the sign by the parity of $m/d$ — for $m$ even,
$m/d$ is even iff $d \mid m/2$ — leaves $\ell(m) = f(m) - 2f(m/2)$, and $\ell(m) = f(m)$ for
$m$ odd. Then $z P' = P \cdot \sum_m \ell(m) z^m$ gives the recurrence. $\square$

This is worth more than tidiness: $n_d$ is on the order of $\binom{n}{p}/d$, so the binomial
expansion of $(1+z^d)^{n_d}$ asked for binomial coefficients of enormous arguments, whereas
$\ell(k)$ is bounded by $\binom{n}{p}$.

The recurrence wants $f(k)$ for **every** $k \le M$, not only divisors — but $\sigma^k$ fixes
the same subsets as $\sigma^{\gcd(k,L)}$ where $L = \operatorname{lcm}$ of the cycle lengths,
so $f(k) = f(\gcd(k,L))$ and the number of distinct evaluations is unchanged.

**This recurrence does not lift.** It holds inside a single cycle type. Burnside averages a
*product* over $\sigma$, and the average of a product is not the product of averages, so there
is no corresponding recurrence for $A$ or $U$ — see §7.

---

## 3. Worked example: $A(2,3,4)$

| $\lambda$ | $z_\lambda$ | #perms | $f(e)$ | $n_d$ | $\prod_d(1+z^d)^{n_d}$ | $[z^3]$ | weighted |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $\{4\}$ | 4 | 6 | $f(1){=}0,\ f(2){=}2$ | $n_2{=}1$ | $1+z^2$ | 0 | 0 |
| $\{3,1\}$ | 3 | 8 | $f(1){=}0,\ f(3){=}6$ | $n_3{=}2$ | $1+2z^3$ | 2 | 16 |
| $\{2,2\}$ | 8 | 3 | $f(1){=}2,\ f(2){=}6$ | $n_1{=}2,\ n_2{=}2$ | $1+2z+3z^2+4z^3$ | 4 | 12 |
| $\{2,1,1\}$ | 4 | 6 | $f(1){=}2,\ f(2){=}6$ | $n_1{=}2,\ n_2{=}2$ | $1+2z+3z^2+4z^3$ | 4 | 24 |
| $\{1^4\}$ | 24 | 1 | $f(1){=}6$ | $n_1{=}6$ | $(1+z)^6$ | 20 | 20 |

$A(2,3,4) = 72/24 = 3$ and $A(2,3,3) = 1$, so $U(2,3,4) = 2$ — the path $P_4$ and the star
$K_{1,3}$. (The $n_d$ column is shown for exposition; the shipped code goes from $f$ straight
to $[z^3]$ by §2.6.)

---

## 4. Implementation map

All in `Kernel/Subpackages/PureComplexes.wl`. Private helpers are named `NumULPC*`; the
derivation above appears as the header comment at line 1915.

| line | symbol | role |
| --- | --- | --- |
| 1972 | `NumULPCPolyMul[a,b,deg]` | product of two coefficient lists, truncated at `deg`, via `ListConvolve` |
| 1977 | `NumULPCPowPoly[r,c,deg]` | $(1+z^r)^c$ truncated, stocking only multiples of $r$ |
| 1986 | `NumULPCFixedSubsets[parts,e,p]` | $f(e)$ of §2.5; `parts` is `Tally[λ]`, i.e. `{length, multiplicity}` pairs |
| 2001 | `NumULPCFixSets[parts,p,M]` | $|\mathrm{Fix}(\sigma)|$ by the Newton recurrence of §2.6 |
| 2027 | `NumULPCA[p,M,n]` | $A(p,M,n)$: the cycle-type sum of §2.3. **Memoized** |
| 2041 | `NumULPCCount[p,M,n]` | $A(n) - A(n-1)$ |
| 2046 | `NumUnlabeledPureComplexes[p,M,n]` | guards, then `NumULPCCount` |
| 2062 | `NumUnlabeledPureComplexes[p,M]` | summed over $n$ |
| 2075 | catch-all | `::argerr` and `$Failed` |

Two implementation choices that are exact, not heuristic:

- **Only $d \le M$.** Justified in §2.4; it also caps which $f(e)$ are needed. Checked to
  agree with the unrestricted version on 240 parameter sets.
- **`Divisible[L, d]` over `Range[M]`, not `Divisors[L]`.** The lcm $L$ of the cycle lengths
  is bounded by Landau's function $g(n)$, which grows fast enough to have many more divisors
  than a typical $M$: $g(20) = 420$ (24 divisors), $g(30) = 4620$ (48). Testing $M$ candidates
  beats enumerating every divisor of $L$ and filtering.

`NumULPCA` is memoized for the session and released by the shared
`ECGrav`Private`NumPCClearCache[]`, alongside the vertex-labeled row cache and the sampler
weight tables. Each entry is a single integer.

---

## 5. Conventions and guards

```wolfram
Which[
  p < 0 || M < 0,            0,
  M == 0,                    If[n == 0, 1, 0],
  n < 0 || n < p || n > p*M, 0,
  Binomial[n, p] < M,        0,
  True,                      NumULPCCount[p, M, n]]
```

- $M = 0$ is the empty complex: $1$ at $n = 0$, else $0$.
- Zero unless $p \le n \le pM$ — fewer than $p$ vertices carry no facet, and $M$ facets of
  $p$ vertices cover at most $pM$.
- Zero when $\binom{n}{p} < M$: there are not enough distinct $p$-subsets to supply the facets.

These match `NumVertexLabeledPureComplexes` and `NumFacetLabeledPureComplexes` exactly; the
test asserts agreement with both rather than restating values.

**The guards are load-bearing, not an optimisation.** Falling through to the cycle-type sum on
an out-of-range $n$ calls `IntegerPartitions` on it. `[3,4,100]` must be rejected *before* the
sum is entered — `PartitionsP[100]` is $1.9 \times 10^8$ — which is why the range checks
precede the dispatch rather than living inside `NumULPCCount`. The degenerate-input test
includes $(3,4,100)$ deliberately for this reason.

**The two-argument form** telescopes rather than summing: $U$ vanishes past $n = pM$, so
$\sum_n \bigl(A(n) - A(n-1)\bigr) = A(p,M,pM)$. It therefore costs the three-argument form at
its largest $n$, and is only practical while $\mathrm{PartitionsP}(pM)$ is.

---

## 6. Cost

The outer loop is $\mathrm{PartitionsP}(n)$ cycle types, plus $\mathrm{PartitionsP}(n-1)$ for
the differencing. Per cycle type: $O(M)$ divisibility tests; one $f(e)$ per divisor of $L$
that is $\le M$, each $O(p \cdot \#\text{distinct parts})$; then the $O(M^2)$ Newton
convolution. **$n$ is the limiting argument, not $p$ or $M$.**

Measured cold (cache cleared before each), current implementation:

| $\{p,M,n\}$ | $\mathrm{PartitionsP}(n)$ | seconds |
| --- | --- | --- |
| $\{3,4,7\}$ | 15 | < 0.01 |
| $\{3,6,12\}$ | 77 | 0.02 |
| $\{3,8,18\}$ | 385 | 0.12 |
| $\{2,10,20\}$ | 627 | 0.21 |
| $\{3,10,22\}$ | 1002 | 0.59 |
| $\{3,12,25\}$ | 1958 | 1.19 |
| $\{4,8,28\}$ | 3718 | 1.73 |
| $\{3,12,30\}$ | 5604 | 3.13 |
| $\{2,20,40\}$ | 37338 | 27.3 |
| $\{3,50,10\}$ | 42 | 0.02 |
| $\{3,100,12\}$ | 77 | 0.07 |

Large $M$ at small $n$ is cheap — $\{3,100,12\}$ returns a 56-digit value in 0.07 s. Large $n$
is what hurts, and no amount of tuning changes that: the partition sum is the algorithm.

---

## 7. What is *not* available: a recurrence

`NumVertexLabeledPureComplexes` is computed by a genuine two-index recurrence
(`NumPCAdvance`), which for $C(n) = \binom{n}{p}$ reads

$$q\,N(p,q,n) = \bigl(C(n) - (q-1)\bigr) N(p,q-1,n) + C(n) \sum_{k=1}^{p} \binom{p}{k} N(p,q-1,n-k).$$

There is no analogue here, and it is worth recording why, so it is not re-attempted.

**The structural reason.** That recurrence is a deletion argument: remove a facet, and the
number of ways to put one back depends only on $(q,n)$. On *isomorphism classes* the fibre of
the deletion map is the number of $\mathrm{Aut}(K)$-orbits on candidate facets, which varies
from class to class. Burnside averaging is not compatible with a sequential decomposition.
(`NumFacetLabeledPureComplexes` has no such recurrence either — it is already a cycle-type
sum. Its Stirling relation $B = \sum_k S(M,k) F$ is an identity between two counting problems,
used for verification, not a means of computing $F$.)

**Empirical confirmation.** A linear-algebra ansatz search for
$\sum_{i,j} c_{ij}(M,n)\,U(M-i,n-j) = 0$ with polynomial $c_{ij}$ returned nothing at every
shape where the system was properly overdetermined: up to $I = 2$ previous facet orders,
$J = 4$ vertex shifts, coefficient degrees to $(4,4)$, at $p=2$ and $p=3$, and for the
non-covering $A$ as well. Three guards made the null result meaningful — the search was
**calibrated** on the vertex-labeled count first (returning nullity exactly $1$, surviving 385
held-out points, and decoding to *identically* the recurrence above); systems were required to
be **overdetermined $\ge 3{:}1$**; and any candidate was **fit on half the points and verified
on the other half**.

**A proof, but a narrower one than it looks.** For $T(n) = \sum_M U(2,M,n)$ — unlabeled graphs
on $n$ nodes with no isolated vertex — $\log T(n) = \Theta(n^2)$, while every P-recursive
sequence satisfies $\log|a_n| = O(n \log n)$. So no polynomial-coefficient recurrence **in $n$
alone** exists. But the *vertex-labeled* count has the same $\Theta(n^2)$ growth and a
perfectly good two-index recurrence, so **this argument does not distinguish the two cases**.
It rules out pure-$n$ recurrences and nothing more; the structural reason above is the real
one.

**What does exist** is the recurrence of §2.6, one level down, inside a single cycle type.

---

## 8. Verification

Ten tests in `Tests/PureComplexes.wlt`:

| TestID | what it pins |
| --- | --- |
| `-known-values` | 22 values from the brute-force oracle table (§9) |
| `-brute-force` | in-suite enumeration + canonicalisation over the $n!$ relabellings, 7 parameter sets |
| `-burnside-identity` | the Burnside average from **explicitly enumerated** fixed sets, bypassing §2.4–2.6 |
| `-fix-recurrence` | §2.6 at its own level: `NumULPCFixSets` vs a direct count of $\sigma$-invariant $M$-sets, over 435 (case, cycle type) pairs |
| `-unlabeled-graphs` | external: differenced unlabelled-graph sequence, $n = 0..7$ |
| `-dominated-by-labelled` | $U \le \min(\mathrm{VL}, \mathrm{FL})$ over a 308-point grid, 93 nonzero |
| `-degenerate` | agreement with both siblings on 14 degenerate inputs, including the $(3,4,100)$ guard |
| `-2arg-equals-sum` | the telescoping of §5 |
| `-clear-cache` | memo is reclaimed and the value survives |
| `-argerr` | `$Failed` on malformed input |

Three of these are independent of the derivation rather than of its endpoints:

1. **Brute-force oracle.** Enumerate the covering $M$-sets and group them into $S_n$ orbits
   explicitly. The same enumeration reproduces the vertex- and facet-labeled counts from
   $n!/|\mathrm{Stab}|$ and $M!/|H|$ per class, which is what ties the three together; internal
   guards check that orbit $\times$ stabiliser $= n!$ and that the orbits partition the object
   set. 29 parameter sets, all matching.

2. **Burnside from enumerated fixed sets.** Averages $|\mathrm{Fix}(\sigma)|$ counted by direct
   enumeration, not via the orbit polynomial or the differencing — so it exercises the
   derivation, not just its endpoints.

3. **External, at $p = 2$.** These are the unlabelled graphs with no isolated vertex.
   Differencing the unlabelled-graph sequence $1,1,2,4,11,34,156,1044,\dots$ gives

   $$1,\ 0,\ 1,\ 2,\ 7,\ 23,\ 122,\ 888,\ 11302,\ 262322,\ 11730500,\ 1006992696,\ 164072174728,\ 50336940195360$$

   for $n = 0..13$, matched by $\sum_M U(2,M,n)$. The suite checks to $n = 7$; the full 14
   terms were confirmed during the recurrence work.

Additionally, the Newton rewrite was differentially tested against the Möbius route it
replaced: **61,172 helper-level comparisons** across $p = 1..4$, $M = 0..40$ and every cycle
type of $n = 0..13$, zero mismatches. Since `NumULPCFixSets` is reached only from `NumULPCA`,
that covers every public value.

---

## 9. Reference table

Brute-force values, all three labellings, from one enumeration.

| $p$ | $M$ | $n$ | VL | FL | **U** |
| --- | --- | --- | --- | --- | --- |
| 2 | 2 | 3 | 3 | 1 | **1** |
| 2 | 2 | 4 | 3 | 1 | **1** |
| 2 | 3 | 3 | 1 | 1 | **1** |
| 2 | 3 | 4 | 16 | 4 | **2** |
| 2 | 3 | 5 | 30 | 3 | **1** |
| 2 | 3 | 6 | 15 | 1 | **1** |
| 2 | 4 | 4 | 15 | 15 | **2** |
| 2 | 4 | 5 | 135 | 29 | **4** |
| 2 | 4 | 6 | 330 | 19 | **3** |
| 2 | 5 | 4 | 6 | 30 | **1** |
| 2 | 5 | 5 | 222 | 222 | **5** |
| 2 | 5 | 6 | 1581 | 301 | **9** |
| 2 | 6 | 4 | 1 | 30 | **1** |
| 2 | 6 | 5 | 205 | 1230 | **5** |
| 2 | 6 | 6 | 3760 | 3850 | **15** |
| 3 | 2 | 4 | 6 | 1 | **1** |
| 3 | 2 | 5 | 15 | 1 | **1** |
| 3 | 2 | 6 | 10 | 1 | **1** |
| 3 | 3 | 4 | 4 | 1 | **1** |
| 3 | 3 | 5 | 100 | 7 | **3** |
| 3 | 3 | 6 | 480 | 10 | **3** |
| 3 | 4 | 5 | 205 | 43 | **5** |
| 3 | 4 | 6 | 3600 | 154 | **15** |
| 3 | 4 | 7 | 22820 | 207 | **17** |
| 4 | 2 | 5 | 10 | 1 | **1** |
| 4 | 2 | 6 | 45 | 1 | **1** |
| 4 | 3 | 6 | 395 | 8 | **4** |
| 1 | 3 | 3 | 1 | 1 | **1** |
| 1 | 4 | 4 | 1 | 1 | **1** |

Hand-checks: $U(2,3,4) = 2$ is $P_4$ and $K_{1,3}$; $U(2,4,5) = 4$ is the three trees on five
vertices plus (triangle $\sqcup$ $K_2$) — the latter confirms disconnected complexes are not
being dropped.

---

## 10. Carrying this into the sampler

- **The group collapse of §2.1 should carry verbatim.** Pair-sample $(\sigma, S)$ with $S$ a
  $\sigma$-invariant covering $M$-set: every orbit contributes exactly $n!$ such pairs whatever
  its size, so discarding $\sigma$ leaves the orbit uniform. No $\tau$, no rejection, no
  automorphism weighting. This is the same argument that made
  `RandomUniformFacetLabeledPureSimplicialComplex` exact.
- **The one piece of arithmetic still owed.** This count reaches covering by differencing, so
  it never forms a *per-cycle-type covering count* — but the sampler's step-1 weights need
  exactly that. Getting it means inclusion–exclusion over sub-multisets of $\lambda$: deleting
  cycles leaves the $\sigma$-invariant sets of the restricted permutation.
- **Rows are nearly free.** `NumULPCFixSets` already computes the whole row $c(0..M)$
  internally, at the same $O(M^2)$, and returns only the last entry. A row form — every facet
  order at one vertex count — is a small change from here. Measured at $\{3,\cdot,20\}$ with
  $M = 0..12$: a single call costs 0.32 s, thirteen separate calls 1.55 s, and the row in one
  pass over the partitions 0.35 s — **4.4x**, and barely more than one value.
- **The identity to test the sampler against**: the weights must sum to
  $n! \cdot U(p,M,n)$ exactly. A wrong cycle-type weight still emits plausible-looking
  complexes, so verify that identity, not just the outputs.
- **`RandomUniformUnlabeledPureSimplicialComplex`** (currently at `PureComplexes.wl:3716` and
  `:3778`) does rejection with acceptance $|\mathrm{Aut}|/n!$, which collapses factorially —
  15 s per 100 samples at $\{3,4,6\}$. It is verified uniform over isomorphism classes, so keep
  it as the slow-but-correct oracle to differential-test the replacement against.
