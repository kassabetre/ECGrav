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
that $n!/z_\lambda$ is the number of permutations of cycle type $\lambda$. $\mathrm{Sym}(\Omega)$
is the symmetric group on a finite set $\Omega$, so $S_n = \mathrm{Sym}([n])$. Square brackets carry two unrelated meanings in what follows, both
standard: $[n]$ as above, and $[x^k]f$ — likewise $[z^k]f$ — for the coefficient of $x^k$ in a
series $f$. The three counts of
§1 are written $U(p,M,n)$ for the fully unlabeled one specified here
(`NumUnlabeledPureComplexes`), $\mathrm{VL}(p,M,n)$ for the vertex-labeled
(`NumVertexLabeledPureComplexes`) and $\mathrm{FL}(p,M,n)$ for the facet-labeled
(`NumFacetLabeledPureComplexes`, written $s_F$ in its own specification).

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
not $M!/|\mathrm{Aut}(S)|$, which is only correct when that action is faithful, and it often
is not. The kernel is the set of vertex permutations fixing every facet **setwise**, and it is
nontrivial already at $p = 2$: for $S = \{\{1,2\},\{3,4\}\}$ the kernel
$\langle(1\,2),(3\,4)\rangle$ has order 4, so $|\mathrm{Aut}(S)| = 8$ while $|H| = 2$ — and
$M!/|\mathrm{Aut}(S)| = 1/4$ is not even an integer. At $p = 2$ this forces a disconnected
complex: if $\sigma$ swaps the ends $u,v$ of an edge and $uw$ is another edge at $u$, then
$\sigma(uw) = \{v,\sigma(w)\} \neq \{u,w\}$, so $u$ and $v$ both have degree 1. From
$p = 3$ the kernel survives connectivity: $(1\,2)$ fixes both facets of
$\{1,2,3\},\{1,2,4\}$ setwise. The vertex-labeled contribution is $n!/|\mathrm{Aut}(S)|$, which *is*
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
$m$ odd. Then $z P' = P \cdot \sum_m \ell(m) z^m$, with $P = \prod_d (1+z^d)^{n_d}$ the product whose
coefficients are the $c(m)$, gives the recurrence. $\square$

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
derivation above appears as the header comment at line 1934.

| line | symbol | role |
| --- | --- | --- |
| 1991 | `NumULPCPolyMul[a,b,deg]` | product of two coefficient lists, truncated at `deg`, via `ListConvolve` |
| 1996 | `NumULPCPowPoly[r,c,deg]` | $(1+z^r)^c$ truncated, stocking only multiples of $r$ |
| 2005 | `NumULPCFixedSubsets[parts,e,p]` | $f(e)$ of §2.5; `parts` is `Tally[λ]`, i.e. `{length, multiplicity}` pairs |
| 2020 | `NumULPCFixSets[parts,p,M]` | $|\mathrm{Fix}(\sigma)|$ by the Newton recurrence of §2.6 |
| 2046 | `NumULPCA[p,M,n]` | $A(p,M,n)$: the cycle-type sum of §2.3. **Memoized** |
| 2060 | `NumULPCCount[p,M,n]` | $A(n) - A(n-1)$ |
| 2065 | `NumUnlabeledPureComplexes[p,M,n]` | guards, then `NumULPCCount` |
| 2081 | `NumUnlabeledPureComplexes[p,M]` | summed over $n$ |
| 2094 | catch-all | `::argerr` and `$Failed` |

Two implementation choices that are exact, not heuristic:

- **Only $d \le M$.** Justified in §2.4; it also caps which $f(e)$ are needed. Checked to
  agree with the unrestricted version on 240 parameter sets.
- **`Divisible[L, d]` over `Range[M]`, not `Divisors[L]`.** The lcm $L$ of the cycle lengths
  is bounded by Landau's function $g(n)$, which grows fast enough to have many more divisors
  than a typical $M$: $g(20) = 420$ (24 divisors), $g(30) = 4620$ (48). Testing $M$ candidates
  beats enumerating every divisor of $L$ and filtering.

`NumULPCA` is memoized for the session and released by the shared
``ECGrav`Private`NumPCClearCache[]``, alongside the vertex-labeled row cache and the sampler
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
(`NumFacetLabeledPureComplexes` **as shipped** has no such recurrence either — it is already a
cycle-type sum. $\mathrm{FL}$ itself does have one, the partition-state deletion recursion of
`FacetLabeledCount.md` §9.3, but over partition states rather than in $(M,n)$; §7.1 below asks
whether it transfers to $U$, and it does not. Its Stirling relation $B(p,M,n) = \sum_k S(M,k)\,\mathrm{FL}(p,k,n)$ — with $S(M,k)$ the
Stirling numbers of the second kind and $B$ the number of *covering* tableaux with repeated columns
allowed (`FacetLabeledCount.md` §3.5) — is an identity between two counting problems, used for
verification, not a means of computing $\mathrm{FL}$.)

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

### 7.1 Why the facet-labeled deletion recursion does not transfer

`FacetLabeledCount.md` §9.3 gives $\mathrm{FL}$ a Burnside-free recursion, linear in $M$, whose state
is an integer partition. Since it sidesteps cycle types entirely, the obvious question is whether it
transfers to $U$. It does not, and the obstruction is the one named above, now with a measurement.

That recursion rests on two facts about a configuration $C$ of $i$ facets: the orbits of one-facet
extensions depend only on the state, and the number of **forbidden** branches — extensions repeating
an existing facet — is exactly $i-1$, a constant. Unlabeled, $\mathrm{Aut}(C)$ may permute the
facets, so the forbidden count becomes the number of $\mathrm{Aut}(C)$-orbits on the facets of $C$,
written $\varphi(C)$ below, and it is neither $i$ nor a function of the partition. At $p = 3$,
$i = 3$, $n = 6$, two configurations with the same $\lambda = (2,2,2,1,1,1)$:

| $C$ | $\lvert\mathrm{Aut}(C)\rvert$ | $\varphi(C)$ | extension orbits at $t = 0,1,2$ |
| --- | --- | --- | --- |
| $\{123,\,124,\,356\}$ | $4$ | $3$ | $7,\ 8,\ 4$ |
| $\{123,\,145,\,246\}$ | $6$ | $1$ | $5,\ 4,\ 2$ |

The second has $\mathrm{Aut}$ transitive on its three facets, collapsing them to one orbit. The group
responsible is $H_C$, the image of $\mathrm{Aut}(C)$ in $\mathrm{Sym}\{F_1,\dots,F_M\}$ — the same
group that gives each class $M!/|H_C|$ facet-labeled lifts. Hence

$$\mathrm{FL}(p,M,n) \;=\; \sum_{C} \frac{M!}{|H_C|}, \qquad\qquad U(p,M,n) \;=\; \sum_{C} 1,$$

the sums running over isomorphism classes. An ordered facet-by-facet construction computes the left
sum natively; the weights are the entire difference between the two problems, not bookkeeping at the
margin.

**Averaging over $S_M$ instead is worse, not better.** $U$ is the set of $S_M$-orbits on the
$\mathrm{FL}$ facet-labeled classes, so Burnside gives

$$U(p,M,n) \;=\; \sum_{\mu \,\vdash\, M} \frac{|\mathrm{Fix}(\mu)|}{z_\mu}, \qquad |\mathrm{Fix}(1^M)| = \mathrm{FL}(p,M,n),$$

with $\mathrm{Fix}(\mu)$ the classes admitting a vertex permutation realising a facet permutation of
type $\mu$. Verified exactly at $(3,3,6)$, $(2,4,5)$ and $(3,4,6)$, where the sums return $3$, $4$
and $15$ and the identity terms are $\mathrm{FL} = 10$, $29$ and $154$. It needs only $P(M)$ terms
against §2.3's $P(n)$, which looks attractive when $n \gg M$ — but every $|\mathrm{Fix}(\tau)|$ with
$\tau \neq \mathrm{id}$ is itself a count of $S_n$-orbits, so the total is $P(M)\,P(n)$. This is the
$S_n \times S_M$ route that §2.1 discards, reached from the other direction.

**Nor is $\mathrm{FL}/M!$ a usable approximation.** It is exact when every class is asymmetric, and
the ratio $(\mathrm{FL}/M!)/U$ might be expected to approach $1$ with $n$. It does the opposite, at
$p = 3$, $M = 6$:

| $n$ | $6$ | $7$ | $9$ | $12$ | $15$ | $18$ |
| --- | --- | --- | --- | --- | --- | --- |
| $(\mathrm{FL}/M!)/U$ | $0.60$ | $0.70$ | $0.67$ | $0.45$ | $0.12$ | $0.0014$ |

It peaks below $0.7$ and decays to $1/M!$. Sparse complexes are **more** symmetric, not less: at
$n = pM$ every facet is disjoint from every other, $H_C = S_M$, and all $M!$ orderings coincide. The
usual "generic objects are asymmetric" intuition points the wrong way at large $n$.

### 7.2 The block-vector signature: invariant, but not complete

A natural attempt to enumerate classes directly. Given a tableau, group its rows by **length** into
blocks $b_1,\dots,b_s$, with $\ell(b_j)$ rows and $|b_j|$ boxes; let $v_{ij}$ be the number of
block-$j$ vertices lying in facet $i$, and take the signature to be $\lambda$ together with the
multiset $\{v_1,\dots,v_M\}$ of the resulting length-$s$ vectors.

**It is genuinely $S_M$-invariant.** Blocks are defined by row length, which relabeling facets
cannot change, so $\tau$ permutes the vectors without altering any of them. The signature therefore
descends to isomorphism classes.

**It is not complete.** The smallest failure is $(p,M,n) = (2,5,6)$, where $U = 9$ but only $8$
signatures occur. The two tableaux sharing one, at $\lambda = (2,2,2,2,1,1)$:

$$\bigl[[1,2],[1,3],[2,3],[4,5],[4],[5]\bigr] \qquad\text{and}\qquad \bigl[[1,2],[1,3],[2,4],[3,5],[4],[5]\bigr]$$

Four length-$2$ rows and two length-$1$ rows in each, and both give
$\{(2,0),(2,0),(2,0),(1,1),(1,1)\}$. At $p = 2$ these are graphs: the first is a **triangle** on
$F_1F_2F_3$ plus a path, the second a **path** on all six vertices. One carries a cycle and the other
does not. The signature records how many vertices of each block a facet takes and never which, and
that is exactly the distinction it cannot see.

| $(p,M,n)$ | $(2,5,5)$ | $(2,5,6)$ | $(3,4,6)$ | $(3,5,6)$ |
| --- | --- | --- | --- | --- |
| $U$ | $5$ | $9$ | $15$ | $37$ |
| signatures | $5$ | $8$ | $13$ | $27$ |

Every case with $n \le 5$ in $p \le 3$, $M \le 5$ is complete, which is what makes small hand-checks
agree; the gap opens at $n = 6$ and widens.

**Which signatures are realisable is also not settled by local conditions.** The natural necessary
conditions — $\sum_j v_{ij} = p$, $\sum_i v_{ij} = |b_j|$, and $v_{ij} \le \min(p, \ell(b_j))$ — are
sound, in that no realisable signature is ever rejected by them, but they are far from sufficient:
at $(3,5,6)$ they admit $58$ candidates against $27$ realisable. Two refinements were tested. A
vector that is **all-or-nothing** ($v_{ij} \in \{0, \ell(b_j)\}$ for every $j$) determines its facet
exactly and so may not repeat, since the facets are distinct; imposing that cuts $58$ to $31$ and
does most of the available work. Per-block Gale–Ryser realisability of the bipartite degree sequence
adds **nothing** beyond it at any tested point, the block column sums being already forced. A
residue survives — $4$ of the $31$ — failing for reasons that couple blocks to one another rather
than constraining any one of them.

**And completing that characterisation would still not give the count.** The realisable signatures
are by definition the *image* of the map from classes, so

$$\#\{\text{realisable signatures}\} \;\le\; U(p,M,n),$$

with equality exactly when the signature is complete. Perfecting the realisability test yields a
lower bound that is tight only in the range where the answer is already reachable by enumeration.

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

### 8.1 A second oracle: the joint canonical form

Write $X$ for the set of ordered $M$-tuples of pairwise distinct $p$-subsets of $[n]$ covering
$[n]$, on which $S_n$ acts by relabelling vertices and $S_M$ by permuting the facet labels. The
objects of §1 are $X/(S_n \times S_M)$ — tuples modulo both — and that double quotient can be taken
in either order. The brute-force oracle above takes it as $(X/S_M)/S_n$: it forms *sets* of facets,
which has already forgotten the
facet order, then groups those under $S_n$. The oracle here takes the other order,
$(X/S_n)/S_M$: start from the **facet-labeled** classes, which are $X/S_n$ and are exactly the
incidence tableaux of `FacetLabeledCount.md` §2, then quotient those by $S_M$.

**The construction.** A facet-labeled class is a multiset of rows, each row the set of facets
containing one vertex. A facet relabelling $\tau \in S_M$ acts *inside* the rows,
$R \mapsto \tau(R) = \{\tau(x) : x \in R\}$. Order the rows of a tableau by decreasing size, then
lexicographically within each size block — ties can only arise between identical rows, so this is a
total order — and define

$$\mathrm{canon}(S) \;=\; \min_{\tau \in S_M}\ \bigl(\text{rows of } \tau\!\cdot\!S,\ \text{size-then-lex sorted}\bigr),$$

the minimum being lexicographic on the resulting tableaux. Ranging over the whole orbit is what
makes the result independent of the representative. Then

$$U(p,M,n) \;=\; \#\{\,\mathrm{canon}(S) \;:\; S \text{ a facet-labeled class}\,\}.$$

**What it adds over leg 1.** Not a new enumeration — it is still bounded by the number of objects.
What it removes is *arithmetic*: leg 1 recovers the three counts through $n!/|\mathrm{Stab}|$ and
$M!/|H|$ and has to guard that orbit $\times$ stabiliser $= n!$, whereas this computes no
automorphism group, no stabiliser and no index. It counts distinct values of a canonical form. The
two failure modes are disjoint, and it exercises $\mathrm{FL}$ and $U$ against each other through the
tableau picture rather than through orbit–stabiliser bookkeeping.

**Worked example.** At $p=2$, $M=4$, $n=5$, take the tableau
$\{1,2\},\{1,3\},\{2,4\},\{3\},\{4\}$ — facets $\{v_1v_2, v_1v_3, v_2v_4, v_3v_5\}$, the path
$v_4 v_2 v_1 v_3 v_5$. Applying $\tau = (1{\to}2,2{\to}3,3{\to}4,4{\to}1)$ sends the rows to
$\{2,3\},\{2,4\},\{1,3\},\{4\},\{1\}$: the same path with its edges renamed. All $4! = 24$
relabellings collapse to $24/2 = 12$ distinct facet-labeled classes — the $2$ being $|H|$, the
image of the path's flip in $\mathrm{Sym}(S)$, as in §1 — and all 12 share the canonical form
$\{1,2\},\{1,3\},\{2,4\},\{3\},\{4\}$. Grouping all $29$ facet-labeled classes at $(2,4,5)$ this way
leaves exactly $4$, which is $U(2,4,5)$.

> **Sorting the rows and then discarding the facet labels does not work, and the wrong answer is
> plausible.** It is tempting to canonicalise the rows once and then read off the *set* of column
> occupation vectors, on the grounds that forgetting which column is which is what unlabelling the
> facets means. But the row order is computed lexicographically **from the facet numbers**, so it is
> not $\tau$-equivariant: relabelling the facets reorders the rows and permutes every column
> vector's coordinates. The result is a quotient by no group at all, and it lands strictly between
> the two counts — at $(2,4,5)$ it gives $10$ against $U = 4$ and $\mathrm{FL} = 29$, at $(3,4,5)$ it gives
> $11$ against $5$ and $43$. It both merges genuinely distinct facet-labeled classes and splits a
> single unlabeled one: at $(2,4,5)$ the classes
> $\{1,2\},\{1,3\},\{2,3\},\{4\},\{4\}$ and $\{1,2\},\{1,4\},\{2,4\},\{3\},\{3\}$ share a column set,
> while one unlabeled class with 12 facet-labelings produces 5 different ones. The minimisation over
> $\tau$ is not a convenience; it is what makes the construction a quotient.

**Verified** on 17 parameter sets across $p \in \{1,2,3,4\}$, $M \le 6$, $n \le 7$ — zero
mismatches, on values to $U(4,4,7) = 29$ from $\mathrm{FL} = 342$, and $U(2,6,6) = 15$ from $\mathrm{FL} = 3850$.

**Cost, and the range it is usable in.** The work is (number of facet-labeled classes) $\times\ M!$,
so it is enumeration-bound twice over: $0.013$ s at $(3,4,5)$, $0.115$ s at $(4,4,7)$, $0.405$ s at
$(2,5,6)$, but $32$ s at $(2,6,6)$ where $\mathrm{FL} = 3850$ and $M! = 720$. The $M!$ factor can be cut by
refining the facets on cheap invariants before branching — the standard canonical-augmentation
move — but the $\mathrm{FL}$ factor cannot: this is an oracle for small parameters, not a counting method.
§6's Burnside sum returns a 105-digit $U(5,50,16)$ from 231 cycle types in $0.16$ s, and no
canonical form competes with that, because Burnside never touches an object.

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
