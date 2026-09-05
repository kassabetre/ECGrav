# Counting facet-labeled pure complexes

> **What this is.** A specification of `NumFacetLabeledPureComplexes`, matching the
> implementation in `Kernel/Subpackages/PureComplexes.wl`. Companion to `UnlabeledCount.md`,
> which does the same for the fully unlabeled count, and to `FacetLabeledSampler.md`, which
> samples the objects counted here. It records the incidence-tableau transpose the whole count
> rests on, the Burnside argument over row slots, the single substitution that imposes distinct
> facets, the conventions, the cost, and what it was verified against.

Notation: $p$ is the **purity** (every facet has $p$ vertices), $M$ the **facet order** (number
of facets), $n$ the **vertex count**. The argument order is $(p, M, n)$ throughout, matching
`NumVertexLabeledPureComplexes` and `NumUnlabeledPureComplexes`. $[n]$ is $\{1,\dots,n\}$, $S_n$
the symmetric group, $\lambda \vdash n$ an integer partition, $m_k$ the number of parts of
$\lambda$ equal to $k$, and $z_\lambda = \prod_k k^{m_k} m_k!$, so $n!/z_\lambda$ permutations
have cycle type $\lambda$. Write $x^{(r)} = x(x-1)\cdots(x-r+1)$ for the falling factorial.

The count itself is written $s_F(p,M,n)$, as in `FacetLabeledSampler.md`; the letter $F$ is
reserved for facets.

---

## 1. The object being counted

A **facet-labeled pure complex** is an ordered $M$-tuple of pairwise **distinct** $p$-subsets of
$[n]$ covering $[n]$, taken **up to relabelling the vertices**. The facets carry labels
$1,\dots,M$; the vertices do not. Formally $s_F(p,M,n) = |X/S_n|$ with

$$X \;=\; \{(F_1,\dots,F_M) : F_i \subseteq [n],\ |F_i| = p,\ F_i \neq F_j,\ \textstyle\bigcup_i F_i = [n]\}.$$

Three counts sit on the same underlying objects, differing only in what carries a label:

| function | vertices | facets | counts |
| --- | --- | --- | --- |
| `NumVertexLabeledPureComplexes` | labelled | unlabelled | the sets $S$ themselves |
| `NumFacetLabeledPureComplexes` | unlabelled | labelled | $S_n$-orbits of ordered $M$-tuples |
| `NumUnlabeledPureComplexes` | unlabelled | unlabelled | $S_n$-orbits of the sets $S$ |

**Only the vertex-labeled count has a recurrence.** It is computed by a genuine two-index
deletion recursion, because with labelled vertices "how many complexes put $k$ new vertices in
the next facet" is a plain product of binomials. Nothing of the kind survives quotienting by
$S_n$ — see §9.

### 1.1 How $s_F$ decomposes over isomorphism classes

Per isomorphism class the facet-labeled contribution is $M!/|H|$, where $H$ is the **image** of
$\mathrm{Aut}(S)$ in $\mathrm{Sym}(S)$ — not $M!/|\mathrm{Aut}(S)|$, which is correct only when
that action is faithful, and it often is not. The kernel is the set of vertex permutations fixing
every facet **setwise**, and it is nontrivial already at $p = 2$: for
$S = \{\{1,2\},\{3,4\}\}$ the kernel $\langle(1\,2),(3\,4)\rangle$ has order 4, so
$|\mathrm{Aut}(S)| = 8$ while $|H| = 2$ — and $M!/|\mathrm{Aut}(S)| = 1/4$ is not even an
integer. At $p = 2$ this forces a disconnected complex: if $\sigma$ swaps the ends $u,v$ of an
edge and $uw$ is another edge at $u$, then $\sigma(uw) = \{v,\sigma(w)\} \neq \{u,w\}$, so
$u$ and $v$ both have degree 1. From $p = 3$ the kernel survives connectivity: $(1\,2)$ fixes both
facets of $\{1,2,3\},\{1,2,4\}$ setwise. Hence

$$s_F(p,M,n) \;=\; \sum_{\text{classes}} \frac{M!}{|H|},\qquad
\mathrm{VL}(p,M,n) \;=\; \sum_{\text{classes}} \frac{n!}{|\mathrm{Aut}(S)|},$$

the second with the plain stabiliser. This is a useful independent handle on small cases and was
checked directly (§10). It also explains why $s_F$ can **exceed** the vertex-labeled count: at
$(2,6,4)$ the only complex is all six edges of $K_4$, giving $\mathrm{VL}=1$, $s_F=30$, $U=1$.

---

## 2. The incidence tableau

Transposing the incidence matrix gives the picture the entire derivation is stated in, and the
one `FacetLabeledSampler.md` §1 and §11 refer back to.

A complex becomes a **multiset of $n$ nonempty subsets of $[M]$** — one row per vertex, listing
the facets that vertex lies in — in which **every label occurs in exactly $p$ rows**. The rows
form a multiset because the vertices are unlabelled. "The $M$ facets are pairwise distinct"
becomes "the $M$ columns are pairwise distinct", the **separating** condition. Covering becomes
"no row is empty".

Two properties make this the right representation:

- It is a **complete invariant**: two facet-labeled complexes are isomorphic iff their row
  multisets agree. This is what the test suite canonicalises with.
- It moves the symmetry group onto a set of $n$ slots that can be Burnside-averaged, whereas the
  original description has $S_n$ acting on subsets of $[n]$.

Counting tableaux directly over $n$ rows is hopeless; everything below is the Burnside average.

---

## 3. Derivation

### 3.1 Two conditions, four families, and the padding bijection

Two **independent** binary conditions are in play, and keeping them apart is the whole of this
subsection. On the columns: **separating** (pairwise distinct) or arbitrary. On the rows:
**covering** (no empty row) or empty rows allowed — call the latter **padded**. That gives four
families, and the derivation touches all four:

| | separating columns | arbitrary columns |
| --- | --- | --- |
| **covering rows** | $s_F(p,M,n)$ — the object being counted | $B(p,M,n)$ — §3.5, §10 |
| **padded rows** | $A(p,M,n)$ — what (3.1) differences | $\tilde B(p,M,n)$ — what the Burnside sum produces |

Formally, $A(p,M,n)$ is the number of multisets of exactly $n$ subsets of $[M]$, **empty subsets
allowed**, in which every label lies in exactly $p$ of them and **no two labels occupy the same
set of rows**; $\tilde B(p,M,n)$ is the same count without that last clause.

**Padding lemma.** Adding empty rows does not touch the columns. An empty row contributes a $0$ to
every column, so it can neither merge two distinct columns nor separate two equal ones. Padding a
covering tableau on $k \le n$ rows with $n-k$ empty rows is therefore a bijection onto the padded
tableaux on $n$ rows having exactly $k$ nonempty rows — **within either column class, separately**.
Summing over $k$ and differencing,

$$A(n) \;=\; \sum_{k \le n} s_F(k), \qquad s_F(p,M,n) \;=\; A(p,M,n) - A(p,M,n-1), \tag{3.1}$$

and identically $\tilde B(n) = \sum_{k \le n} B(k)$ with $B(n) = \tilde B(n) - \tilde B(n-1)$.

> **The column classes must not be mixed, and the risk is real.** The Burnside average of
> §3.2–§3.3 naturally produces $\tilde B$, not $A$; $A$ is reached from it only by the
> substitution of §3.5. Differencing $\tilde B$ instead yields $B$ — the covering count *with
> repeated columns allowed* — which is not $s_F$: at $(2,3,3)$ it gives $4$ against $s_F = 1$, at
> $(3,3,4)$ it gives $4$ against $1$, and at $(3,3,6)$ it gives $13$ against $10$. It is $A$ that
> (3.1) differences, and the shipped `NumFLPCCount` accordingly passes the falling-factorial
> weight (§5).

This is the same padding device `UnlabeledCount.md` §2.2 uses, but **that derivation cannot go
wrong in this particular way and this one can**. There the objects are *sets* of $M$ distinct
$p$-subsets, so facet-distinctness is inherent in the object and its $A$ is "the same count with
covering not required" — one facet class throughout, nothing to mismatch. Here the facets are an
ordered tuple and distinctness has to be *imposed*, by the falling-factorial substitution of §3.5;
that substitution is what creates the second column class, and with it the opportunity to
difference the wrong family.

Doing covering by differencing rather than imposing it inside the sum is what keeps the summand a
plain product; the price is two cycle-type sums instead of one (§8 shows they collapse to one in
closed form).

### 3.2 Burnside over the $n$ row slots

The tableaux of §3.1 are **multisets** of rows — the vertices are unlabelled, so nothing
distinguishes one row slot from another. Burnside needs a set carrying a genuine group action, and
on the multisets themselves $S_n$ would act **trivially**: reordering the elements of a multiset
returns the same multiset, every orbit is a singleton, and the average would degenerate to the
tautology "the number of multisets is the number of multisets".

So put the labels back on first. Let

$$\tilde X \;=\; \bigl\{(R_1,\dots,R_n) \;:\; R_j \subseteq [M],\ \ \#\{\,j : i \in R_j\,\} = p \ \text{ for every } i \in [M]\bigr\}$$

be the set of **ordered** $n$-tuples of rows meeting the label condition — the same tableaux, but
with the $n$ row slots temporarily numbered $1,\dots,n$ — and let $S_n$ act by permuting those
positions,

$$\sigma\cdot(R_1,\dots,R_n) \;=\; (R_{\sigma^{-1}(1)},\dots,R_{\sigma^{-1}(n)}).$$

**This action is not trivial**, and $\tilde X$, not the set of multisets, is what Burnside is
applied to. Two tuples lie in the same orbit exactly when one is a rearrangement of the other,
which is exactly when they carry the same multiset of rows. So the orbits of $S_n$ on $\tilde X$
correspond one-to-one with the tableaux counted by $\tilde B$, the label condition being
$S_n$-invariant, and counting the latter means counting the former:

$$\tilde B(p,M,n) \;=\; \frac{1}{n!}\sum_{\sigma \in S_n} |\mathrm{Fix}(\sigma)| \;=\; \frac{1}{n!}\sum_{\lambda \vdash n} \frac{n!}{z_\lambda}\,|\mathrm{Fix}(\lambda)|. \tag{3.2}$$

A tuple in $\tilde X$ fixed by $\sigma$ is **constant on each cycle** of $\sigma$, so choosing one
amounts to choosing a subset of $[M]$ per cycle — which is what §3.3 counts.

Labelling and then averaging the labelling away is the same move the object itself is built on:
§1 defines $s_F$ as $|X/S_n|$, an orbit count on labelled tuples, for exactly the same reason.
Note also the left-hand side: the Burnside average imposes no condition on the columns, so what
(3.2)–(3.3) compute is $\tilde B$. The separating family $A$ that (3.1)
needs comes from it by the single substitution of §3.5, and everything between here and there is
stated for $\tilde B$.

### 3.3 The labels decouple

Write $R_c \subseteq [M]$ for the subset attached to cycle $c$. Label $i$ occupies exactly the
rows of those cycles $c$ with $i \in R_c$, so the condition "label $i$ lies in exactly $p$ rows"
reads

$$\sum_{c\,:\,i \in R_c} |c| \;=\; p,$$

**separately for each label**. The labels therefore impose no joint constraint, and

$$|\mathrm{Fix}(\lambda)| \;=\; N(\lambda,p)^M, \qquad N(\lambda,p) \;=\; [x^p] \prod_{k \ge 1} (1 + x^k)^{m_k}, \tag{3.3}$$

$N(\lambda,p)$ being the number of sub-multisets of the cycle lengths summing to $p$ — each cycle
of length $k$ contributes a factor $(1+x^k)$, taken or not. This is a truncated cycle-index-style
generating function, the same one the sampler calls $N(\nu)$.

**The decoupling is the step that makes the count tractable**, and it is worth saying why it is
not obvious: the rows are shared between labels, so a joint condition would have been the
default. It survives only because the constraint is stated per label and the $R_c$ are recovered
from the per-label choices without interference.

### 3.4 Long cycles enter in bulk — they are not excluded

Since every part is positive, a cycle of length $k > p$ lies in no sub-multiset summing to $p$.
Hence $N(\lambda,p)$ depends only on $m_1,\dots,m_p$.

**It does not follow that such cycle types contribute nothing.** They contribute whatever the
short cycles give them; the long cycles simply go unused, and their rows come out empty — which
is legitimate here precisely because §3.1 allowed empty rows. So (3.2) groups the sum as: choose
$m_1,\dots,m_p$; let $s = n - \sum_{k \le p} k\,m_k$ be what is left; arrange the remainder in
$g(s)$ ways, $g(s)$ being the number of permutations of $[s]$ with **every** cycle longer than
$p$. The number of permutations with prescribed short multiplicities is then

$$\frac{n!}{s!\ \prod_{k \le p} k^{m_k} m_k!}\; g(s),$$

giving the sum that `NumFLPCTSum` walks. $g$ is tabulated by `NumFLPCLongCycleTable` off the
recurrence obtained by reading the cycle through the point $1$: it has some length $k>p$, laid
down in $(s-1)!/(s-k)!$ ways, leaving the same problem on $s-k$ points. (Equivalently
$\sum_s g(s) z^s/s! = \exp\bigl(\sum_{k>p} z^k/k\bigr)$; the two agree through $s=8$.)

> **This is a real difference from the sampler, and the direction is easy to get backwards.**
> `FacetLabeledSampler.md` §4 restricts to `IntegerPartitions[n, All, Range[p]]` because *there*
> the tuples must cover, so a long cycle forces $|\mathrm{Fix}(\sigma)| = 0$. Here covering is
> imposed by the differencing (3.1) instead, so every cycle type participates. Dropping the
> long-cycle types from both terms of (3.1) is **not** a valid simplification: it disagrees with
> the true count on 80 of the in-range triples with $p \le 4$, $M \le 7$, $n \le 16$. The smallest
> is $(2,3,6)$, where it returns $\tfrac{2}{3}$ against the true $1$ — **not even an integer,
> and in range**; further out it gives $1$ instead of $6$ at $(2,4,7)$ and $165$ instead of
> $175$ at $(2,5,7)$. §8 says exactly how much of the factor it drops.

### 3.5 Separating is one substitution

Everything so far counts tableaux with **arbitrary** columns: (3.2)–(3.3) with the weight $N^M$
give $\tilde B(p,M,n)$, the padded arbitrary-column family of §3.1. The padded **separating**
family $A$ — the one (3.1) needs — is obtained by replacing that weight with the falling
factorial:

$$N^M \;\longrightarrow\; N^{(M)} \;=\; N(N-1)\cdots(N-M+1). \tag{3.5}$$

The reason is the standard identity $x^M = \sum_k S(M,k)\, x^{(k)}$ [2]: merging equal columns of
a tableau partitions the $M$ labels into $k$ blocks and leaves a separating tableau with $k$
columns. Merging columns neither creates nor destroys an empty row, so the correspondence holds
inside each **row** class separately:

$$B(p,M,n) \;=\; \sum_{k} S(M,k)\, s_F(p,k,n), \qquad
\tilde B(p,M,n) \;=\; \sum_{k} S(M,k)\, A(p,k,n). \tag{3.6}$$

Both hold, and the **row condition must match across the equation**. The first is the form the
suite checks (§10, leg 3), whose brute force enumerates covering tableaux with repeated columns
allowed. Pairing $\tilde B$ with $s_F$ is the mismatch to avoid: at $(2,3,5)$, $\tilde B = 15$
while $\sum_k S(3,k)\,s_F(2,k,5) = 3$.

Substituting (3.5) cycle type by cycle type **is** the Stirling inversion of (3.6) — and unlike
performing that inversion afterwards, it has **no cancellation**: every summand stays
non-negative, because $N^{(M)} \ge 0$ for integer $N \ge 0$ and $N^{(M)} = 0$ exactly when
$N < M$, which is the honest statement that a block collection offering fewer than $M$ admissible
facets can supply no separating tableau at all.

`NumFLPCTSum` takes the weight as a parameter `wf` for this reason; passing `#^M &` returns
$\tilde B$ and passing `FactorialPower[#, M] &` returns $A$. The shipped code always
passes the latter, but the hook is what makes (3.6) directly testable (§10).

---

## 4. Worked example: $s_F(2,3,5) = 3$

Three labelled edges on five unlabelled vertices. By (3.1) this needs $A(5)$ and $A(4)$ — the
**padded separating** family, so the weight column below is the falling factorial $N^{(3)}$ and
not $N^3$.

$A(2,3,5)$, over all seven cycle types of $S_5$:

| $\lambda$ | $z_\lambda$ | #perms | $N(\lambda,2)$ | $N^{(3)}$ | weighted |
| --- | --- | --- | --- | --- | --- |
| $\{5\}$ | 5 | 24 | 0 | 0 | 0 |
| $\{4,1\}$ | 4 | 30 | 0 | 0 | 0 |
| $\{3,2\}$ | 6 | 20 | 1 | 0 | 0 |
| $\{3,1,1\}$ | 6 | 20 | 1 | 0 | 0 |
| $\{2,2,1\}$ | 8 | 15 | 2 | 0 | 0 |
| $\{2,1,1,1\}$ | 12 | 10 | 4 | 24 | 240 |
| $\{1^5\}$ | 120 | 1 | 10 | 720 | 720 |

$T(5) = 960$, so $A(5) = 960/120 = 8$.

$A(2,3,4)$, over the five cycle types of $S_4$:

| $\lambda$ | $z_\lambda$ | #perms | $N(\lambda,2)$ | $N^{(3)}$ | weighted |
| --- | --- | --- | --- | --- | --- |
| $\{4\}$ | 4 | 6 | 0 | 0 | 0 |
| $\{3,1\}$ | 3 | 8 | 0 | 0 | 0 |
| $\{2,2\}$ | 8 | 3 | 2 | 0 | 0 |
| $\{2,1,1\}$ | 4 | 6 | 2 | 0 | 0 |
| $\{1^4\}$ | 24 | 1 | 6 | 120 | 120 |

$T(4) = 120$, so $A(4) = 5$ and $s_F(2,3,5) = 8 - 5 = 3$.

**Hand-check.** Three edges covering five vertices force degree sequence $(2,1,1,1,1)$, so the
only isomorphism class is $P_3 \sqcup K_2$. Its automorphism group has order $4$ (flip the path,
flip the edge), but the edge-flip acts trivially on the edge set, so the image $H$ has order $2$
and the class contributes $3!/2 = 3$ by §1.1. This is also a live check on §3.4: the rows
$\{3,2\}$, $\{3,1,1\}$ and $\{2,2,1\}$ carry $N > 0$ and are only killed by the falling factorial,
not by being long.

---

## 5. Implementation map

All in `Kernel/Subpackages/PureComplexes.wl`. Private helpers are named `NumFLPC*`; the
derivation above appears as the header comment at line 1801.

| line | symbol | role |
| --- | --- | --- |
| 1834 | `NumFLPCLongCycleTable[p,smax]` | $g(0..s_{\max})$, permutations with every cycle $> p$ (§3.4) |
| 1853 | `NumFLPCNCoeff[mvec,jvecs]` | $N(\lambda,p)$ of (3.3), expanded over the partitions of $p$ |
| 1859 | `NumFLPCTSum[p,n,wf,gg,fac,jvecs]` | the cycle-type sum: $n!\,A(n)$ with `wf` $= N^{(M)}$, or $n!\,\tilde B(n)$ with `wf` $= N^M$ (§3.1) |
| 1884 | `NumFLPCCount[p,M,n]` | $\bigl(T(n) - n\,T(n-1)\bigr)/n!$, the differencing (3.1) |
| 1897 | `NumFacetLabeledPureComplexes[p,M,n]` | guards, then `NumFLPCCount` |
| 1913 | `NumFacetLabeledPureComplexes[p,M]` | summed over the vertex count |
| 1926 | catch-all | `::argerr` and `$Failed` |

Three implementation details worth naming:

- **$N$ is expanded over the partitions of $p$, not built as a polynomial.** `jvecs` holds, for
  each $\nu \vdash p$, its part counts, and $N = \sum_{\nu \vdash p} \prod_k \binom{m_k}{\nu_k}$ —
  choosing $\nu_k$ cycles of length $k$. This costs $P(p)$ products of $p$ binomials, independent
  of $n$ and $M$. The sampler's `RandFLPCWeightCounts` instead builds the whole coefficient list
  by truncated convolution, because it needs every coefficient out to $x^p$, not just the last.
- **The differencing is done on $T = n!A$, not on $A$.** `NumFLPCCount` returns
  $(T(n) - n\,T(n-1))/n!$, keeping the arithmetic in integers throughout; $T(n-1) = (n-1)!A(n-1)$,
  so the factor $n$ rescales it. Nothing rational is ever formed.
- **Nothing here is memoized.** Unlike `NumULPCA`, which caches per $(p,M,n)$, the `NumFLPC*`
  helpers hold no `DownValues` between calls, and `NumPCClearCache[]` correctly does not list
  them. A repeated call repeats the work. At the measured costs (§7) that is a defensible
  trade — a call is milliseconds and the natural callers ask once — but it is a choice, not an
  oversight, and it differs from every other counter in the file.

---

## 6. Conventions and guards

```wolfram
Which[
  p < 0 || M < 0,            0,
  M == 0,                    If[n == 0, 1, 0],
  n < 0 || n < p || n > p*M, 0,
  Binomial[n, p] < M,        0,
  True,                      NumFLPCCount[p, M, n]]
```

- $M = 0$ is the empty complex: $1$ at $n = 0$, else $0$.
- Zero unless $p \le n \le pM$ — fewer than $p$ vertices carry no facet, and $M$ facets of $p$
  vertices cover at most $pM$.
- Zero when $\binom{n}{p} < M$: there are not enough distinct $p$-subsets to supply $M$ distinct
  facets. This one is specific to the labelled-and-distinct reading and has no analogue in a
  multiset count.

These match `NumVertexLabeledPureComplexes` and `NumUnlabeledPureComplexes` exactly, and the
degenerate-input test asserts agreement rather than restating values.

**The guards are load-bearing, not an optimisation** — the same point as `UnlabeledCount.md` §5,
for a different reason. Falling through on an out-of-range $n$ enters the cycle-type sum at that
$n$; $(3,4,100)$ would walk every multiplicity vector for $n = 100$. It is in the degenerate test
for exactly this reason.

**The two-argument form sums rather than telescoping.** It finds the least $n$ with
$\binom{n}{p} \ge M$ and adds $s_F(p,M,n)$ up to $n = pM$. Since (3.1) telescopes, the whole sum
is just $A(p,M,pM)$ — verified on $p \in \{2,3\}$, $M \in \{2,\dots,5\}$ — so the form currently
pays $O(pM)$ cycle-type sums where one would do. Measured: $0.038$ s at $(3,6)$, $0.098$ s at
$(3,8)$, $0.213$ s at $(3,10)$. Not a correctness issue, and not yet changed.

---

## 7. Cost

The sum is over **multiplicity vectors $(m_1,\dots,m_p)$ with $\sum_k k\,m_k \le n$** — that is,
partitions of at most $n$ into parts at most $p$ — not over all $P(n)$ partitions. Long cycles
are absorbed into $g(s)$ rather than enumerated. Per vector the work is one $N$ (that is $P(p)$
binomial products) and one falling factorial.

**$M$ is essentially free**, entering only as the length of a falling factorial. At
$(p,n) = (3,12)$:

| $M$ | 4 | 6 | 10 | 20 | 50 | 100 | 200 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| seconds | 0.0024 | 0.0023 | 0.0023 | 0.0023 | 0.0024 | 0.0026 | 0.0032 |
| digits in the answer | 1 | 5 | 15 | 38 | 106 | 214 | 395 |

A $50\times$ increase in $M$ costs $35\%$, while the answer grows from one digit to 395. What
does cost is $n$, and then $p$:

| $\{p,M,n\}$ | mult. vectors | seconds | digits |
| --- | --- | --- | --- |
| $\{2,6,8\}$ | 25 | 0.0005 | 4 |
| $\{3,4,7\}$ | 31 | 0.0007 | 3 |
| $\{3,6,12\}$ | 102 | 0.0023 | 5 |
| $\{3,8,18\}$ | 274 | 0.0062 | 7 |
| $\{4,6,15\}$ | 295 | 0.0084 | 7 |
| $\{3,10,22\}$ | 458 | 0.0105 | 10 |
| $\{4,8,20\}$ | 717 | 0.0203 | 11 |
| $\{3,12,30\}$ | 1041 | 0.0239 | 9 |
| $\{3,20,40\}$ | 2282 | 0.0528 | 30 |
| $\{4,10,30\}$ | 2724 | 0.0769 | 12 |
| $\{5,10,30\}$ | 5326 | 0.1897 | 20 |
| $\{6,8,30\}$ | 8547 | 0.4187 | 15 |
| $\{3,100,12\}$ | 102 | 0.0027 | 214 |
| $\{3,200,40\}$ | 2282 | 0.0662 | 751 |

Time tracks the vector count almost exactly — $\approx 23\,\mu$s each at $p \le 3$, $28\,\mu$s at
$p = 4$, drifting to $49\,\mu$s at $p = 6$ as $P(p)$ grows. Asymptotically the count is $\Theta(n^p)$, but in the
usable range it grows far more gently than that suggests: at $n = 18$ it is $100, 274, 515, 769,
996$ for $p = 2,\dots,6$.

**This is the sharpest contrast with `NumUnlabeledPureComplexes`**, whose sum is over all $P(n)$
cycle types and which therefore costs $27$ s at $\{2,20,40\}$ where this returns in
milliseconds. The reason is structural and is stated in `UnlabeledCount.md` §2.3: with the facets
unlabelled a long cycle *can* be covered, by a facet whose orbit walks around it, so no cycle type
can be summarised away.

---

## 8. What the differencing and the long-cycle table amount to

Not shipped; an observation made while writing this specification, verified on ten parameter sets
including out-of-range zeros. It is the cleanest statement of why §3.4 cannot be shortcut.

Let $h_u = \sum_{\lambda \vdash u,\ \text{parts} \le p} N(\lambda,p)^{(M)}/z_\lambda$ and
$H(z) = \sum_{u \ge 0} h_u z^u$ — the short-cycle part alone. With
$G(z) = \sum_s g(s) z^s/s! = \exp\bigl(\sum_{k>p} z^k/k\bigr)$ the long-cycle EGF, §3.4 says
$A(n) = [z^n]\,H(z)G(z)$ — note $h_u$ is built from $N^{(M)}$, so this is the padded separating
family of §3.1. Since $\exp\bigl(\sum_{k \ge 1} z^k/k\bigr) = 1/(1-z)$, the differencing
(3.1) contributes the factor $(1-z)$ and the two collapse:

$$s_F(p,M,n) \;=\; [z^n]\; H(z)\,\exp\!\Bigl(-\sum_{k=1}^{p} \frac{z^k}{k}\Bigr). \tag{8.1}$$

The naive "restrict to parts $\le p$" of §3.4 is exactly the substitution of $(1-z)$ for that
exponential. The two series agree through $z^p$ and **first differ at $z^{p+1}$** — for $p=2$,
$\exp(-z-z^2/2) = 1 - z + 0\,z^2 + z^3/3 - \cdots$ against $1 - z$ — which is precisely the
discrepancy §3.4 reports. It is not a small correction: the coefficients do not decay.

(8.1) also suggests a refinement, unmeasured: one pass over the restricted multiplicity vectors
building $H$ to degree $n$, convolved once with a fixed length-$(n{+}1)$ series, would replace the
two cycle-type sums of `NumFLPCCount` with one and give the whole row $n = 0..N$ at once.

---

## 9. What is *not* available: a recurrence

The vertex-labeled count has a genuine two-index deletion recurrence. **The facet-labeled and
unlabeled counts do not, and the reason is structural**: on isomorphism classes the deletion fibre
is the number of $\mathrm{Aut}(K)$-orbits on candidate facets, which varies from class to class,
so the deletion argument has nothing uniform to count. A *calibrated* ansatz search confirmed
this — the same search recovers the vertex-labeled recurrence exactly when pointed at it, so a
null result is informative rather than a failure of the search.

**Do not over-read the growth argument.** $\log T(n) = \Theta(n^2)$ rules out P-recursiveness in
$n$ alone, but the vertex-labeled count has the same growth *and* a perfectly good two-index
recurrence, so that argument excludes pure-$n$ recurrences only. It is not evidence against a
recurrence in $(M,n)$ or $(p,n)$. `UnlabeledCount.md` §7 states this at length; the same caveat
applies verbatim here.

One recurrence does exist a level down, though the facet-labeled case does not need it: Newton's
identity applied to $\log \prod_d (1+z^d)^{n_d}$ turns the coefficient extraction inside a single
cycle type into $m\,c(m) = \sum_k l(k)\,c(m-k)$. That is what carries the unlabeled count's
$O(M^2)$; here the per-type work is $P(p)$ binomial products and there is nothing to accelerate.

---

## 10. Verification

Five independent legs, in `Tests/PureComplexes.wlt` lines 224–293 unless noted.

1. **Against the source notebook.** Eleven separating counts from the author's
   `Facet-labeled-count.nb` (`sFDirect`), plus its 93-digit $s_F(3,50,10)$.
   `NumFacetLabeledPureComplexes-known-values`, `-large-M`.
2. **Against direct enumeration of the definition.** All $M$-tuples of distinct $p$-subsets
   covering $[n]$, canonicalised over the $n!$ relabellings, on five parameter sets.
   `NumFacetLabeledPureComplexes-brute-force`. This tests the *description*, not the derivation.
3. **Against the derivation itself, via (3.6).** $B(p,M,n) = \sum_k S(M,k)\,s_F(p,k,n)$ checked
   against independently brute-forced counts of all **covering** $(p,M,n)$ tableaux — separating or
   not — on four parameter sets. Covering on both sides is what makes the identity true; the brute
   force enforces it by drawing rows from `Rest[Subsets[Range[M]]]`, which omits the empty row. `NumFacetLabeledPureComplexes-stirling-identity`. This is the strongest
   leg: it exercises §3.5, the one step where a plausible-looking wrong answer is easiest to
   produce, rather than only the endpoints.
4. **Against the automorphism decomposition of §1.1.** $s_F = \sum_{\text{classes}} M!/|H|$ with
   $H$ the image of $\mathrm{Aut}$ in $\mathrm{Sym}(S)$, computed by explicit orbit enumeration,
   on $(2,3,4)$, $(2,3,5)$, $(2,4,5)$, $(3,3,5)$, $(2,4,6)$, $(3,4,6)$ — the last with 15 classes
   summing to 154. Not currently in the suite; worth adding, since it is the only check that
   touches §1.1 and it catches the $|\mathrm{Aut}|$-vs-$|H|$ error specifically.
5. **Structure and guards.** Degenerate inputs including the deliberately out-of-range
   $(3,4,100)$ of §6, the two-argument/sum agreement, and argument-error handling.

**On method.** Everything above is an exact integer identity. There is no statistical leg here
and none is wanted — contrast the samplers, where chi-square is available and is still the
weakest evidence (`FacetLabeledSampler.md` §9). A counter that is wrong is wrong by an integer.

---

## 11. Reference table

All three labellings, from one brute-force enumeration.

| $p$ | $M$ | $n$ | VL | **FL** | U |
| --- | --- | --- | --- | --- | --- |
| 2 | 3 | 4 | 16 | **4** | 2 |
| 2 | 3 | 5 | 30 | **3** | 1 |
| 2 | 3 | 6 | 15 | **1** | 1 |
| 2 | 4 | 4 | 15 | **15** | 2 |
| 2 | 4 | 5 | 135 | **29** | 4 |
| 2 | 4 | 6 | 330 | **19** | 3 |
| 2 | 5 | 4 | 6 | **30** | 1 |
| 2 | 5 | 5 | 222 | **222** | 5 |
| 2 | 5 | 6 | 1581 | **301** | 9 |
| 2 | 6 | 4 | 1 | **30** | 1 |
| 2 | 6 | 5 | 205 | **1230** | 5 |
| 2 | 6 | 6 | 3760 | **3850** | 15 |
| 3 | 3 | 5 | 100 | **7** | 3 |
| 3 | 3 | 6 | 480 | **10** | 3 |
| 3 | 4 | 5 | 205 | **43** | 5 |
| 3 | 4 | 6 | 3600 | **154** | 15 |
| 3 | 4 | 7 | 22820 | **207** | 17 |
| 3 | 5 | 6 | 13992 | **2472** | 37 |
| 4 | 3 | 6 | 395 | **8** | 4 |
| 4 | 4 | 7 | 42910 | **342** | 29 |
| 1 | 3 | 3 | 1 | **1** | 1 |

The rows $(2,5,4)$, $(2,6,4)$ and $(2,6,5)$ are the ones to keep in view: FL exceeds VL by a wide
margin there, which is §1.1 at work — few classes, each highly symmetric under vertex relabelling
and so contributing few vertex labellings, but each supplying up to $M!$ facet labellings.

---

## 12. Relation to the sampler

`FacetLabeledSampler.md` is not built on this count; it re-derives the same Burnside sum under a
covering constraint. The correspondence is worth stating because the two look interchangeable and
are not:

| | this count | the sampler |
| --- | --- | --- |
| covering | by differencing (3.1) | imposed inside the draw, by inclusion–exclusion |
| cycle types summed | all $\lambda \vdash n$ | only parts $\le p$ |
| per-type quantity | $N(\lambda,p)^{(M)}$ | $C(\lambda,M,0,\lambda)$, a covering count |
| $N$ computed as | last coefficient only, over $\nu \vdash p$ | full list to $x^p$, by convolution |
| memoized | no | yes, via `NumPCClearCache[]` |

The one place they meet is the identity the sampler is tested by: its step-1 weights must sum to
$n! \cdot s_F(p,M,n)$ exactly, which is (3.2) read backwards. That identity is the sampler's
strongest correctness leg, and this counter is its right-hand side — so an error here would be
masked there, not caught. Legs 2 and 3 of §10 exist for that reason.

`RandomUniformFacetLabeledPureSimplicialComplex[{p,M}, k]` also draws its vertex count from
$\Pr[n] \propto s_F(p,M,n)$ directly, and its emptiness guard tests $s_F \neq 0$ rather than the
vertex-labeled count — the $\binom{n}{p} \ge M$ condition of §6 being the reason the two differ.

---

## 13. References

1. W. Burnside, *Theory of Groups of Finite Order*, 2nd ed., Cambridge University Press, 1911,
   §145. — The orbit-counting lemma of (3.2). Due to Cauchy and Frobenius, not Burnside; see
   P. M. Neumann, "A lemma that is not Burnside's", *Math. Sci.* **4** (1979), 133–141.
2. R. P. Stanley, *Enumerative Combinatorics*, Vol. 1, 2nd ed., Cambridge University Press, 2012,
   Ch. 1–2. — The falling factorial, the identity $x^M = \sum_k S(M,k)\,x^{(k)}$ of §3.5, and the
   exponential formula behind the long-cycle EGF of §3.4.
3. G. Pólya, "Kombinatorische Anzahlbestimmungen für Gruppen, Graphen und chemische
   Verbindungen", *Acta Math.* **68** (1937), 145–254; English translation in G. Pólya and
   R. C. Read, *Combinatorial Enumeration of Groups, Graphs, and Chemical Compounds*, Springer,
   1987. — Cycle-index enumeration; $N(\lambda,p) = [x^p]\prod_k(1+x^k)^{m_k}$ is a truncated
   instance. See also Stanley, Vol. 2, §7.24.
4. J. H. Redfield, "The theory of group-reduced distributions", *Amer. J. Math.* **49** (1927),
   433–455. — The earlier, independent statement of the same theory.
5. G.-C. Rota, "On the foundations of combinatorial theory I: Theory of Möbius functions",
   *Z. Wahrscheinlichkeitstheorie* **2** (1964), 340–368. — Möbius inversion over the partition
   lattice. The falling-factorial substitution of §3.5 is that inversion in closed form; see
   `FacetLabeledSampler.md` §11 for what it costs to re-derive it by brute force instead.
6. H. S. Wilf, *generatingfunctionology*, 3rd ed., A K Peters, 2006, Ch. 3. — The EGF manipulation
   of §8, and $\exp\bigl(\sum_k z^k/k\bigr) = 1/(1-z)$.

Companion specifications in this repository: `Theory.md` (the model), `FacetLabeledSampler.md`
(sampling the objects counted here), `UnlabeledCount.md` and `UnlabeledSampler.md` (the fully
unlabeled case, where the long-cycle collapse of §3.4 is unavailable and $n$ becomes the limiting
argument).
