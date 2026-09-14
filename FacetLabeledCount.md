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

> **Why the inverse, and why it costs nothing.** It is what makes this a *left* action:
> $(\sigma\cdot(\tau\cdot R))_j = R_{\tau^{-1}(\sigma^{-1}(j))} = R_{(\sigma\tau)^{-1}(j)}$, so
> $\sigma\cdot(\tau\cdot R) = (\sigma\tau)\cdot R$. Writing $R_{\sigma(j)}$ instead gives
> $\sigma*(\tau*R) = (\tau\sigma)*R$ — the factors reverse, which is an anti-homomorphism, and the
> sentence above would then be claiming something false of it. The inverse is also the convention
> that reads as a push-forward, "the row in slot $i$ moves to slot $\sigma(i)$", which is what makes
> this the *same* $\sigma$ that relabels vertices in §1 rather than its inverse.
>
> **Nothing downstream depends on the choice**, and not merely up to a relabelling of the sum: the
> two conventions have the *same fixed-point set* for every $\sigma$. Fixed under the inverse form
> means $R$ is constant on the cycles of $\sigma^{-1}$, fixed under the other means constant on the
> cycles of $\sigma$ — and $\sigma$ and $\sigma^{-1}$ have the same cycles, traversed the other way.
> Both are the condition §3.3 uses. Checked at $M=2$, $n=4$: $|\mathrm{Fix}(\sigma)|$ agrees for
> every one of the $24$ permutations individually, and both averages give the $35$ orbits that
> $\binom{4+4-1}{4}$ predicts.

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

**A fixed tuple is a choice of one subset per cycle.** $\sigma$ fixes $(R_1,\dots,R_n)$ exactly
when $R_{\sigma(j)} = R_j$, that is, when the tuple is constant on each cycle of $\sigma$. So the
$n$ row slots collapse to the $\ell(\lambda)$ cycles, and a fixed tuple is precisely a choice of
one subset $R_c \subseteq [M]$ for each cycle $c$.

Row $j$ then contains label $i$ iff $i \in R_{c(j)}$, and all $|c|$ rows of a cycle carry the same
subset, so label $i$ occupies $\sum_{c\,:\,i \in R_c}|c|$ rows in total. The condition "label $i$
lies in exactly $p$ rows" reads

$$\sum_{c\,:\,i \in R_c} |c| \;=\; p,$$

**separately for each label**, and

$$|\mathrm{Fix}(\lambda)| \;=\; N(\lambda,p)^M, \qquad N(\lambda,p) \;=\; [x^p] \prod_{k \ge 1} (1 + x^k)^{m_k}. \tag{3.3}$$

#### Why it is a power, and what would destroy it

The step that turns a shared constraint into a product is a **transpose**. Specifying the family
$\{R_c\}_c$ is the same data as specifying, for each label $i$, the set

$$S_i \;=\; \{\,c \;:\; i \in R_c\,\}$$

of cycles that contain it — the same cycle-by-label incidence read the other way round. Under that
reading the constraints are $\sum_{c \in S_i}|c| = p$, **one per label, with nothing linking
different labels**. Each of the $M$ labels therefore chooses independently from the same menu of
admissible $S$, and the count is that menu's size raised to the $M$.

This is worth pausing on, because it is fragile and it is not the default: the rows are shared
between labels, so a joint condition is what one would expect. It survives only because **nothing
constrains the $R_c$ themselves** — no cap on $|R_c|$, and, decisively, no requirement that a row
be nonempty.

> **Covering would destroy the decoupling.** A row is empty exactly when *no label* chose its
> cycle, which is a condition across all $M$ labels at once. Imposing covering here would couple
> them and the product would collapse. **This is the real reason §3.1 removes covering first**, and
> why the price of two cycle-type sums is worth paying: the differencing is not a convenience, it
> is what makes (3.3) a power at all.

#### What $N(\lambda,p)$ counts

For a single label, choose a sub-collection of the **cycles** whose lengths total $p$. Each cycle
is either taken or not, and one of length $k$ contributes $x^k$ when taken; with $m_k$ cycles of
length $k$ that is a factor $(1+x^k)^{m_k}$, and $[x^p]$ selects total length exactly $p$.
Equivalently, expanded over the partitions $\nu \vdash p$,

$$N(\lambda,p) \;=\; \sum_{\nu \,:\, \sum_k k\,\nu_k = p}\ \prod_k \binom{m_k}{\nu_k},$$

which is the form `NumFLPCNCoeff` evaluates. It is a truncated cycle-index-style generating
function, the same one the sampler calls $N(\nu)$.

> **It counts sub-collections of cycles, not sub-multisets of cycle lengths.** Cycles of equal
> length are distinguishable, which is exactly what the $\binom{m_k}{\nu_k}$ record. For
> $\lambda = \{2,1,1,1\}$ and $p = 2$ there are only two sub-multisets of *lengths* summing to 2,
> namely $\{2\}$ and $\{1,1\}$ — but $N = 4$, because the three 1-cycles offer $\binom{3}{2} = 3$
> distinct pairs. (Earlier revisions of this section described it the wrong way; the generating
> function and the code were always right.)

#### Examples

Taking $p = 2$, as in §4:

| $\lambda$ | generating function | $[x^2]$ | directly |
| --- | --- | --- | --- |
| $\{2,1,1,1\}$ | $(1+x)^3(1+x^2)$ | **4** | the 2-cycle, or two of the three 1-cycles: $1 + \binom{3}{2}$ |
| $\{1^5\}$ | $(1+x)^5$ | **10** | any two of the five 1-cycles: $\binom{5}{2}$ |
| $\{2,2,1\}$ | $(1+x)(1+x^2)^2$ | **2** | either 2-cycle; two 1-cycles are unavailable, there is one |
| $\{3,1,1\}$ | $(1+x)^2(1+x^3)$ | **1** | the two 1-cycles; the 3-cycle overshoots |
| $\{3,2\}$ | $(1+x^2)(1+x^3)$ | **1** | the 2-cycle |
| $\{4,1\}$ | $(1+x)(1+x^4)$ | **0** | nothing sums to 2 |

**A fixed tuple in full.** Take $n = 5$, $p = 2$, $M = 3$ and $\sigma = (1\,2)$, so the cycles are
$c_1 = \{1,2\}$ of length 2 and $c_2 = \{3\}$, $c_3 = \{4\}$, $c_4 = \{5\}$. The menu has the
$N = 4$ entries $\{c_1\}, \{c_2,c_3\}, \{c_2,c_4\}, \{c_3,c_4\}$. Let label 1 choose $\{c_1\}$,
label 2 choose $\{c_2,c_3\}$ and label 3 choose $\{c_3,c_4\}$. Transposing gives
$R_{c_1} = \{1\}$, $R_{c_2} = \{2\}$, $R_{c_3} = \{2,3\}$, $R_{c_4} = \{3\}$, so

$$(R_1,\dots,R_5) \;=\; \bigl(\{1\},\ \{1\},\ \{2\},\ \{2,3\},\ \{3\}\bigr),$$

constant on $c_1$ as required, with each label occupying exactly two rows. Read as a complex — row
$j$ lists the facets containing vertex $j$ — this is $F_1 = \{1,2\}$, $F_2 = \{3,4\}$,
$F_3 = \{4,5\}$: three edges on five vertices, fixed by the swap of vertices 1 and 2. All three
labels choose independently, so this cycle type contributes $4^3 = 64$ fixed tuples.

#### A consequence worth carrying to §3.5

Two columns of a fixed tuple coincide exactly when two labels chose the **same** menu entry, since
$S_i$ determines the set of rows label $i$ occupies. So, on $\sigma$-fixed tuples,

$$\text{separating} \iff \text{the } M \text{ labels pick pairwise \emph{distinct} menu entries},$$

which is $N(N-1)\cdots(N-M+1) = N^{(M)}$. That is why the substitution of §3.5 is **exact cycle
type by cycle type**, not a global identity that happens to come out right.

**Verified by brute force** — enumerating every $\sigma$-fixed tuple and comparing, for
$\lambda = \{2,1,1,1\}$ at $p=2, M=3$: $|\mathrm{Fix}| = 64 = 4^3$ and $24 = 4\cdot3\cdot2$
separating; for $\{1^5\}$: $1000 = 10^3$ and $720 = 10\cdot9\cdot8$; for $\{2,2,1\}$: $8 = 2^3$
and $0 = 2\cdot1\cdot0$. These are the $N$ and $N^{(3)}$ columns of §4's tables.

### 3.4 Long cycles enter in bulk — they are not excluded

Since every part is positive, a cycle of length $k > p$ lies in no sub-collection summing to $p$.
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

§3.3 already shows why this is the right weight and not merely a convenient one: on the fixed
tuples of a single cycle type, two columns coincide exactly when two labels chose the same
admissible $S$, so separating means the $M$ labels pick pairwise **distinct** menu entries, which
is $N^{(M)}$ on the nose. The substitution is therefore exact cycle type by cycle type. The
identity below is the same fact stated globally.

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
| 1859 | `NumFLPCTSum[p,n,wf,gg,fac,jvecs]` | the cycle-type sum: $n!\,A(n)$ with `wf` $= N^{(M)}$, or $n!\,\tilde B(n)$ with `wf` $= N^M$ (§3.1). Kept as the reference form; the shipped path now goes through the table below |
| — | `NumFLPCWeightTable[p,n]` | the $\{N,\text{weight}\}$ pairs of that sum with equal $N$ merged. **Memoized**, cleared by `NumPCClearCache[]` |
| — | `NumFLPCTFromTable[tab,M]` | $n!\,A(n)$ for one $M$, one pass over the table |
| — | `NumFLPCCount[p,M,n]` | $\bigl(T(n) - n\,T(n-1)\bigr)/n!$, the differencing (3.1) |
| — | `NumFacetLabeledPureComplexes[p,M,n]` | guards, then `NumFLPCCount` |
| — | `NumFacetLabeledPureComplexes[p,Mlist,n]` | the same for a list of facet orders, guards applied per entry |
| — | `NumFacetLabeledPureComplexes[p,M]` | summed over the vertex count |
| — | catch-all | `::argerr` and `$Failed` |

Three implementation details worth naming:

- **$N$ is expanded over the partitions of $p$, not built as a polynomial.** `jvecs` holds, for
  each $\nu \vdash p$, its part counts, and $N = \sum_{\nu \vdash p} \prod_k \binom{m_k}{\nu_k}$ —
  choosing $\nu_k$ cycles of length $k$. This costs $P(p)$ products of $p$ binomials, independent
  of $n$ and $M$. The sampler's `RandFLPCWeightCounts` instead builds the whole coefficient list
  by truncated convolution, because it needs every coefficient out to $x^p$, not just the last.
- **The differencing is done on $T = n!A$, not on $A$.** `NumFLPCCount` returns
  $(T(n) - n\,T(n-1))/n!$, keeping the arithmetic in integers throughout; $T(n-1) = (n-1)!A(n-1)$,
  so the factor $n$ rescales it. Nothing rational is ever formed.
- **The $(N,\text{weight})$ table is memoized per $(p,n)$** (since 1.14.1's successor; before that
  nothing here was cached and a repeated call repeated the work). Everything the cycle-type sum
  computes is independent of $M$ except the falling factorial, so the table serves every facet
  order at that vertex count, and equal $N$ are merged when it is built — different multiplicity
  vectors often reach the same $N$, which shrinks the sum by 1.5× at $p=2$, 2.3× at $p=3$ and
  2.9× at $p=4$. `NumPCClearCache[]` lists `NumFLPCWeightTable` and releases it.

  The natural
  caller does *not* ask once: Track A of `ExpansionRoadmap.md` sweeps $M$ at fixed $(p,n)$, and
  that is what this is for (§7).


**The clearing call is private, and getting the context wrong fails silently.**

```wolfram
ECGrav`Private`NumPCClearCache[]   (* works: returns the bytes reclaimed *)
ECGrav`NumPCClearCache[]           (* undefined: returns unevaluated, no message *)
```

An undefined symbol in Wolfram evaluates to itself, with no message, so the second form is a no-op
that looks exactly like success. The memo then survives and every timing taken after it is warm.
This corrupted a first pass at §9.3.8's Map 2: the shipped column read tens of microseconds instead
of milliseconds, and carried a spurious discontinuity precisely where the sweep first reached an
$n$ that no earlier row had already cached. The working call returns a byte count — **check for a
number, not for silence.**

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

**Sweeping $M$ at fixed $(p,n)$** is the access pattern that motivated the memoized table of §5,
and it is where the cost above is misleading: the figures are per call, and before the table each
call rebuilt $g$, the factorials and `jvecs` from scratch. Over $M = 1\ldots60$:

| $\{p,n\}$ | before | after | |
| --- | --- | --- | --- |
| $\{2,30\}$ | 0.252 s | 0.020 s | 12.6× |
| $\{3,30\}$ | 1.181 s | 0.059 s | 20× |
| $\{4,30\}$ | 3.939 s | 0.145 s | 27× |

A repeated scalar call at the same $(p,n)$ goes from 0.0397 s to 0.00068 s, a factor of 59. The
list form `NumFacetLabeledPureComplexes[p, Mlist, n]` is **not faster than the memoized scalar
form** — it is the same machinery — but it states the intent and does not depend on the memo
surviving a `NumPCClearCache[]`.
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

## 9. Recurrences: what is and is not available

The vertex-labeled count has a genuine two-index deletion recurrence. **The facet-labeled and
unlabeled counts do not, and the reason is structural**: on isomorphism classes the deletion fibre
is the number of $\mathrm{Aut}(K)$-orbits on candidate facets, which varies from class to class,
so the deletion argument has nothing uniform to count. A *calibrated* ansatz search confirmed
this — the same search recovers the vertex-labeled recurrence exactly when pointed at it, so a
null result is informative rather than a failure of the search.

**Read that as a statement about $(p,M,n)$ alone.** It says no recurrence closes in those three
indices; it does **not** say deletion is unavailable. §9.3 gives a deletion recursion that works,
by carrying a partition alongside them — a state rich enough to determine the very orbit count the
paragraph above calls non-uniform.

**Do not over-read the growth argument.** $\log T(n) = \Theta(n^2)$ rules out P-recursiveness in
$n$ alone, but the vertex-labeled count has the same growth *and* a perfectly good two-index
recurrence, so that argument excludes pure-$n$ recurrences only. It is not evidence against a
recurrence in $(M,n)$ or $(p,n)$. `UnlabeledCount.md` §7 states this at length; the same caveat
applies verbatim here.

One recurrence does exist a level down, though the facet-labeled case does not need it: Newton's
identity applied to $\log \prod_d (1+z^d)^{n_d}$ turns the coefficient extraction inside a single
cycle type into $m\,c(m) = \sum_k l(k)\,c(m-k)$. That is what carries the unlabeled count's
$O(M^2)$; here the per-type work is $P(p)$ binomial products and there is nothing to accelerate.

### 9.1 A constant-coefficient recurrence in $M$ — and why it is not an algorithm

(3.2) and (3.3) already give a closed form that is an **exponential sum in $M$**:

$$\tilde B(p,M,n) \;=\; \sum_{\lambda \vdash n} \frac{N(\lambda,p)^M}{z_\lambda},
\qquad B(p,M,n) \;=\; \tilde B(p,M,n) - \tilde B(p,M,n-1).$$

A finite sum of $M$-th powers satisfies a linear recurrence with **constant** coefficients whose
characteristic roots are the bases. Let $R$ be the distinct nonzero values of $N(\lambda,p)$ over
$\lambda \vdash n$ *and* $\lambda \vdash n-1$, and write $\prod_{N \in R}(x - N) = \sum_j c_j x^j$.
Then for every $M \ge 1$,

$$\sum_{j=0}^{|R|} c_j\, B(p,M+j,n) \;=\; 0.$$

| $p$ | $n$ | roots $R$ | order |
| --- | --- | --- | --- |
| 2 | 4 | $\{1,2,3,6\}$ | 4 |
| 2 | 6 | $\{1,2,3,4,7,10,15\}$ | 7 |
| 3 | 5 | $\{1,2,4,10\}$ | 4 |
| 3 | 6 | $\{1,2,4,8,10,20\}$ | 6 |

**It does not pass to $s_F$.** Feeding $s_F$ the same recurrence leaves large nonzero residuals —
at $p=2$, $n=6$: $\{24863,\,319407,\,1198102,\,767830\}$. The reason is structural: $s_F$ has the
same closed form with the falling factorial in place of the power,

$$s_F(p,M,n) \;=\; \sum_{\lambda \vdash n}\frac{N(\lambda,p)^{(M)}}{z_\lambda}
\;-\; \sum_{\lambda \vdash n-1}\frac{N(\lambda,p)^{(M)}}{z_\lambda},$$

and constant-coefficient recurrences are *exactly* the exponential sums. $N^{(M)}$ is not one. It
is P-recursive — $N^{(M+1)} = (N-M)\,N^{(M)}$, first order with a polynomial coefficient — so
$s_F$ is a sum of P-recursive sequences in $M$, which is the weaker statement.

**And it is not an algorithm.** The coefficients $c_j$ are elementary symmetric functions of the
roots, so knowing them means computing $N(\lambda,p)$ for every cycle type of $S_n$ and $S_{n-1}$
— exactly the work `NumFLPCTSum` already does, after which $\tilde B$ can be evaluated at any $M$
directly from the exponential sum. The order also grows with $n$. At fixed $n$ the sequence is
finite in any case: $s_F(p,M,n) = 0$ once $M > \binom{n}{p}$, and indeed
$\max_\lambda N(\lambda,p) = N(1^n,p) = \binom{n}{p}$, attained at the identity.

### 9.2 A dynamic programme on content vectors — Burnside-free

The recurrence of §9.1 is the wrong kind: it is read off the answer. A genuine recurrence for $B$
does exist, it never mentions cycle types, and it generalises the count to arbitrary content.

Let $B_m(n; r_1,\dots,r_m)$ be the number of incidence tableaux on labels $[m]$ with $n$ rows in
which label $i$ appears exactly $r_i$ times, so that $B(p,M,n) = B_M(n; p,\dots,p)$. Peel off the
last label. Every row containing $m$ has the form $R \cup \{m\}$ with $R \subseteq [m-1]$; let
$a_R$ count the rows of that shape, so $\sum_R a_R = r_m$, and write
$\rho_i(a) = \sum_{R \ni i} a_R$ for the occurrences of label $i$ consumed by those rows. The rows
*not* containing $m$ are then a tableau on $[m-1]$ with $n - r_m$ rows and the residual content,
and the split is a bijection because a row is sorted by whether it contains $m$. Hence

$$B_m(n; r_1,\dots,r_m) \;=\;
\sum_{\substack{a_R \ge 0\ (R \subseteq [m-1])\\ \sum_R a_R = r_m}}
B_{m-1}\bigl(n - r_m;\ r_1 - \rho_1(a),\ \dots,\ r_{m-1} - \rho_{m-1}(a)\bigr),$$

terms with a negative residual being zero, and specialising to $r_1 = \dots = r_M = p$ gives $B$.

> **The base case is the whole covering/padded distinction.** $B_0(n;()) = [\,n = 0\,]$ gives the
> **covering** family $B$: with no labels every row would be empty, so only $n = 0$ survives.
> Replacing it by $\tilde B_0(n;()) = 1$ for every $n \ge 0$ gives the **padded** family
> $\tilde B$, whose rows may be empty. Nothing else in the recurrence changes. Given that §3.1
> exists because those two axes were once conflated, this line is worth stating rather than leaving
> implicit. Equivalently, for one label $B_1(n;r) = [\,n = r\,]$.
>
> **The $n \ge 0$ is load-bearing and its omission is silent.** $\tilde B_0(n;()) = 1$ written for
> *all* $n$ lets a negative residual row count return $1$ instead of $0$, and the recursion then
> quietly stops enforcing the row count from below: what it computes is the saturated
> $\sum_{k} B(k)$ rather than $\sum_{k \le n} B(k)$. It is invisible from $s_F$, because the
> covering base $[\,n=0\,]$ is already zero at negative $n$ and the shipped route never touches
> $\tilde B$ — so the error surfaces only once $\tilde B$ is used in its own right, as the
> convolution below uses it. Caught by brute force against the definition: at $p=2$, $m=2$,
> $\tilde B$ must run $1, 2, 3$ over $n = 2,3,4$ and the unguarded form returns $3, 3, 3$.

**Verified** against the shipped counter, through (3.6), on 48 parameter sets across
$p \in \{2,3\}$, $M \le 4$, $n \le 6$ — zero mismatches, on values up to $B(3,4,6) = 221$.
Re-run 2026-09-10: reproduces exactly, and extends to 120 sets across $p \le 4$, $M \le 5$,
$n \le 8$, still zero, to $B(4,5,8) = 47986$.

**It yields a Burnside-free route to $s_F$.** Inverting the first form of (3.6),

$$s_F(p,M,n) \;=\; \sum_k s(M,k)\, B(p,k,n)$$

with $s(M,k)$ the signed Stirling numbers of the first kind. Composed with the recurrence above
this computes the facet-labeled count without ever forming a cycle type — verified to agree with
`NumFacetLabeledPureComplexes` on twelve parameter sets, up to $s_F(3,5,7) = 6561$; re-run
2026-09-10 over 80 sets ($p \in \{2,3\}$, $M \le 5$, $n \le 7$), zero mismatches, same maximum. As an independent
derivation it is a stronger check on the shipped algorithm than any of §10's legs. **It was the
first Burnside-free route recorded here and is no longer the only one** — §9.3 gives a second, on
a different state and much faster; where both apply they agree.

#### Cost, measured

| $p$ | $M$ | $n$ | DP + Stirling | shipped |
| --- | --- | --- | --- | --- |
| 2 | 3 | 5 | 0.004 s | 0.0004 s |
| 2 | 5 | 7 | 0.040 s | 0.0005 s |
| 3 | 5 | 5 | 0.098 s | 0.0005 s |
| 3 | 5 | 7 | 0.248 s | 0.0008 s |

As written it is 10–500× slower and diverging in $M$, where the shipped cost is nearly flat in
$M$. Two things are responsible: the state space is $(p+1)^m$ content vectors, and the transition
sums over $(a_R)$ across **all $2^{m-1}$ subsets** — compositions of $r_m$ into $2^{m-1}$ parts,
$\binom{r_m + 2^{m-1} - 1}{r_m}$ of them.

> **Both factors are singly exponential in $m$, not doubly** — an earlier revision of this section
> said doubly, and the correction is what makes the reductions below worth attempting rather than
> hopeless. Specialising to $r_1 = \dots = r_M = p$ fixes $r_m = p$, so the binomial is a degree-$p$
> polynomial in $2^{m-1}$:
>
> $$\binom{p + 2^{m-1} - 1}{p} \;=\; \Theta\!\left(\frac{2^{p(m-1)}}{p!}\right),$$
>
> singly exponential with base $2^p$. Measured over $m = 1\dots12$: successive ratios climb to
> $2^p$ from below — $3.99$, $7.95$, $15.81$ at $m = 10$ for $p = 2,3,4$ — and
> $\Delta \log_2(\text{count})$ converges to exactly $p$. The doubly exponential signature,
> constant $\Delta\log_2\log_2$, is absent: that quantity decays toward zero here, while the
> control $(p+1)^{2^{m-1}}$ holds it at exactly $1$.
>
> The reason is **sparsity**. With $\sum_R a_R = p$ at most $p$ of the $2^{m-1}$ parts are nonzero,
> so a transition picks a multiset of $p$ subsets out of the $2^{m-1}$ available rather than filling
> $2^{m-1}$ bins; the number of bins therefore enters polynomially, with degree $p$, however many
> there are. Drop that constraint — let each subset take a count in $0..p$ independently — and one
> does get the doubly exponential $(p+1)^{2^{m-1}}$, which is the likely origin of the error.
> Total naive cost is states × transitions $= O\bigl(((p+1)2^p)^m\bigr)$, base $32$ at $p = 3$:
> brutal, but singly exponential. **The base is what that bound gets right; the constant is not.**
> It charges every state the top level's out-degree, while most states sit at levels $m < M$ where
> $2^{m-1}$ is smaller and many have $r_m < p$. Instrumented at $(3,6,8)$: $2558$ states against
> the $(p+1)^M = 4096$ estimate, tight — but $300873$ transitions against the product's
> $1.07 \times 10^9$, loose by $3600\times$, the average out-degree being $118$ where a top-level
> state has $5984$. Successive ratios of the measured transition total do climb toward the
> predicted base — $6.2, 7.5, 9.2$ at $p = 2$ heading to $12$, and $11.4, 15.9, 21.5$ at $p = 3$
> heading to $32$ — so what is established here is the growth rate, not the constant. Even the
> least charitable reading — the unmemoised
> recursion tree, multiplying transition counts down the levels — gives $2^{\Theta(pm^2)}$,
> super-exponential and still not doubly.

#### Two reductions, and what is still open

- **Label symmetry.** $B_m(n;r)$ is symmetric in $r$ — permuting label names is a bijection on
  tableaux. Verified over all permutations of five content vectors, no asymmetric case. States
  therefore collapse to sorted $r$, i.e. $\binom{m+p}{p}$ of them, polynomial in $m$ at fixed $p$.
- **The subset sum is a convolution.** The summand depends on $a$ only through $\rho$, so grouping
  by the induced profile gives

  $$B_m(n;r) \;=\; \sum_{\rho} \tilde B_{m-1}(r_m;\rho)\; B_{m-1}\bigl(n - r_m;\ r_{<m} - \rho\bigr),$$

  with $\tilde B$ the padded variant — the same recurrence under the other base case.

  **Re-verified from scratch rather than trusted**: the original 42-parameter-set check predates
  the $n \ge 0$ correction to the padded base case above, and $\tilde B$ is exactly what that
  correction moves, so the figure it produced described a $\tilde B$ that was wrong. The
  replacement is 80 sets ($p \in \{2,3\}$, $M \le 5$, $n \le 7$) against the naive form in both
  families, plus 36 sets against brute-force enumeration of the definition — all four of
  $B$ and $\tilde B$, naive and convolved — and 80 sets of $s_F$ by both routes against the
  shipped counter. Zero mismatches throughout.

  This removes the $2^{m-1}$ **subsets** from the enumeration,
  replacing the $\binom{p + 2^{m-1} - 1}{p}$ compositions with the $(p+1)^{m-1}$ profiles $\rho$
  — but not the exponential in $m$: the base drops from $2^p$ to $p+1$, by a margin that itself
  grows with $m$.

  | $p$ | $m$ | compositions | profiles | ratio |
  | --- | --- | --- | --- | --- |
  | 2 | 6 | 528 | 243 | 2.2× |
  | 2 | 10 | 131328 | 19683 | 6.7× |
  | 3 | 6 | 5984 | 1024 | 5.8× |
  | 3 | 10 | 22500864 | 262144 | 86× |
  | 4 | 6 | 52360 | 3125 | 17× |
  | 4 | 10 | 2896986240 | 1953125 | 1483× |

  Measured on a fresh re-implementation of both forms, cold cache, verified against brute-force
  enumeration of the definition (so the naive column is not the same code as the table above and
  its figures differ from it):

  | $\{p,M,n\}$ | naive DP | convolution DP | shipped |
  | --- | --- | --- | --- |
  | $\{2,3,5\}$ | 0.0048 s | 0.0005 s | 0.0004 s |
  | $\{2,5,7\}$ | 0.064 s | 0.0061 s | 0.0006 s |
  | $\{3,5,5\}$ | 0.462 s | 0.021 s | 0.0005 s |
  | $\{3,5,7\}$ | 0.443 s | 0.022 s | 0.0008 s |
  | $\{3,6,8\}$ | 17.99 s | 0.140 s | 0.0011 s |

  10.5× at $M = 3$ rising to 129× at $M = 6$: the widening is the base change from $2^p$ to
  $p+1$, not a constant factor. It does not close the gap to the shipped counter, which is still
  129× ahead of the convolution at the last row.

**Open: whether the two compose.** The convolution is componentwise in $\rho$, so storing states
on sorted classes needs alignment bookkeeping — a symmetric-function product rather than a free
win. The corrected classes sharpen what is at stake: label symmetry alone makes the state count
$\binom{m+p}{p}$, polynomial, and the convolution alone makes the transition $(p+1)^{m-1}$, still
exponential — so the question is whether composing them reaches a polynomial total in $m$, or
whether the alignment cost puts the exponential back. Whether the combination beats the shipped
cost of $P(n)$ cycle types at $P(p)$ work each is unmeasured, and it would matter most at large
$n$, where the partition count is what hurts.


### 9.3 A deletion recursion on partition states — Burnside-free, and linear in $M$

§9.2's DP peels *labels* and carries a content vector. This one places *facets* in order and carries
a **partition**. It is the second Burnside-free route, it is much the faster of the two, and unlike
either it is linear in $M$. Contributed 2026-09-14; the derivation below is reconstructed and
verified, not taken on trust.

It also sharpens §9's opening. The obstruction there — the deletion fibre is the number of
$\mathrm{Aut}(K)$-orbits on candidate facets and varies from class to class — is a statement about
recurrences in $(p,M,n)$ **alone**. It is not an obstruction once the recursion carries a state
rich enough to determine that orbit count, and §9.3.2 says exactly which state does.

#### 9.3.1 The new-vertex profile

Order the facets and record how many vertices each one introduces that no earlier facet had:

$$n_i \;=\; \bigl|F_i \setminus (F_1 \cup \dots \cup F_{i-1})\bigr|.$$

Every isomorphism class has exactly one profile $(n_1,\dots,n_M)$, so the profiles **partition** the
classes and $s_F = \sum_c W(c)$ with $W(c)$ the number of classes carrying profile $c$. Three
conditions cut the candidates down:

1. $n_1 = p$ — facet 1 is entirely new;
2. $n_k \le p$ — a facet cannot introduce more than its own size;
3. $\binom{S_i}{p} \ge i$ with $S_i = \sum_{j \le i} n_j$ — after $i$ facets the pool must be able to
   supply $i$ *distinct* $p$-subsets.

Write $C(p,M,n)$ for the survivors. At $(p,M,n) = (3,4,5)$ there are three, and their weights sum
to the right answer:

| $c$ | $W(c)$ |
| --- | --- |
| $(3,1,0,1)$ | 6 |
| $(3,1,1,0)$ | 22 |
| $(3,2,0,0)$ | 15 |
| **total** | **43** $= s_F(3,4,5)$ |

Condition 3 is doing real work here and is exactly tight: the three profiles it rejects —
$(3,0,2,0)$, $(3,0,1,1)$, $(3,0,0,2)$ — all have $\binom{3}{p} = 1 < 2$ at $i=2$, and direct
enumeration confirms each supports **zero** classes. It rejects nothing that contributes.

#### 9.3.2 The state is a partition, and why that suffices

After $i$ facets, give each vertex its **type**, the set of facets containing it — the row of the
incidence tableau of §2. Two vertices of the same type are interchangeable, and
$\mathrm{Aut}(F_1,\dots,F_i)$ is precisely the group permuting vertices within types. So:

> **Orbits of $k$-subsets under $\mathrm{Aut}$ are sub-multisets of the type multiset.** Distinct
> orbits give non-isomorphic extensions, because an isomorphism of facet-*labelled* tuples must fix
> each $F_j$ setwise and is therefore an automorphism of the earlier structure.

That is what licenses throwing the type *identities* away and keeping only the multiset of type
**multiplicities** — a partition $\lambda$ of $S_i$. Choosing which old vertices the next facet
reuses is choosing a sub-multiset, and both the number of choices and the resulting partition depend
on $\lambda$ alone:

$$\bigl|MS_k(\lambda)\bigr| \;=\; [x^k] \prod_i \bigl(1 + x + \dots + x^{\lambda_i}\bigr).$$

| $\lambda$ | $\vert MS_2\vert$ | $\vert MS_3\vert$ |
| --- | --- | --- |
| $(3)$ | 1 | 1 |
| $(2,2,1)$ | 5 | 5 |
| $(2,1,1,1)$ | 7 | 7 |
| $(1^5)$ | 10 | 10 |
| $(3,2,1)$ | 5 | 6 |

Having chosen a sub-multiset, each partially-taken type **splits** into its taken and untaken halves
— they are no longer interchangeable, one being in the new facet — and the $n_i$ new vertices enter
as a fresh part. Fully-taken and untaken types do not split. That is the whole state transition.

#### 9.3.3 Distinctness, and why multiplicities survive it

Facet $i$ must differ from $F_1,\dots,F_{i-1}$. When $n_i > 0$ it contains a brand-new vertex and
cannot coincide with any of them, so nothing is needed. When $n_i = 0$ the correction is exactly
$i-1$: each $F_j$ is an available sub-multiset of size $p$, and $F_j \neq F_{j'}$ forces their
profiles apart, since $F_j$ takes all of every type containing $j$ and none of the rest.

That last sentence is also the reason the partition state survives a correction that seems to need
type identities:

> **The forbidden branches are all-or-nothing.** $F_j$ consumes entire types, so nothing splits and
> the resulting partition is **unchanged**. Multiplicities cannot say *which* $i-1$ branches are
> forbidden — and do not have to, because all $i-1$ land on the same state $\lambda$. The correction
> collapses to a single term, $-(i-1)\,W(\lambda;\,\cdot)$.

#### 9.3.4 Worked example: $W(3,2,0,0) = 15$

Profile $c = (3,2,0,0)$ at $p=3$, $M=4$. Partitions written with multiplicities descending.

**Stage 2** ($n_2 = 2$). State $(3)$. Facet 2 reuses $p - n_2 = 1$ old vertex; $|MS_1((3))| = 1$.
That type splits $1 + 2$, and the two new vertices join as a part: state $(2,2,1)$, five vertices
$= S_2$. In types: $\{1\}^2,\ \{2\}^2,\ \{1,2\}^1$.

**Stage 3** ($n_3 = 0$). Facet 3 is three old vertices; $|MS_3((2,2,1))| = 5$. Writing profiles in
the coordinates $(\{1\}, \{2\}, \{1,2\})$:

| profile | grows to | |
| --- | --- | --- |
| $(2,1,0)$ | $(2,1,1,1)$ | |
| $(1,2,0)$ | $(2,1,1,1)$ | |
| $(1,1,1)$ | $(1^5)$ | |
| $(2,0,1)$ | $(2,2,1)$ | **$= F_1$, forbidden** |
| $(0,2,1)$ | $(2,2,1)$ | **$= F_2$, forbidden** |

Both forbidden branches are all-or-nothing and leave the state at $(2,2,1)$, exactly as §9.3.3 says.
Three branches survive.

**Stage 4** ($n_4 = 0$). Count the leaves, less $M - 1 = 3$:

$$W(3,2,0,0) \;=\; \underbrace{(7-3)}_{(2,1,1,1)} + \underbrace{(7-3)}_{(2,1,1,1)} + \underbrace{(10-3)}_{(1^5)} \;=\; 4+4+7 \;=\; 15. $$

#### 9.3.5 Folding the profile away

Enumerating $C(p,M,n)$ and running a tree per profile duplicates work — profiles sharing a suffix
share every subtree below it. Memoising on $(\lambda, \text{remaining tail})$ recovers that, and is
worth 3.3–16×. But once the memo is keyed on the tail, the enumeration is doing nothing the memo is
not, so fold the choice of $n_i$ **into** the recursion and drop $C(p,M,n)$ entirely. With $j$
facets still to place,

$$V(\lambda, j) \;=\; \sum_{t=0}^{\min(p,\,r)}\ \sum_{a \,\in\, MS_{p-t}(\lambda)}
V\bigl(\mathrm{grow}(\lambda,a,t),\ j-1\bigr)\ -\ (M-j)\,V(\lambda,\ j-1), \tag{9.3}$$

the subtraction applying to the $t = 0$ term only, with base $V(\lambda, 0) = [\,|\lambda| = n\,]$
and

$$s_F(p,M,n) \;=\; V\bigl((p),\ M-1\bigr).$$

Here $t$ is the new-vertex count of the facet being placed, $r = n - |\lambda|$ is the budget left,
$\mathrm{grow}$ splits the partially-taken parts and appends $t$, and $M - j = i - 1$ is the number
of facets already down. Two things are worth pulling out:

- **The budget is not a free coordinate.** $r = n - |\lambda|$ is forced, so the state is
  $(\lambda, j)$ — two coordinates, not three.
- **Condition 3 becomes unnecessary.** It cannot over-subtract: the $M-j$ earlier facets are always
  distinct available orbits, so the bracket never goes negative, and infeasible profiles contribute
  zero on their own.

#### 9.3.6 Where it meets the Burnside derivation

Split the $t = 0$ term by whether the choice splits a type. The non-splitting choices are the
all-or-nothing ones — sub-collections of *parts* summing to $p$ — and that count is

$$A(\lambda) \;=\; [x^p]\prod_k (1+x^k)^{m_k} \;=\; N(\lambda, p),$$

the very quantity (3.3) is built on, here with $\lambda$ a partition of the vertex count rather than
a cycle type. Verified identical on every $\lambda$ tested. So (9.3)'s $t=0$ term reads

$$\sum_{a\ \text{splitting}} V\bigl(\mathrm{grow}(a),\, j-1\bigr) \;+\; \bigl(N(\lambda,p) - (M-j)\bigr)\,V(\lambda,\, j-1),$$

and $N(\lambda,p) - (M-j)$ is "how many all-or-nothing extensions are genuinely new facets". The two
routes meet at $N(\lambda,p)$ from opposite directions: §3.5 raises it to a falling factorial
$N^{(M)}$ in one stroke, and (9.3) decrements it one facet at a time. That is also the cleanest
statement of why the shipped counter is flat in $M$ and this is not — it never iterates over facets
at all.

#### 9.3.7 Cost

**The reachable $\lambda$ are partitions of integers $\le n$, a set that does not depend on $M$.**
Hence the state count grows linearly in $M$, and so does the method. At $(p,n) = (3,9)$:

| $M$ | 4 | 8 | 12 | 16 | 20 |
| --- | --- | --- | --- | --- | --- |
| states | 20 | 146 | 274 | 402 | 530 |
| states / $M$ | 5.0 | 18.3 | 22.8 | 25.1 | 26.5 |

Measured, same $(p,n)$, seconds:

| $M$ | 4 | 6 | 8 | 10 | 12 |
| --- | --- | --- | --- | --- | --- |
| per-profile | 0.0015 | 0.100 | 0.849 | 3.640 | 10.759 |
| shared-tail | 0.0013 | 0.023 | 0.090 | 0.266 | 0.651 |
| **folded (9.3)** | 0.0036 | 0.018 | 0.031 | 0.043 | **0.055** |
| shipped | 0.0014 | 0.0014 | 0.0014 | 0.0014 | 0.0014 |

195× over the per-profile form and 11.8× over shared-tail at $M = 12$, the folded cost rising about
$0.006$ s per extra facet while the others multiply. Against the other routes:

| $\{p,M,n\}$ | §9.2 naive DP | §9.2 convolution | **(9.3) folded** | shipped |
| --- | --- | --- | --- | --- |
| $\{2,5,7\}$ | 0.067 s | 0.0062 s | 0.0021 s | 0.00055 s |
| $\{3,5,5\}$ | 0.460 s | 0.0209 s | 0.0012 s | 0.00048 s |
| $\{3,5,7\}$ | 0.443 s | 0.0218 s | 0.0040 s | 0.00082 s |
| $\{3,6,8\}$ | 17.77 s | 0.142 s | 0.0104 s | 0.00108 s |

**323× faster than §9.2's DP at $(3,6,8)$**, and linear where that one is exponential.

In $n$ the state count *saturates* rather than tracking $P(n)$ — 65, 97, 119, 124, 125 for
$n = 8\dots16$ at $(p,M) = (3,6)$, against $\sum_{k \le n} P(k) = 67, 139, 272, 508, 915$ — because a
partition needs enough facets to be reachable, so $M$ caps it. The corollary is that the folded form
is **not** uniformly best: at $\{3,6,12\}$, small $M$ against large $n$, shared-tail wins
($0.040$ s against $0.064$ s). §9.3.8 maps that boundary.

It does not threaten the shipped counter in the regime that counter is built for — sweeping $M$ at
fixed $(p,n)$, where shipped is flat and this is linear, 20–40× behind. Its value is as the fast
independent check §9.2 wanted to be.

#### 9.3.8 Two crossover maps

Neither of the two fast forms dominates, and neither does the shipped counter uniformly. Both
boundaries were measured rather than reasoned, on a cold cache throughout (see the caution in §5
about which symbol actually clears it).

**Map 1: folded against shared-tail.** Entries are the ratio (folded / shared-tail); below $1$ the
folded form wins. A dot marks an infeasible $(M,n)$.

$p = 2$:

| $M \backslash n$ | 6 | 8 | 10 | 12 | 14 | 16 | 18 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 0.09 | . | . | . | . | . | . |
| 4 | 1.25 | 3.29 | . | . | . | . | . |
| 5 | 0.82 | 1.65 | 3.57 | . | . | . | . |
| 6 | 0.58 | 0.82 | 1.69 | 2.55 | . | . | . |
| 7 | 0.47 | 0.48 | 0.85 | 1.46 | 1.33 | . | . |
| 8 | 0.38 | 0.28 | 0.40 | 0.73 | 0.89 | 0.60 | . |
| 9 | 0.32 | 0.17 | 0.20 | 0.35 | 0.56 | 0.49 | 0.29 |
| 10 | 0.27 | 0.12 | 0.10 | 0.13 | 0.29 | 0.33 | 0.27 |

$p = 3$:

| $M \backslash n$ | 6 | 8 | 10 | 12 | 14 | 16 | 18 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 1.23 | 2.38 | . | . | . | . | . |
| 4 | 1.17 | 1.97 | 4.79 | 11.07 | . | . | . |
| 5 | 0.91 | 1.10 | 1.94 | 4.50 | 14.61 | . | . |
| 6 | 0.81 | 0.73 | 0.93 | 1.61 | 3.58 | 9.90 | 13.26 |
| 7 | 0.68 | 0.52 | 0.55 | 0.77 | 1.30 | 2.57 | 5.91 |
| 8 | 0.61 | 0.38 | 0.33 | 0.40 | 0.63 | 1.10 | 1.96 |
| 9 | 0.58 | 0.28 | 0.20 | 0.20 | 0.29 | 0.53 | 0.98 |
| 10 | 0.53 | 0.22 | 0.12 | 0.11 | 0.13 | 0.22 | 0.47 |

The boundary is close to linear and **nearly independent of $p$**: folded wins above roughly

$$M^\ast \;\approx\; \tfrac{n}{2} + 2,$$

and the useful way to read that is against the feasibility floor $M \ge n/p$. At $p = 2$ the floor
*is* $n/2$, so the shared-tail region is a sliver one or two facets wide and folded wins almost
everywhere. At $p = 3$ the floor drops to $n/3$ while the boundary stays near $n/2$, opening a real
band — every $M$ from $6$ to $8$ at $n = 18$ — where shared-tail is the right choice. Higher $p$
widens that band further. Deep in the folded region the margin is large: $8\times$ at $(2,10,10)$,
$9\times$ at $(3,10,12)$.

**Map 2: folded against shipped, in $n$ at fixed $M$.** $p = 3$, ratio folded / shipped, cold:

| $n$ | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $M = 6$ | 2.3 | 6.6 | 13.8 | 21.6 | **26.8** | 24.3 | 18.9 | . | . |
| $M = 8$ | 3.2 | 10.8 | 22.0 | 42.4 | 74.0 | 141.4 | 235.3 | **274.4** | 251.6 |

**There is essentially no crossover.** The shipped counter wins across the whole feasible range; the
single exception found is the degenerate corner $(3,4,4)$, ratio $0.5$. What the map does show is
that the ratio is **not monotone** — it rises, peaks, and falls, the peak sitting near
$n \approx 0.7\,pM$ to $0.8\,pM$ ($n = 8$ of $12$ at $M = 4$, $13$ of $18$ at $M = 6$, $19$ of $24$
at $M = 8$).

The decline past the peak is §9.3.7's saturation seen from the other side: the folded state count
stops growing once $M$ caps the reachable partitions, while the shipped cost keeps climbing with the
multiplicity-vector count. So the gap narrows as $n \to pM$ — but within $n \le pM$ it never closes,
and the honest summary is that this route is a fast independent check, not a replacement.

Two cautions on reading these numbers. The shipped column is a **cold** single call, which is the
fair comparison for one evaluation but not the access pattern §7 is built for; warm, its memoised
table answers in $10$–$30\,\mu$s and the ratios above multiply by roughly $40\times$. And the
folded form was measured with its own memo dropped between points, so neither side is being flattered.

#### 9.3.9 Four ways to get it wrong

All four were live in the first draft of the method, and none is caught by $s_F$ coming out right on
a single small case.

1. **Pad by $n_i$, not $n_{i+1}$.** The new vertices belong to the facet being placed. Padding one
   step late silently loses vertices — at $c = (3,2,0,0)$ it leaves 3 where $S_2 = 5$ — and every
   later stage is then built on a short state. Invisible whenever $n_i = n_{i+1}$, which is why a
   profile like $(3,1,1,0)$ hides it.
2. **Subtract $i-1$, not $i$.** There are $i-1$ earlier facets at stage $i$, and $M-1$ at the last.
3. **Drop the forbidden branches; do not subtract from each surviving term.** Subtracting $i-1$ from
   every term of the sum gives $\sum_{s} W(s) - (i-1)|MS|$, which at stage 3 of $W(3,2,0,0)$ is
   $19 - 15 = 4$ against the true $15$. The two coincide only at the leaves, where $W \equiv 1$ —
   which is why a terminal-only statement of the rule looks right.
4. **The all-or-nothing branches are the forbidden ones.** Not merely equinumerous with them: it is
   because they leave $\lambda$ fixed that the correction can be written without type identities at
   all (§9.3.3).

#### Verification

Against `NumFacetLabeledPureComplexes` on **180 parameter sets** across $p \le 4$, $M \le 6$,
$n \le 10$ — zero mismatches, to $s_F(4,6,10) = 20612880$ — for all three forms (per-profile,
shared-tail, folded). Plus the degenerate inputs of §6, and all 21 rows of §11's reference table.
The per-profile $W(c)$ were checked individually at $(3,4,5)$ against direct enumeration of the
isomorphism classes, which is what fixed $W(3,1,1,0) = 22$ and confirmed condition 3 rejects only
empty profiles.


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
