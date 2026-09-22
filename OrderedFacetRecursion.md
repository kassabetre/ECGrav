# Counting facet-labeled pure complexes by an ordered deletion recursion

> **What this is.** A self-contained derivation of a recursion that counts facet-labeled pure
> simplicial complexes by placing their facets in order, carrying an integer partition as its state.
> It is written for **P1** and assumes nothing from the package specifications; the same material
> appears as §9.3 of `FacetLabeledCount.md`, which states it against that document's notation.
>
> Two properties make it worth reporting rather than merely implementing. It is **linear in the
> facet order $M$**, where the natural dynamic programme over label content is exponential in $M$.
> And it **never forms a cycle type**, so it is independent of the Burnside/cycle-index route by
> which this count is normally obtained — which makes it usable as a check on that route rather
> than a restatement of it.

---

## 1. The object, and notation

Throughout, $p \ge 1$ is the **purity** (every facet has exactly $p$ vertices), $M \ge 0$ the
**facet order** (the number of facets), and $n \ge 0$ the **vertex count**. Write
$[n] = \{1,\dots,n\}$ and let $S_n$ be the symmetric group on $[n]$; for an arbitrary finite set
$\Omega$ write $\mathrm{Sym}(\Omega)$ for the symmetric group on $\Omega$, so
$S_n = \mathrm{Sym}([n])$. For a
partition $\lambda$ let $m_k(\lambda)$ be its number of parts equal to $k$, $\ell(\lambda)$ its number
of parts, and $|\lambda|$ the sum of its parts; $P(k)$ denotes the number of partitions of $k$. We write
partitions as weakly decreasing tuples, $(2,2,1)$, with $(1^5)$ abbreviating $(1,1,1,1,1)$. Square brackets carry three unrelated meanings in what follows, all
standard: $[n]$ as above; $[x^k]f$ for the coefficient of $x^k$ in a series $f$; and the **Iverson
bracket** $[\,P\,]$, equal to $1$ when the proposition $P$ holds and $0$ when it does not, so that
$[\,|\lambda| = n\,]$ is $1$ exactly on the partitions of $n$.

Let

$$X \;=\; \bigl\{(F_1,\dots,F_M) \;:\; F_i \subseteq [n],\ |F_i| = p,\ F_i \neq F_j \text{ for } i \neq j,\ \textstyle\bigcup_i F_i = [n]\bigr\}$$

be the set of ordered $M$-tuples of pairwise **distinct** $p$-subsets of $[n]$ whose union is all of
$[n]$. The group $S_n$ acts on $X$ by relabelling vertices, $\sigma\cdot(F_1,\dots,F_M) =
(\sigma F_1,\dots,\sigma F_M)$. A **facet-labeled pure complex** is an orbit of this action: the
facets carry the labels $1,\dots,M$, the vertices carry none. Write

$$s_F(p,M,n) \;=\; |X / S_n|$$

for the number of such orbits. Two conventions are worth stating once. *Distinct* is a condition on
the facets as sets, so $X$ excludes repeats; and *covering* means no vertex of $[n]$ is isolated, so
$n$ is the number of vertices actually used.

An **isomorphism** of facet-labeled complexes is a vertex bijection $\phi$ with
$\phi(F_i) = F_i'$ for **each** $i$ — the labels are matched, not permuted. This is the notion of
sameness used throughout; permuting the facet labels as well gives the fully unlabeled count, which
is a different and harder problem not treated here.

---

## 2. Incidence tableaux, types, and the state

### 2.1 Types

Fix a partial configuration $(F_1,\dots,F_i)$ — a tuple of $i$ distinct $p$-subsets, not yet
required to cover — and let $V_i = F_1 \cup \dots \cup F_i$ be the vertices used so far. Give each
$v \in V_i$ its **type**

$$\mathrm{ty}(v) \;=\; \{\, j \le i \;:\; v \in F_j \,\} \;\subseteq\; [i],$$

a **subset of $[i]$**, not a number. Types are nonempty on $V_i$ by construction.

Two different partitions now appear, and keeping them apart matters. Grouping $V_i$ by type — two
vertices in the same block exactly when they have the *same* type, so the blocks are the fibres of
the map $v \mapsto \mathrm{ty}(v)$ — gives a **set partition of $V_i$**, whose blocks we call the
**type classes**. Recording only the *sizes* of those blocks gives an **integer partition of
$|V_i|$**, written $\lambda(F_1,\dots,F_i)$. The set partition is the structure; the integer
partition is the state the recursion carries.

*Example.* Take $p = 3$ and $(F_1,F_2) = (\{1,2,3\},\{1,2,4\})$, so $V_2 = \{1,2,3,4\}$. Then
$\mathrm{ty}(1) = \mathrm{ty}(2) = \{1,2\}$, $\mathrm{ty}(3) = \{1\}$ and
$\mathrm{ty}(4) = \{2\}$, so the set partition is

$$V_2 \;=\; \underbrace{\{1,2\}}_{\text{type }\{1,2\}} \ \sqcup\ \underbrace{\{3\}}_{\text{type }\{1\}} \ \sqcup\ \underbrace{\{4\}}_{\text{type }\{2\}},$$

three blocks, and $\lambda(F_1,F_2) = (2,1,1)$, an integer partition of $4$.

What $\lambda$ discards is which class is which. For $(F_1,F_2) = (\{1,2,3\},\{1,4,5\})$ the
classes are $\{1\}$ of type $\{1,2\}$, $\{2,3\}$ of type $\{1\}$ and $\{4,5\}$ of type
$\{2\}$, giving $\lambda = (2,2,1)$ — two classes of equal size lying in different facets, which
$\lambda$ cannot tell apart. Proposition 1 is exactly the statement that it need not.

Listing the types of all vertices, as a multiset, is the **incidence tableau**: one row per vertex,
each row a subset of $[i]$. It is a complete isomorphism invariant — two configurations are
isomorphic exactly when their row multisets agree — since an isomorphism is precisely a bijection
matching rows.

### 2.2 The automorphism group, and why the partition is enough

Write $\mathrm{Aut}(F_1,\dots,F_i) = \{\sigma \in \mathrm{Sym}(V_i) : \sigma(F_j) = F_j \text{ for
all } j \le i\}$.

**Why an automorphism group appears in a count of isomorphism classes.** We never enumerate the
classes directly. They are built one facet at a time, and the only quantity each step needs is *how
many inequivalent ways the next facet can be added*. That is an orbit count, because two extensions
of the **same** configuration are isomorphic exactly when some automorphism of that configuration
carries one added facet to the other: an isomorphism of $(i{+}1)$-tuples must match $F_j$ to $F_j$
for every $j \le i$, so its restriction to $V_i$ lies in $\mathrm{Aut}(F_1,\dots,F_i)$. So
$\mathrm{Aut}$ of the **parent** classifies the **children**.

> **The group is the one being extended, not the one produced.** Take $F_1 = \{1,2,3\}$ at $p=3$,
> and its two extensions $A = (\{1,2,3\},\{1,4,5\})$ and $B = (\{1,2,3\},\{2,4,5\})$. These are
> isomorphic, via $(1\,2)$ — and $(1\,2) \in \mathrm{Aut}(F_1) = \mathrm{Sym}(\{1,2,3\})$, which
> is the group governing this step. It is *not* in $\mathrm{Aut}(A)$, and that is no objection:
> $\mathrm{Aut}(A)$ governs the step *after*, when a third facet is added to $A$. The recursion
> rolls, each configuration's automorphism group classifying the next facet's options and the
> extended configuration then supplying its own.
>
> Carried out for $F_1 = \{1,2,3\}$, the three orbits are $|F_1 \cap F_2| = 2, 1, 0$ —
> representatives $(\{1,2,3\},\{1,2,4\})$, $(\{1,2,3\},\{1,4,5\})$, $(\{1,2,3\},\{4,5,6\})$ —
> which is precisely the branch list out of $\lambda = (3)$ in §4, the fourth possibility
> $|F_1 \cap F_2| = 3$ being barred by distinctness. Both $A$ and $B$ represent the middle orbit.

**Proposition 1.** *Let $T_1,\dots,T_\ell$ be the type classes of $(F_1,\dots,F_i)$, so $\ell = \ell(\lambda)$. Then*

1. $\mathrm{Aut}(F_1,\dots,F_i) \;=\; \mathrm{Sym}(T_1) \times \dots \times \mathrm{Sym}(T_\ell)$,
   *the permutations preserving each type class setwise;*
2. *the orbits of $\mathrm{Aut}$ on $k$-element subsets of $V_i$ are in bijection with the vectors*
   $$MS_k(\lambda) \;=\; \bigl\{\, a = (a_1,\dots,a_\ell) \;:\; 0 \le a_u \le |T_u|,\ \textstyle\sum_u a_u = k \,\bigr\},$$
   *a $k$-element subset $S \subseteq V_i$ corresponding to the vector with $a_u = |S \cap T_u|$,
   its intersection sizes with the type classes;*
3. *two extensions of $(F_1,\dots,F_i)$ by one further facet are isomorphic as facet-labeled
   $(i{+}1)$-tuples if and only if they have the same $a$ and reuse the same number of new vertices.*

*Proof.* (1) $\sigma$ fixes every $F_j$ setwise iff for every $v$ and every $j$, $v \in F_j
\Leftrightarrow \sigma(v) \in F_j$, i.e. iff $\mathrm{ty}(\sigma(v)) = \mathrm{ty}(v)$ for all $v$ —
which is exactly preservation of the type classes. (2) Immediate from (1): a product of symmetric
groups acting on subsets has orbits classified by the intersection sizes, and every vector $a$ in
range is realised. (3) Let $F_{i+1} = S \sqcup D$ with $S = F_{i+1} \cap V_i$ and $D$ the new
vertices, and similarly $F_{i+1}' = S' \sqcup D'$. Any isomorphism $\phi$ of the $(i{+}1)$-tuples
maps $F_j \to F_j$ for $j \le i$, hence maps $V_i$ onto itself, so $\phi|_{V_i} \in \mathrm{Aut}$;
and $\phi(S) = S'$, so $S$ and $S'$ lie in one $\mathrm{Aut}$-orbit and $a = a'$. Conversely, given
$a = a'$ and $|D| = |D'|$, compose an element of $\mathrm{Aut}$ carrying $S$ to $S'$ with any
bijection $D \to D'$. $\square$

**Corollary 1.1.** $\displaystyle |\mathrm{Aut}(F_1,\dots,F_i)| \;=\; \prod_{u=1}^{\ell} |T_u|! \;=\; \prod_u \lambda_u!$

The factors are independent because the classes are disjoint, so an automorphism is one free
permutation per class. Classes of size $1$ contribute $1!= 1$ and impose no choice at all; a
configuration whose classes are all singletons has trivial $\mathrm{Aut}$. Note the order depends
only on $\lambda$ — a first indication that the integer partition is the right state. The recursion itself uses only the **orbits** of
$\mathrm{Aut}$ and never its order, so this corollary is orientation rather than machinery.

*The two examples of §2.1.* For $(\{1,2,3\},\{1,2,4\})$, with classes $\{1,2\}, \{3\}, \{4\}$,
only the first is free:
$$\mathrm{Aut} \;=\; \mathrm{Sym}(\{1,2\}) \times \mathrm{Sym}(\{3\}) \times \mathrm{Sym}(\{4\}) \;=\; \{\,e,\ (1\,2)\,\}, \qquad 2!\,1!\,1! = 2 .$$
Read from the constraints instead: $\sigma$ must fix $F_1 \cap F_2 = \{1,2\}$,
$F_1 \setminus F_2 = \{3\}$ and $F_2 \setminus F_1 = \{4\}$ setwise, and the last two being
singletons force $\sigma(3) = 3$, $\sigma(4) = 4$. For $(\{1,2,3\},\{1,4,5\})$, with classes
$\{1\}, \{2,3\}, \{4,5\}$, two classes are free and combine independently:
$$\mathrm{Aut} \;=\; \{\,e,\ (2\,3),\ (4\,5),\ (2\,3)(4\,5)\,\}, \qquad 1!\,2!\,2! = 4 .$$

> **This is not the automorphism group of the underlying unlabeled complex.** $\mathrm{Aut}$ here
> fixes each $F_j$ *individually*, because the facets carry labels. In the second example a vertex
> permutation swapping the two triangles exists, but it sends $F_1 \to F_2$ and is therefore
> excluded. The distinction is the same one that makes the facet-labeled count per isomorphism class
> $M!/|H|$ with $H$ the image of the full automorphism group in $\mathrm{Sym}\{F_1,\dots,F_M\}$,
> rather than $M!$ divided by that group's own order.

**What the action is.** $\mathrm{Aut}$ acts on $\binom{V_i}{k}$, the set of all $k$-element
subsets, *elementwise*: each $\sigma$ sends each individual subset $S$ to $\sigma(S)$, so $\sigma$
induces a permutation of those $\binom{|V_i|}{k}$ objects, and the orbits are the cycles of that
induced permutation group — a partition of $\binom{V_i}{k}$ into classes of mutually reachable
subsets. It is a different and vacuous statement that $\mathrm{Aut}$ carries the *collection*
$\binom{V_i}{k}$ onto itself; every permutation group does, which is what makes the action well
defined rather than what its orbits are.

*Example, $(F_1,F_2) = (\{1,2,3\},\{1,2,4\})$ at $k = 2$.* Here $\mathrm{Aut} = \{e,(1\,2)\}$, and
$(1\,2)$ moves the six $2$-subsets by
$\{1,3\} \leftrightarrow \{2,3\}$ and $\{1,4\} \leftrightarrow \{2,4\}$, fixing $\{1,2\}$ and
$\{3,4\}$. The four orbits are therefore

| orbit | $\{1,2\}$ | $\{1,3\},\{2,3\}$ | $\{1,4\},\{2,4\}$ | $\{3,4\}$ |
| --- | --- | --- | --- | --- |
| $a$ | $(2,0,0)$ | $(1,1,0)$ | $(1,0,1)$ | $(0,1,1)$ |
| size | $1$ | $2$ | $2$ | $1$ |

summing to $6 = \binom{4}{2}$. Two features are worth pausing on. First, $\{1,3\}$ and $\{1,4\}$
lie in *different* orbits although both take one vertex from the large class and one from a
singleton: $a$ records **which** class, not merely how big it was, and no identity is discarded at
this stage. Second, orbits need not have equal size — for $(\{1,2,3\},\{1,4,5\})$ at $k = 2$ they
are $1,1,2,2,4$ — since each has size $|\mathrm{Aut}|$ divided by a stabiliser that varies from
orbit to orbit. The recursion counts orbits and never weights them, so this does not intrude.

**When is there a single orbit?** Exactly when $\mathrm{Aut}$ is transitive on $k$-subsets. Setting
aside the trivial ends $k = 0$ and $k = |V_i|$, where there is only one $k$-subset to be transitive
on, a Young subgroup is transitive only when $\ell = 1$, i.e. $\mathrm{Aut} = \mathrm{Sym}(V_i)$:
the orbit count $|MS_k(\lambda)|$ equals $\ell$ at $k = 1$ and, its coefficient sequence being
symmetric and unimodal, is no smaller anywhere between. In the example
$|\mathrm{Aut}| = 2$ forces every orbit to have size $1$ or $2$ and hence at least $6/2 = 3$ of
them; concretely no automorphism can carry $\{1,2\} = F_1 \cap F_2$ to $\{3,4\}$, which meets it
not at all. One facet earlier the picture is the opposite: for $(F_1) = (\{1,2,3\})$ there is a
single type class, $\mathrm{Aut} = \mathrm{Sym}(\{1,2,3\})$, and $\binom{V_1}{2}$ *is* one orbit,
matching $MS_2\bigl((3)\bigr) = \{(2)\}$. That is the general shape — each facet refines the type
classes, shrinks $\mathrm{Aut}$ and splits orbits:

$$\lambda = (3),\ \ell = 1,\ |\mathrm{Aut}| = 6,\ 1 \text{ orbit}
\quad\xrightarrow{\ +\,F_2 = \{1,2,4\}\ }\quad
\lambda = (2,1,1),\ \ell = 3,\ |\mathrm{Aut}| = 2,\ 4 \text{ orbits.}$$

The recursion starts at the one-orbit end, with $V\bigl((p),\,M-1\bigr)$, and the growth in the
number of orbits is its branching factor.

Part (2) is what licenses discarding the *identities* of the types and keeping only the multiset of
their sizes. The count of choices is

$$\bigl| MS_k(\lambda) \bigr| \;=\; [x^k] \prod_{u=1}^{\ell} \bigl(1 + x + \dots + x^{\lambda_u}\bigr),$$

a function of $\lambda$ alone. Part (3) says distinct choices give non-isomorphic extensions, so no
correction for over- or under-counting is needed when we branch over $MS_k(\lambda)$.

### 2.3 The transition

Having chosen $a \in MS_{p-t}(\lambda)$ — reuse $a_u$ vertices from type class $T_u$ — and $t$ new
vertices, the new type classes are: for each $u$, the $a_u$ reused vertices (whose type gains the new
facet) and the $\lambda_u - a_u$ untouched ones (whose type does not), plus one class of the $t$ new
vertices. A class **splits** exactly when $0 < a_u < \lambda_u$; a fully-taken or wholly-untaken
class does not. So define, writing $\uplus$ for **multiset** union — union that adds multiplicities,
so that $\{1\} \uplus \{1\} = \{1,1\}$ and not $\{1\}$ —

$$\mathrm{grow}(\lambda, a, t) \;=\; \biguplus_{u \,:\, a_u > 0} \{\, a_u \,\} \;\;\uplus\;\; \biguplus_{u \,:\, \lambda_u - a_u > 0} \{\, \lambda_u - a_u \,\} \;\;\uplus\;\; \{\, t : t > 0 \,\},$$

a partition of $|\lambda| + t$. The multiplicities carry real weight, which is why the first two
groups are indexed by $u$ rather than written as plain sets: equal-sized classes are still distinct
classes and each owes its own part. At $\lambda = (2,1,1)$ with $a = (1,1,1)$ and $t = 0$, the three
taken parts and the one leftover part give $\mathrm{grow} = (1,1,1,1)$, a partition of $4$; collapsing
duplicates would return a single $1$ and lose three vertices. The last group needs no such care,
holding one element or none.

Note that $\mathrm{grow}$ depends only on $\lambda$, $a$ and $t$, which together with Proposition 1
is the whole reason a partition suffices as state.

---

## 3. The ordered decomposition

### 3.1 The new-vertex profile

For a facet-labeled complex, set

$$n_i \;=\; \bigl|\,F_i \setminus (F_1 \cup \dots \cup F_{i-1})\,\bigr|, \qquad i = 1,\dots,M,$$

the number of vertices facet $i$ introduces that no earlier facet used, and call
$c = (n_1,\dots,n_M)$ the **new-vertex profile**. It is an isomorphism invariant, since a vertex
bijection matching each $F_i$ preserves every $F_i \setminus \bigcup_{j<i}F_j$. Hence:

**Proposition 2.** *Every isomorphism class has exactly one profile, so the profiles partition the
classes and*
$$s_F(p,M,n) \;=\; \sum_{c} W(c),$$
*where $W(c)$ is the number of classes with profile $c$. Moreover any realised profile satisfies*

1. $n_1 = p$ *(facet 1 is entirely new);*
2. $0 \le n_i \le p$ *(a facet introduces at most its own size);*
3. $\sum_i n_i = n$ *(the facets cover $[n]$, and the sets $F_i \setminus \bigcup_{j<i}F_j$ are
   disjoint with union $[n]$);*
4. $\binom{|V_i|}{p} \ge i$ *with $|V_i| = \sum_{j \le i} n_j$ (after $i$ facets the pool $V_i$ must
   supply $i$ pairwise distinct $p$-subsets).*

Write $C(p,M,n)$ for the profiles meeting (1)–(4). Condition (4) is a genuine restriction and is
sharp in the sense that it discards only empty cases: at $(p,M,n) = (3,4,5)$ it rejects exactly
$(3,0,2,0)$, $(3,0,1,1)$ and $(3,0,0,2)$, each of which has $\binom{3}{3} = 1 < 2$ at $i = 2$, and
direct enumeration confirms each supports no classes. The three survivors carry

| $c$ | $(3,1,0,1)$ | $(3,1,1,0)$ | $(3,2,0,0)$ | total |
| --- | --- | --- | --- | --- |
| $W(c)$ | 6 | 22 | 15 | **43** $= s_F(3,4,5)$ |

### 3.2 Distinctness

Placing facets in order makes the distinctness condition local, and cheap.

**Proposition 3.** *Let $(F_1,\dots,F_{i-1})$ have type classes $T_1,\dots,T_\ell$ and partition
$\lambda$, and consider extending it by $F_i$ with $n_i$ new vertices.*

1. *If $n_i > 0$ then $F_i \neq F_j$ automatically for every $j < i$.*
2. *If $n_i = 0$ then exactly $i - 1$ of the branches $a \in MS_p(\lambda)$ are forbidden, namely
   those with $F_i = F_j$ for some $j < i$.*
3. *Each forbidden branch is **all-or-nothing** — $a_u \in \{0, \lambda_u\}$ for every $u$ — and
   therefore satisfies $\mathrm{grow}(\lambda, a, 0) = \lambda$.*

*Proof.* (1) $F_i$ contains a vertex outside $V_{i-1} \supseteq F_j$. (2) and (3): for $j < i$ the
set $F_j$ is a $p$-subset of $V_{i-1}$, so it is one of the branches, and its selection vector is
$a_u = \lambda_u$ when $j \in \mathrm{ty}(T_u)$ and $a_u = 0$ otherwise — a vertex of $V_{i-1}$ lies
in $F_j$ iff its type contains $j$, and type is constant on a class. That is all-or-nothing, and an
all-or-nothing selection splits no class, so the partition is unchanged. Distinctness of the
forbidden branches: if $F_j \neq F_{j'}$ then some vertex lies in one and not the other, so its type
class has $j$ in its type and not $j'$ (or conversely), and the two selection vectors differ at that
class. Hence the $i-1$ earlier facets give $i-1$ distinct forbidden branches. $\square$

Part (3) is the key to the whole construction. Knowing *which* branches are forbidden would require
the type identities, which the partition state has discarded — but we never need to know. All $i-1$
forbidden branches lead to the **same** successor state $\lambda$, so their removal is a single
subtraction of $(i-1)$ copies of one term.

### 3.3 The profile weight

Proposition 2 reduces $s_F$ to the weights $W(c)$ without saying how to obtain one. Propositions 1
and 3 now supply that, and the method is complete.

Fix $c \in C(p,M,n)$ and place the facets in order. Facet 1 may be taken to be any $p$-set, all
choices isomorphic, giving the state $\lambda = (p)$. At stage $i \ge 2$ the profile already fixes how
many vertices are new, so the only choice left is which old ones to reuse: $p - n_i$ of them,
selected by some $a \in MS_{p-n_i}(\lambda)$, after which the state becomes
$\mathrm{grow}(\lambda,a,n_i)$. Note that $|\lambda| = \sum_{j \le i} n_j$ after stage $i$, fixed by
the profile — so unlike (4.1) there is no vertex budget to carry, which is the one respect in which
the unfolded form is the simpler of the two.

Write $W_c(\lambda, i)$ for the number of ways to place facets $i, i+1, \dots, M$ from state
$\lambda$, counted up to isomorphism fixing the existing configuration.

**Proposition 3.1.** *Define $W_c$ at every pair by $W_c(\lambda,\, M+1) = 1$ and, for
$2 \le i \le M$,*

$$W_c(\lambda,\, i) \;=\; \sum_{a \,\in\, MS_{p-n_i}(\lambda)} W_c\bigl(\mathrm{grow}(\lambda,a,n_i),\ i+1\bigr) \;-\; [\,n_i = 0\,]\;(i-1)\,W_c(\lambda,\ i+1). \tag{3.1}$$

*Then $W_c(\lambda,i)$ is the count just described at every realisable $(\lambda,i)$, and*

$$W(c) \;=\; W_c\bigl((p),\ 2\bigr), \qquad\qquad s_F(p,M,n) \;=\; \sum_{c \,\in\, C(p,M,n)} W_c\bigl((p),\ 2\bigr). \tag{3.2}$$

*Proof.* Downward induction on $i$, over states realisable by $i-1$ pairwise distinct facets with
profile prefix $(n_1,\dots,n_{i-1})$. At $i = M+1$ all $M$ facets are placed; condition (3) of
Proposition 2 gives $|\lambda| = \sum_{j \le M} n_j = n$, so the configuration is complete, covers by
construction, and counts once.

For $i \le M$, let $\mathcal{F}$ realise $(\lambda, i)$. Facet $i$ brings $n_i$ new vertices and
reuses $p - n_i$ old ones, so by Proposition 1(3) the isomorphism classes of such extensions are
exactly the $a \in MS_{p-n_i}(\lambda)$, the extension by $a$ having partition
$\mathrm{grow}(\lambda,a,n_i)$. A **legal** extension leaves $i$ pairwise distinct facets, so its
successor state is realisable and the induction hypothesis applies to it; the count sought is the sum
over legal $a$ of $W_c(\mathrm{grow}(\lambda,a,n_i), i+1)$.

That (3.1) computes it is Proposition 3. When $n_i > 0$ no extension violates distinctness and the
correction is absent. When $n_i = 0$ exactly $i-1$ do, every one with successor state $\lambda$, so
together they contribute $(i-1)\,W_c(\lambda, i+1)$ to the sum — which is what the correction
removes, term by term, using nothing about whether $(\lambda, i+1)$ is itself realisable. $\square$

The same caution as in §4 applies for the same reason: off the realisable states $W_c$ is an
intermediate that exists to be cancelled, and it may be negative.

§5.1 runs (3.1) through for $c = (3,2,0,0)$ at $(p,M,n) = (3,4,5)$, and §5.2 does the same complex by
(4.1) for comparison. The three profiles of $C(3,4,5)$ give $W = 6, 22, 15$ as tabulated in §3.1, and
the three that condition (4) rejects all return $0$ — as §4 explains they must.

**What this costs.** One tree per profile, with no sharing between them, and $|C(p,M,n)|$ grows
exponentially in $M$: profiles are compositions of $n - p$ into $M - 1$ parts each at most $p$. §4
removes the enumeration altogether by folding the choice of $n_i$ into the recursion, at which point
the profile stops being an object the algorithm handles at all.

---

## 4. The recursion

Fix $p$, $M$ and $n$. For a partition $\lambda$ with $|\lambda| \le n$ and an integer $0 \le j \le
M$, let

$$V(\lambda, j) \;=\; \text{the number of ways to place } j \text{ further facets},$$

counted up to isomorphism fixing the existing configuration, so that the resulting complex has
exactly $n$ vertices, covers them, and has all $M$ facets pairwise distinct — starting from any
configuration of $M - j$ facets whose type-multiplicity partition is $\lambda$. Propositions 1 and 3
say this is well defined: the branching and the size of the forbidden set depend on $(\lambda, j)$
alone, not on the configuration realising them.

Call $(\lambda, j)$ **realisable** when some configuration of $M - j$ pairwise distinct facets does
have partition $\lambda$. The description above defines $V$ there and nowhere else, while (4.1)
below, read as a definition, returns a number at every pair; that the two agree wherever both mean
anything is part of what Theorem 4 asserts.

**Theorem 4.** *Define $V$ at every pair by $V(\lambda, 0) = [\,|\lambda| = n\,]$ and, for $j > 0$,
with $r = n - |\lambda|$ the remaining vertex budget,*

$$\begin{aligned}
V(\lambda, j) \;=\;\ & \Biggl(\ \sum_{a \,\in\, MS_{p}(\lambda)} V\bigl(\mathrm{grow}(\lambda,a,0),\ j-1\bigr)\ -\ (M-j)\,V(\lambda,\ j-1)\ \Biggr) \\
& +\ \sum_{t=1}^{\min(p,\,r)}\ \sum_{a \,\in\, MS_{p-t}(\lambda)} V\bigl(\mathrm{grow}(\lambda,a,t),\ j-1\bigr),
\end{aligned} \tag{4.1}$$

*the parenthesised group being the $t = 0$ term — it alone carries the subtraction — and the inner
sum empty unless $0 \le p - t \le |\lambda|$. Then $V(\lambda, j)$ is the count described above at
every realisable $(\lambda, j)$, and*

$$s_F(p,M,n) \;=\; V\bigl((p),\ M-1\bigr). \tag{4.2}$$

*Proof.* Induction on $j$, over realisable states. At $j = 0$ nothing more is placed, so the
configuration is complete and is counted iff it already has $n$ vertices; covering holds because
every vertex of $V_M$ lies in some facet by construction.

For $j > 0$, let $\mathcal{F}$ realise $(\lambda, j)$, with its $i = M - j$ pairwise distinct
facets. The next facet reuses $p - t$ old vertices and brings $t$ new ones, $t$ ranging over
$0,\dots,\min(p,r)$ since it cannot exceed the budget. By Proposition 1(3) the isomorphism classes
of such extensions are exactly the pairs $(t, a)$ with $a \in MS_{p-t}(\lambda)$, the extension by
$(t,a)$ having partition $\mathrm{grow}(\lambda,a,t)$. By Proposition 3 the pairs violating
distinctness occur only at $t = 0$, number exactly $M - j$, and all have successor state $\lambda$.
A **legal** pair extends $\mathcal{F}$ to $i + 1$ pairwise distinct facets, so its successor state
is realisable and the induction hypothesis applies to it. The count sought is therefore

$$\sum_{(t,a)\ \mathrm{legal}} V\bigl(\mathrm{grow}(\lambda,a,t),\ j-1\bigr).$$

It remains to see that (4.1) computes that sum. Split its double sum into legal and forbidden pairs.
The forbidden pairs number $M - j$ and every one has successor state $\lambda$, so together they
contribute exactly $(M-j)\,V(\lambda, j-1)$ — which is what the subtraction removes, leaving the sum
over legal pairs alone. This cancellation is term by term and uses nothing about the *number*
$V(\lambda, j-1)$, in particular not that $(\lambda, j-1)$ be realisable. Finally (4.2): facet 1 may
be taken to be any $p$-set, all choices isomorphic, giving the single state $(p)$ on $p$ vertices
with $M-1$ facets left, which is realisable. $\square$

> **Off the realisable states $V$ means nothing, and can be negative.** The recursion does reach
> such pairs: a forbidden branch's successor is $(\lambda, j-1)$, the state of a configuration with
> a repeated facet, which need not be realisable by distinct ones. At $(p,M,n) = (3,4,5)$ it
> evaluates $V\bigl((3), 2\bigr)$ — partition $(3)$ with two facets already down, impossible, since
> three vertices admit only one $3$-subset. There $|MS_3((3))| = 1$ while $M - j = 2$, so the
> $t = 0$ term is $1 - 2 = -1$: strictly more is subtracted than is present, and Proposition 3 does
> not apply, having assumed a configuration that does not exist. Nothing is wrong, because that
> value is consumed only in the cancelling pair at the realisable parent $\bigl((3), 3\bigr)$, where
> the $t = 0$ term is $V((3),2) - V((3),2) = 0$. It is the cancellation, not any non-negativity,
> that makes the recursion safe — with a consequence for implementations, §9 item 5.

Three consequences are worth stating separately.

**The vertex budget is not a free coordinate.** It is $r = n - |\lambda|$, forced by the state, so
$V$ has two arguments and not three. This is what keeps the state space small.

**Condition (4) of Proposition 2 is not needed here.** There is no profile list left to prune. Nor
was the condition ever a matter of correctness in §3: the profiles it rejects carry $W(c) = 0$
anyway, and infeasible branches starve of their own accord — at $\lambda = (3)$ with one facet down,
$MS_3((3))$ has a single element and that element is the repeat of $F_1$, so the term is $1 - 1 = 0$.
Condition (4) buys speed, by not descending into subtrees that will cancel to nothing.

**Folding $t$ into the recursion removes the profile enumeration entirely.** The two organisations
walk the same objects: unrolled over the $M$ levels, the $t$ chosen at level $i$ is $n_i$, so a
root-to-leaf path of (4.1) *is* a profile. What changes is that a profile is a **path**, while the
memo is keyed on a **state** $(\lambda, j)$. Profiles arriving at a common state have identical
futures, which the profile-first organisation recomputes once per profile and (4.1) computes once —
and there are exponentially many profiles against polynomially many states. Profiles are the
*derivation*, not the algorithm.

---

## 5. Worked examples

### 5.1 One profile weight: $W(3,2,0,0) = 15$

Take $p = 3$, $M = 4$, and the profile $c = (3,2,0,0)$, so $n = 5$. The stages below are (3.1)
unrolled, one per facet.

**Stage 1.** Facet 1 is three new vertices. State $(3)$; types $\{1\}^3$.

**Stage 2**, $n_2 = 2$. Facet 2 reuses $p - n_2 = 1$ old vertex, and $|MS_1((3))| = 1$: the only
choice is one vertex from the single class. That class splits $1 + 2$, and the two new vertices form
a class of their own, giving state $(2,2,1)$ on $|V_2| = 5$ vertices. In types:
$\{1\}^2,\ \{2\}^2,\ \{1,2\}^1$.

**Stage 3**, $n_3 = 0$. Facet 3 is three old vertices; $|MS_3((2,2,1))| = 5$. Writing selection
vectors in the coordinates $(\{1\}, \{2\}, \{1,2\})$:

| $a$ | $\mathrm{grow}(\lambda,a,0)$ | |
| --- | --- | --- |
| $(2,1,0)$ | $(2,1,1,1)$ | |
| $(1,2,0)$ | $(2,1,1,1)$ | |
| $(1,1,1)$ | $(1^5)$ | |
| $(2,0,1)$ | $(2,2,1)$ | $= F_1$, forbidden |
| $(0,2,1)$ | $(2,2,1)$ | $= F_2$, forbidden |

The two forbidden branches are exactly the all-or-nothing ones and leave the state at $(2,2,1)$, as
Proposition 3 predicts. Three branches survive.

**Stage 4**, $n_4 = 0$. Each surviving state contributes its number of $3$-subsets less the $M - 1 =
3$ already-placed facets, using $|MS_3((2,1,1,1))| = 7$ and $|MS_3((1^5))| = \binom{5}{3} = 10$:

$$W(3,2,0,0) \;=\; (7-3) \;+\; (7-3) \;+\; (10-3) \;=\; 4 + 4 + 7 \;=\; 15 .$$

### 5.2 The whole count: $s_F(3,4,5) = 43$

The same parameters, now through (4.1), which never mentions a profile. By (4.2) the answer is
$V\bigl((3), 3\bigr)$; throughout, $n = 5$, $M = 4$, the base is $V(\lambda, 0) = [\,|\lambda| = 5\,]$,
and $M - j$ is the number of facets already placed, so $j = 3$ means one.

**The top level.** $\lambda = (3)$, so $r = n - 3 = 2$ and $M - j = 1$. A single type class admits one
selection vector per $t$:

| $t$ | $a \in MS_{3-t}\bigl((3)\bigr)$ | $\mathrm{grow}$ | contributes |
| --- | --- | --- | --- |
| $0$ | $(3)$ | $(3)$ | $V\bigl((3),2\bigr) - 1\cdot V\bigl((3),2\bigr) \;=\; 0$ |
| $1$ | $(2)$ | $(2,1,1)$ | $V\bigl((2,1,1),2\bigr) \;=\; 28$ |
| $2$ | $(1)$ | $(2,2,1)$ | $V\bigl((2,2,1),2\bigr) \;=\; 15$ |

$$s_F(3,4,5) \;=\; 0 + 28 + 15 \;=\; 43 .$$

Compare the three orbits listed in §2.2 for $F_1 = \{1,2,3\}$. Since $t = p - |F_1 \cap F_2|$, the
rows above are $|F_1 \cap F_2| = 3, 2, 1$; the fourth orbit, $|F_1 \cap F_2| = 0$, would need $t = 3$
and is cut by the budget $r = 2$, not by distinctness. The $t = 0$ row is the repeat $F_2 = F_1$, and
here the subtraction consumes the entire term, $MS_3((3))$ having exactly one element.

**The $t = 1$ subtree.** $\lambda = (2,1,1)$, $r = 1$, $M - j = 2$. The classes are $\{1,2\}$ of type
$\{1,2\}$, one vertex of type $\{1\}$ and one of type $\{2\}$ — the configuration
$(\{1,2,3\},\{1,2,4\})$ of §2.2. In those coordinates:

| $t$ | $a$ | $\mathrm{grow}$ | $V(\cdot,1)$ | |
| --- | --- | --- | --- | --- |
| $0$ | $(2,1,0)$ | $(2,1,1)$ | $4$ | $= F_1$, forbidden |
| $0$ | $(2,0,1)$ | $(2,1,1)$ | $4$ | $= F_2$, forbidden |
| $0$ | $(1,1,1)$ | $(1^4)$ | $6$ | |
| $1$ | $(2,0,0)$ | $(2,1,1,1)$ | $4$ | |
| $1$ | $(1,1,0)$ | $(1^5)$ | $7$ | |
| $1$ | $(1,0,1)$ | $(1^5)$ | $7$ | |
| $1$ | $(0,1,1)$ | $(2,1,1,1)$ | $4$ | |

The $t = 0$ rows sum to $14$, less $2\,V\bigl((2,1,1),1\bigr) = 8$, leaving $6$; the $t = 1$ rows sum
to $22$. So $V\bigl((2,1,1),2\bigr) = 6 + 22 = 28$. The subtraction has removed precisely the two
forbidden rows, each worth $4$ — visible here because they are the all-or-nothing selections and
return to $\lambda$.

**The $t = 2$ subtree.** $\lambda = (2,2,1)$ with $r = 0$, so only $t = 0$ survives, and this node is
§5.1's stage 3 read again. Its five branches carry $V(\cdot,1) = 4, 2, 4, 7, 2$, summing to $19$,
less $2\,V\bigl((2,2,1),1\bigr) = 4$, giving $V\bigl((2,2,1),2\bigr) = 15$. Note the bookkeeping
differs from §5.1 while the answer does not: §5.1 discards the forbidden branches before descending,
whereas (4.1) descends into them and subtracts them afterwards.

**The leaves**, all at $j = 1$, where $M - j = 3$ and a state counts $1$ iff it has five vertices:

| state | $r$ | value |
| --- | --- | --- |
| $V\bigl((1^5),1\bigr)$ | $0$ | $\binom{5}{3} - 3 = 10 - 3 = 7$ |
| $V\bigl((2,1,1,1),1\bigr)$ | $0$ | $\lvert MS_3\bigl((2,1,1,1)\bigr)\rvert - 3 = 7 - 3 = 4$ |
| $V\bigl((2,2,1),1\bigr)$ | $0$ | $5 - 3 = 2$ |
| $V\bigl((1^4),1\bigr)$ | $1$ | $t=0$ gives $0$; $t=1$ gives $\binom{4}{2} = 6$ |
| $V\bigl((2,1,1),1\bigr)$ | $1$ | $t=0$ gives $0$; $t=1$ gives $4$ |

At $\lambda = (1^5)$ the automorphism group is trivial, so orbits *are* subsets and the count is the
plain $\binom{5}{3}$ less the three facets already placed. Where $r = 1$ the $t = 0$ branches all
die at the base, having only four vertices with no further facet able to add a fifth.

**The state §4 warns about.** $V\bigl((3),2\bigr)$ is evaluated here, and its $t = 0$ term is
$1 - 2\,V\bigl((3),1\bigr) = 1 - 2 = -1$, with $t = 1$ and $t = 2$ giving $4$ and $2$, so the entry is
$5$. It is unrealisable — two distinct facets cannot leave one type class — and it enters the answer
only as $V\bigl((3),2\bigr) - V\bigl((3),2\bigr)$ in the top row above.

**The profile table, recovered.** Ten states are visited in all. Grouping §3.1's profiles by $n_2$
reproduces the top-level $t$-sum exactly:

| $t = n_2$ | profiles | $\sum W$ | top-level term |
| --- | --- | --- | --- |
| $0$ | $(3,0,2,0)$, $(3,0,1,1)$, $(3,0,0,2)$ | $0 + 0 + 0$ | $0$ |
| $1$ | $(3,1,0,1)$, $(3,1,1,0)$ | $6 + 22$ | $28$ |
| $2$ | $(3,2,0,0)$ | $15$ | $15$ |

and one level down the split continues: the $6$ and the $22$ are exactly the $t = 0$ and $t = 1$
halves of $V\bigl((2,1,1),2\bigr)$ computed above. This is what "a root-to-leaf path is a profile"
means concretely. The $t = 0$ row is also the claim of §4 in miniature — the three profiles
condition (4) would have discarded contribute $0$ without being identified, let alone excluded.

---

## 6. Complexity

Throughout this section $p$ and $n$ are fixed and $M$ varies; $O_{p,n}$ hides constants depending on
them but not on $M$. Write $\Pi(n) = \sum_{k=0}^{n} P(k)$ for the number of partitions of the
integers $0,\dots,n$, and $B(p,n) = (p+1)\binom{p+n-1}{n-1}$.

**Proposition 5.** *Assume $p \le n$. Then*

1. *every reachable $\lambda$ satisfies $|\lambda| \le n$, and the number of reachable states
   $(\lambda,j)$ lies between $M$ and $\Pi(n)\,M$, hence is $\Theta_{p,n}(M)$;*
2. *evaluating (4.1) at all of them takes $O_{p,n}(M)$ arithmetic operations;*
3. *the value returned has $O_{p,n}(M)$ digits and no intermediate exceeds $O_{p,n}(M \log M)$, so
   the cost in **bit** operations is $\tilde O_{p,n}(M^2)$ — linear in the count of arithmetic
   operations, but not in time.*

*Proof.* (1) The initial state is $\bigl((p), M-1\bigr)$, and $|(p)| = p \le n$. From $(\lambda, j)$
the successors are $\bigl(\mathrm{grow}(\lambda,a,t),\,j-1\bigr)$ with $t \le \min(p,r)$ and
$r = n - |\lambda|$; since $\mathrm{grow}(\lambda,a,t)$ is a partition of $|\lambda| + t$ and
$t \le r$, its size is at most $|\lambda| + r = n$. By induction every reachable $\lambda$ is a
partition of an integer in $[0,n]$, and there are $\Pi(n)$ of those. Each transition decrements $j$
by exactly one from $M-1$, so $j \in \{0,\dots,M-1\}$ takes $M$ values, giving at most $\Pi(n)\,M$
states. For the lower bound, take $t = 0$ and $a = (p)$ at $\lambda = (p)$: this lies in
$MS_p\bigl((p)\bigr)$ and has $\mathrm{grow}\bigl((p),(p),0\bigr) = (p)$, so $\bigl((p),j\bigr)$ is
reachable for every $0 \le j \le M-1$ — at least $M$ distinct states.

(2) At a state with $\ell = \ell(\lambda)$ parts, $MS_k(\lambda)$ is contained in the set of
compositions of $k$ into $\ell$ non-negative parts, so $|MS_k(\lambda)| \le \binom{k+\ell-1}{\ell-1}$.
Since $\ell \le |\lambda| \le n$ by (1) and $k \le p$, summing over $k = 0,\dots,p$ bounds the
branches at any state by $B(p,n)$, which does not involve $M$. Each branch costs one evaluation of
$\mathrm{grow}$ — at most $O(n\log n)$ comparisons — and one memo lookup; the $t = 0$ group costs one
further multiplication by the scalar $M-j$ and one subtraction. So every state is $O_{p,n}(1)$
arithmetic operations, and by (1) there are $\Theta_{p,n}(M)$ states.

(3) Orbits are no more numerous than the set acted on, so
$s_F(p,M,n) \le |X| \le \binom{n}{p}^{M}$ and the answer has $O(Mp\log n) = O_{p,n}(M)$ digits.
(The tempting sharper bound $\binom{\binom{n}{p}}{M}$, counting *sets* of facets, is false here:
$s_F$ counts ordered tuples up to $S_n$ only, and at $(p,M,n) = (2,6,5)$ it is $1230$ against
$\binom{10}{6} = 210$.) For the intermediates, (2) gives
$|V(\lambda,j)| \le \bigl(B(p,n) + M\bigr)\max_{\lambda'}|V(\lambda',j-1)|$ with $|V(\cdot,0)| \le 1$,
whence $|V(\lambda,j)| \le (B+M)^{M}$, of $O_{p,n}(M \log M)$ digits. The $O_{p,n}(M)$ operations of
(2) are therefore performed on integers of that length, and the multiplications are by a scalar below
$M$. $\square$

The gap between (2) and (3) is the point to carry away: **the operation count is linear, the running
time is not.** The answer's digit count grows linearly in $M$ on its own, so the work per operation
grows too.

**Why the $\Pi(n)\,M$ bound is loose.** It counts partitions that no configuration can present.

**Lemma 5.1.** *If $(\lambda,j)$ is reachable and $i = M - j$ facets have been placed, then
$|\lambda| \le \min(n,\,ip)$ and $\ell(\lambda) \le \min(n,\,2^i - 1)$.*

*Proof.* $|\lambda| = |V_i|$ counts the vertices used by $i$ facets of size $p$, so $|\lambda| \le ip$,
and $|\lambda| \le n$ by Proposition 5(1). The parts of $\lambda$ are the type classes, whose types
are non-empty subsets of $[i]$; there are at most $2^i - 1$ of those, and $\ell(\lambda) \le |\lambda|$
gives the other bound. $\square$

Both bounds stop increasing once $i \ge \max\bigl(\lceil n/p\rceil,\ \lceil \log_2(n+1)\rceil\bigr)$,
so from that depth on, the states per level are capped independently of $M$ — which is why the totals
below are linear rather than merely bounded by one. At $(p,n) = (3,9)$ the threshold is $i = 4$, and
the reachable $\lambda$ per level run $1, 4, 15, 30, 32, 32, 32, \dots$, constant from $i = 5$.

Instrumented at $(p,n) = (3,9)$:

| $M$ | 4 | 8 | 12 | 16 | 20 |
| --- | --- | --- | --- | --- | --- |
| states | 20 | 146 | 274 | 402 | 530 |
| states $/\,M$ | 5.0 | 18.3 | 22.8 | 25.1 | 26.5 |

with states$/M$ flattening toward a constant. Timings at the same $(p,n)$, in seconds, against the
per-profile form of §3 (one tree per element of $C(p,M,n)$, no sharing):

| $M$ | 4 | 6 | 8 | 10 | 12 |
| --- | --- | --- | --- | --- | --- |
| per-profile | 0.0015 | 0.100 | 0.849 | 3.640 | 10.759 |
| folded (4.1) | 0.0036 | 0.018 | 0.031 | 0.043 | **0.055** |

The folded cost rises by about $0.006$ s per additional facet where the per-profile form multiplies —
$195\times$ at $M = 12$, and widening.

Over this range the folded column looks linear, and Proposition 5(3) says it cannot stay so. The
evidence is independent of any implementation: at $(p,n) = (3,9)$ the answer $s_F$ has

| $M$ | 4 | 12 | 20 | 28 | 40 |
| --- | --- | --- | --- | --- | --- |
| digits of $s_F$ | 2 | 18 | 32 | 47 | 67 |

digits, growing by roughly $1.8$ per facet, so the $\Theta(M)$ additions are performed on operands
that themselves lengthen with $M$. Carrying the instrumentation out to $M = 40$ in a separate run —
not comparable in absolute terms with the table above, which was timed differently — the cost **per
facet** rises by about $30\%$ between $M = 12$ and $M = 40$, the second-order term becoming visible
once the operands outgrow a machine word. The practical reading is that (4.1) is linear in the work
it schedules and mildly quadratic in the work it performs, and that the comparison against the
per-profile form is unaffected: that form pays the same widening arithmetic on top of a branching
factor that grows.

In $n$ the behaviour is different and worth recording, because it is easy to predict wrongly. The
state count does **not** track the number of partitions of $\le n$: at $(p,M) = (3,6)$ it runs
$65, 97, 119, 124, 125$ for $n = 8,\dots,16$, against $\sum_{k \le n} P(k) = 67, 139, 272, 508, 915$.
A partition needs enough facets to be reachable, so $M$ caps the state space, and the count
saturates. The practical consequence is that (4.1) is the right form when $M$ is the large
parameter; at small $M$ against large $n$, memoising the §3 form on (state, remaining profile tail)
can win instead.

---

## 7. Relation to the cycle-index formula

The standard route to $s_F$ is a Burnside average over $S_n$ acting on incidence tableaux, in which
the central quantity is

$$N(\lambda, p) \;=\; [x^p] \prod_{k \ge 1} \bigl(1 + x^k\bigr)^{m_k(\lambda)},$$

the number of sub-collections of the cycles of a permutation of cycle type $\lambda$ whose lengths
sum to $p$; the count is assembled from $N(\lambda,p)^{(M)} = N(N-1)\cdots(N-M+1)$ summed over cycle
types with the usual weights. The two derivations are independent, but they meet at a single
quantity, and the meeting point explains their different behaviour in $M$.

**Proposition 6.** *Split the $t = 0$ branches of (4.1) by whether they split a type class. The
non-splitting branches are the all-or-nothing selections, sub-collections of the **parts** of
$\lambda$ summing to $p$, and their number is*
$$[x^p]\prod_k (1 + x^k)^{m_k(\lambda)} \;=\; N(\lambda,p),$$
*the same function, evaluated on a partition of the vertex count instead of a cycle type.
Consequently the $t = 0$ term of (4.1) may be written*
$$\sum_{a\ \mathrm{splitting}} V\bigl(\mathrm{grow}(\lambda,a,0),\, j-1\bigr) \;+\; \bigl(N(\lambda,p) - (M-j)\bigr)\,V(\lambda,\, j-1).$$

The bracket $N(\lambda,p) - (M-j)$ reads as *how many all-or-nothing extensions are genuinely new
facets*. So the two routes handle facet-distinctness at the same quantity from opposite directions:
the cycle-index route raises $N$ to a falling factorial in one stroke, closing over all $M$ facets at
once, while (4.1) decrements it one facet at a time. That is also the sharpest statement of why the
cycle-index route is essentially free in $M$ — it never iterates over facets — and why this one
cannot be: the ordering is what buys the elementary derivation, and $M$ is the length of the
ordering.

---

## 8. Verification

The recursion was checked against an independent implementation of the cycle-index formula on
**180 parameter sets** spanning $p \le 4$, $M \le 6$, $n \le 10$, with zero mismatches, reaching
$s_F(4,6,10) = 20\,612\,880$. It reproduces the degenerate conventions ($M = 0$ counted as the empty
complex at $n = 0$; zero unless $p \le n \le pM$; zero when $\binom{n}{p} < M$) and a published
reference table of 21 values across all three labelings.

The per-profile weights of §3 were checked separately, by enumerating isomorphism classes directly
and grouping them by profile — this is what establishes $W(3,1,1,0) = 22$ and the $6 + 22 + 15 = 43$
of §3.1, and what confirms that condition (4) of Proposition 2 rejects only empty profiles.

The negative intermediates of §4 were checked not to disturb the totals. Over $p \le 3$,
$2 \le M \le 6$, $n \le 9$ — 63 parameter sets — a negative $t = 0$ term arises in **35** of them,
and in every one the recursion still agrees with the profile organisation of §3, which never forms
such a term. Six of these were taken further and matched against a direct construction of
$X / S_n$ from the definition of §1: $s_F = 1, 15, 1, 43, 222, 252$ at $(p,M,n) = (2,3,3)$,
$(2,4,4)$, $(3,4,4)$, $(3,4,5)$, $(2,5,5)$ and $(3,5,5)$.

---

## 9. Implementation notes

Five things are easy to get wrong and none is caught by the total coming out right on a single
small case; the first four were all live in a first implementation.

1. **Add the new vertices of the facet being placed, not of the next one.** At stage $i$ the fresh
   class has $n_i$ members. Adding $n_{i+1}$ instead silently loses vertices — at $c = (3,2,0,0)$ it
   leaves three where $|V_2| = 5$ — and every later stage is then built on a short state. The error is
   invisible whenever $n_i = n_{i+1}$, which is why a profile like $(3,1,1,0)$ hides it.
2. **Subtract $i-1$, not $i$.** There are $i-1$ earlier facets at stage $i$, hence $M-1$ at the last.
3. **Remove the forbidden branches; do not subtract from each surviving term.** Subtracting $i-1$
   from every term of the sum gives $\sum_a W(a) - (i-1)|MS_p(\lambda)|$, which at stage 3 of the
   $\lambda = (2,2,1)$ node of §5.2 is $19 - 15 = 4$ against the true $15$. The two agree only at the leaves, where the
   completions are all $1$ — which is exactly why a terminal-only statement of the rule looks
   correct.
4. **The all-or-nothing branches *are* the forbidden ones**, not merely equinumerous with them. It
   is because they leave $\lambda$ unchanged (Proposition 3(3)) that the correction can be applied
   without the type identities at all. An implementation that identifies them by comparing successor
   partitions must compare against the *sorted* $\lambda$, since $\mathrm{grow}$ returns a sorted
   partition while the live type-class list need not be in sorted order.
5. **Do not assert that memo entries are non-negative.** $V$ is a count only at realisable states;
   elsewhere it is an intermediate that exists to be cancelled, and it may be negative. At
   $(p,M,n) = (3,4,5)$ the $t = 0$ term of $V\bigl((3),2\bigr)$ is $-1$, and whole entries go
   negative at $(3,4,4)$ and $(3,6,6)$. A defensive check for non-negativity, or an unsigned
   accumulator, fires — or wraps — on correct code.

For the state to be canonical, keep the type classes sorted by size so that $\lambda$ and the
selection vectors agree in indexing; branches should be grouped by successor partition before
weighting, since the weight depends only on the successor.
