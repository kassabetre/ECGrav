# Sampling facet-labeled pure complexes, exactly

> **What this is.** A specification of `RandomUniformFacetLabeledPureSimplicialComplex`, matching
> the implementation in `Kernel/Subpackages/PureComplexes.wl`. Companion to `Theory.md`, and to
> `UnlabeledCount.md` / `UnlabeledSampler.md`, which do the same for the fully unlabeled case.
> It records what is being sampled, the pair-sampling principle that makes the draw exact, the
> lemma that licenses the sequential facet draw, the two indexing tricks that make it fast, and
> how correctness is established — by exact identities rather than by goodness-of-fit.

Notation: $p$ purity (every facet has $p$ vertices), $M$ facet order (number of facets), $n$
vertex count. $\lambda \vdash n$ is a cycle type, written either as a partition or as its
multiplicity vector $(m_1,\dots,m_p)$, with $z_\lambda = \prod_k k^{m_k} m_k!$, so $n!/z_\lambda$
permutations of $[n]$ have that type. $s_F(p,M,n)$ is `NumFacetLabeledPureComplexes[p,M,n]`, and
$x^{(r)} = x(x-1)\cdots(x-r+1)$ is the falling factorial.

---

## 1. What is being sampled

A **facet-labeled pure complex** is an ordered $M$-tuple of pairwise distinct $p$-subsets of a
vertex set, covering it, **up to relabelling the vertices**. The facets carry labels
$1,\dots,M$; the vertices do not. Formally, let

$$X \;=\; \{(F_1,\dots,F_M) : F_i \subseteq [n],\ |F_i| = p,\ F_i \neq F_j,\ \textstyle\bigcup_i F_i = [n]\},$$

on which $S_n$ acts by relabelling. The objects are the **orbits**, counted by $s_F(p,M,n)$.

Uniform means each of the $s_F(p,M,n)$ orbits has probability exactly $1/s_F(p,M,n)$ — **not**
that each element of $X$ is equally likely. The two differ whenever orbits have different sizes,
i.e. whenever some complex has a nontrivial automorphism.

The output is a facet list on vertex set $[n]$: a representative of its orbit, whose labelling is
itself uniform among the $n!/|\mathrm{Aut}|$ labellings of that class.

### The incidence-tableau picture

Transposing the incidence matrix gives the picture the counter uses, and the one the DP
alternative in §11 is stated in. A complex becomes a **multiset of $n$ nonempty subsets of
$[M]$** — one row per vertex, listing the facets it lies in — in which every label occurs in
exactly $p$ rows. The rows form a multiset because the vertices are unlabelled; "the $M$ facets
are pairwise distinct" becomes "the $M$ columns are pairwise distinct", the **separating**
condition. This is the form in which $s_F$ is derived (see the header comment at
`PureComplexes.wl:1782`), and it is a complete invariant: two facet-labeled complexes are
isomorphic iff their row multisets agree. That is what the test suite canonicalises with.

---

## 2. Why the vertex-labeled sampler cannot be ported

`RandomVertexLabeledPureSimplicialComplex` reads its weights straight off the counting recursion,
because with labelled vertices "how many complexes put $k$ new vertices in the next facet" is a
plain product of binomials, and the count decomposes facet by facet.

Nothing like that exists here. With the vertices unlabelled, $s_F$ is a **Burnside average over
cycle types**, not a sequential decomposition — there is no recursion to walk facet by facet, and
looking for one is wasted effort.

The older `RandomFacetLabeledPureSimplicialComplex` (no `Uniform` in the name) does exactly that
anyway, by iterative vertex addition off a shrinking pool. It is **not uniform**, which its
source comment says and its usage message does not. It is a different function and is not
specified here.

---

## 3. The principle: pair sampling

Let

$$P \;=\; \{(\sigma, x) : \sigma \in S_n,\ x \in X,\ \sigma \cdot x = x\}$$

be the set of **fixed pairs**. For any orbit $O$, the orbit–stabiliser theorem gives

$$\#\{(\sigma,x) \in P : x \in O\} \;=\; \sum_{x \in O} |\mathrm{Stab}(x)| \;=\; |O| \cdot \frac{n!}{|O|} \;=\; n!,$$

**the same number for every orbit, whatever its size.** So drawing $(\sigma,x)$ uniformly from
$P$ and discarding $\sigma$ leaves the orbit exactly uniform. No rejection, and no
automorphism-order weighting to get wrong.

This is the double-counting that proves the Cauchy–Frobenius (Burnside) lemma [1,2], read as a
sampling recipe rather than a counting one; summing the display over orbits gives
$|P| = \sum_\sigma |\mathrm{Fix}(\sigma)| = n!\,s_F(p,M,n)$, which is the identity §9 checks. The
same observation underlies Jerrum's Burnside process [3,4] — the difference is that the Burnside
process is a Markov chain that converges to the orbit-uniform distribution, whereas here $P$ is
simple enough to sample from **directly and exactly**, in one shot.

The identical argument makes `RandomUniformUnlabeledPureSimplicialComplex` exact; only the inner
object changes, from an ordered tuple to a $\sigma$-invariant set (`UnlabeledSampler.md`).

---

## 4. Which $\sigma$ contribute

$\sigma$ fixes a tuple $x$ exactly when every facet $F_i$ is $\sigma$-invariant, i.e. is a union
of cycles of $\sigma$. Calling the cycles **blocks**, a facet is a set of blocks whose sizes sum
to $p$.

Two consequences:

- **Only cycle types with all parts $\le p$ occur.** A cycle longer than $p$ fits inside no facet,
  so it could never be covered, and $|\mathrm{Fix}(\sigma)| = 0$. This is why the type list is
  `IntegerPartitions[n, All, Range[p]]` and not all of $\mathcal{P}(n)$.
- **Only the cycle *supports* matter, never the cyclic order.** Given the block partition, the
  number of $\sigma$ realising it is $\prod_k ((k-1)!)^{m_k}$, a constant depending only on
  $\lambda$. So sampling $\sigma$ uniformly among its type is the same as sampling the block
  partition uniformly, and the implementation only ever materialises blocks.

Writing $N(\nu)$ for the number of admissible facets available to a block collection with $\nu_k$
blocks of size $k$,

$$N(\nu) \;=\; [x^p] \prod_{k \ge 1} (1 + x^k)^{\nu_k},$$

the number of sub-multisets of the block sizes summing to $p$ — a truncated cycle-index-style
generating function [5,6]. `RandFLPCWeightCounts` builds the whole coefficient list out to $x^p$
by truncated convolution, so the degree-$n$ polynomial is never formed.

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
exactly $z_\lambda$ times, hence uniformly.

Step 3 is the recursive method for random generation [9,10]: at each stage choose a branch with
probability proportional to the number of completions below it, so that every complete object is
reached along exactly one path with probability $1/(\text{total count})$.

---

## 6. Completion counts

The quantity step 3 needs is: having spent $j$ subsets and with a set $U$ of blocks still
uncovered, how many ways are there to choose $M - j$ further weight-$p$ subsets, pairwise
distinct, distinct from the $j$ already spent, covering all of $U$?

Inclusion–exclusion over which uncovered blocks are *left* uncovered gives

$$C(\lambda, M, j, \nu) \;=\; \sum_{\omega \le \nu} (-1)^{|\omega|} \left(\prod_k \binom{\nu_k}{\omega_k}\right) \Big(N(\lambda - \omega) - j\Big)^{(M-j)},$$

where $\nu$ is the multiplicity vector of the sizes of the blocks in $U$, $\omega$ runs
componentwise over $0 \le \omega_k \le \nu_k$, and $|\omega| = \sum_k \omega_k$
(`RandFLPCCompletions`). The falling factorial counts ordered choices of distinct subsets, which
is what "the facets are labelled and pairwise distinct" asks for.

### 6.1 The lemma that makes step 3 exact

The $-j$ in $\big(N(\lambda-\omega) - j\big)$ is the whole subtlety, and it is exact for a reason
worth stating:

> **The number of completions depends on the draws so far only through $(j, U)$ — how *many*
> subsets have been used and which blocks remain uncovered — never on *which* subsets were used.**

Every used subset lies inside the already-covered blocks, so it automatically avoids any
$\omega \subseteq U$, and therefore drops out of each inclusion–exclusion term by itself; only its
*count* enters. Without this the state would have to remember the used set and the draw would not
decompose. It is the same observation that lets the vertex-labeled generator subtract $(j-1)$ for
a repeated facet.

Note this lemma is specific to the facet-labeled case. Its analogue **fails** in the unlabeled
sampler, where a naive port produced plausible complexes from a subtly wrong distribution
(`UnlabeledSampler.md` §7).

### 6.2 The $n^p$ indexing

The state is indexed by the **multiplicity vector of uncovered block sizes**, not by the uncovered
block *set*. Blocks of equal size are interchangeable at this point, so nothing is lost, and a
$2^{\#\text{blocks}}$ inclusion–exclusion becomes $\prod_k(\nu_k + 1) = O(n^p)$ terms — the same
polynomial scaling, and the same near-independence of $M$, that
`NumFacetLabeledPureComplexes` has.

### 6.3 Step-1 weights

$$w(\lambda) \;=\; \frac{n!}{z_\lambda}\, C(\lambda, M, 0, \lambda),$$

the number of permutations of type $\lambda$ times the number of tuples each fixes
(`RandFLPCTypeWeights`). By §3 these must satisfy

$$\sum_{\lambda} w(\lambda) \;=\; n!\; s_F(p,M,n). \tag{I1}$$

---

## 7. The grouped draw

Steps 3's naive form weighs every candidate subset separately: at $\{p,M,n\} = \{4,6,15\}$ that is
441 candidates $\times$ 6 facets, and it was ~100% of the cost of a draw.

But a candidate's weight enters **only** through $\nu\big(U \setminus S\big)$, a multiplicity
vector of block sizes with $\sum_k k\,\nu_k \le p$. The number of distinct such vectors is at most

$$\sum_{i=0}^{p} P(i) \;=\; 4,\ 7,\ 12 \quad\text{for } p = 2,3,4,$$

**however many hundreds of candidates there are** ($P$ the partition function). At $\{4,6,15\}$
the 441 candidates carry 3 distinct keys.

So the implementation groups candidates by that key and draws

$$\Pr[\text{group } g] \;\propto\; |g| \cdot C(\lambda, M, j{+}1, \text{key}(g)),$$

then a member of $g$ uniformly. Members of a group have equal weight, so this is the same
distribution — but with a handful of memo lookups per facet instead of one per candidate.

Two supporting details:

- **The key is carried as an integer.** Each $\nu_k \le p$, so base $p+1$ packs the vector without
  carrying: $\mathrm{code} = \sum_k \nu_k (p+1)^{k-1}$, and `Tally` on a packed integer array does
  the grouping in one pass.
- **Codes are updated in place, not recomputed.** When a facet covers block $b$ for the first
  time, exactly the candidates containing $b$ lose it from their intersection with $U$, so
  `codes[[containing[[b]]]] -= radix[[sizes[[b]]]]`. The blocks-to-candidates incidence, the
  candidate list, and the initial codes depend only on the cycle type, so `RandFLPCCandTable`
  memoises them; there are only $P_{\le p}(n)$ cycle types.
- Spent candidates are marked with the code $-1$, which no real key can take.

A wrong key here would still emit perfectly plausible complexes from the wrong distribution, so
the test suite checks the two invariants directly (§10).

---

## 8. Vertex count free

`RandomUniformFacetLabeledPureSimplicialComplex[{p,M}, k]` draws $n$ from the facet-labeled counts
themselves, $\Pr[n] \propto s_F(p,M,n)$ over $p \le n \le pM$, then samples at that $n$. The
previous implementation drew $n$ from the **vertex-labeled** counts and let a rejection step
repair the difference; weighting by $s_F$ is the distribution that was wanted all along, and needs
no repair.

The emptiness guard tests $s_F(p,M,n) \neq 0$, not the vertex-labeled count. The two vanish
together on $p \le n \le pM$, but $s_F$ also needs $\binom{n}{p} \ge M$ distinct facets, so it is
the correct test.

---

## 9. Why it is uniform

Three independent legs, in decreasing order of how much they would catch:

1. **(I1)** $\sum_\lambda w(\lambda) = n!\,s_F(p,M,n)$, exactly, as integers. A wrong cycle-type
   weight still produces plausible-looking complexes, so this identity — not the outputs — is the
   thing to check.
2. **The completion-count invariants** of §7: the candidate table lists exactly the admissible
   facets, and the incrementally maintained codes still decode to $\nu(U \setminus S)$ after any
   sequence of facets has been laid down.
3. **Goodness-of-fit** over exhaustively enumerated classes, as a backstop. At 60,000–100,000
   draws on $\{2,4,5\}$, $\{3,3,5\}$ and $\{4,3,8\}$ the draw is uniform and statistically
   indistinguishable from an independently derived exact sampler.

Chi-square is the weak leg. It is insensitive to exactly the errors that matter, and repeated runs
throw low p-values by chance; prefer legs 1 and 2.

---

## 10. Implementation map

All in `Kernel/Subpackages/PureComplexes.wl`; the derivation is the header comment at line 3785.
Private helpers are named `RandFLPC*`.

| line | symbol | role |
| --- | --- | --- |
| 3832 | `RandFLPCWeightCounts[ν]` | coefficients of $\prod_k (1+x^k)^{\nu_k}$ out to $x^p$; last entry is $N(\nu)$ (§4) |
| 3848 | `RandFLPCCompletions[λ,M,j,ν]` | $C(\lambda,M,j,\nu)$, the inclusion–exclusion of §6 |
| 3863 | `RandFLPCTypeWeights[p,M,n]` | step-1 weights $w(\lambda)$, with the types and their multiplicity vectors |
| 3892, 3895 | `$RandomFacetLabeledParallelThreshold`, `RandFLPCGoParallel` | sample-count gate, also gated on `$KernelCount` |
| 3899 | `RandFLPCCandTable[sizes,p]` | per-cycle-type candidate list, packed keys, place values, block→candidate incidence (§7) |
| 3921 | `RandFLPCOne[p,M,n]` | one sample, steps 1–4 |
| 3972 | primary pattern | `[{p,M,n}, numSamples]` |
| 4016 | overload | `[{p,M}, numSamples]` (§8) |
| 4059, 4068 | overloads | single-sample forms |

The counter it sits on is `NumFacetLabeledPureComplexes` at line 1878, with its own derivation at
line 1782 and helpers `NumFLPC*`.

All memo tables are released by the shared `ECGrav`Private`NumPCClearCache[]` — which must list
every memoised helper, `RandFLPCCandTable` included.

---

## 11. Cost

Measured warm in a fresh kernel, against the previous body extracted from `39dc9d3`. "Before" is
the per-candidate weighing; "after" is the grouped draw of §7.

| $\{p,M,n\}$ | before | after | speedup |
| --- | --- | --- | --- |
| $\{2,3,4\}$ | 0.094 ms | 0.066 ms | 1.4× |
| $\{2,4,6\}$ | 0.189 ms | 0.089 ms | 2.1× |
| $\{3,4,8\}$ | 0.553 ms | 0.099 ms | 5.6× |
| $\{2,8,12\}$ | 1.202 ms | 0.163 ms | 7.4× |
| $\{3,6,12\}$ | 2.688 ms | 0.153 ms | 17.6× |
| $\{3,7,14\}$ | 4.819 ms | 0.190 ms | 25.4× |
| $\{4,6,15\}$ | 14.42 ms | 0.280 ms | 51.5× |
| $\{3,10,20\}$ | 19.72 ms | 0.328 ms | 60.1× |
| $\{2,20,30\}$ | 18.97 ms | 0.516 ms | 36.8× |
| $\{3,12,24\}$ | 46.53 ms | 0.606 ms | 76.8× |
| $\{4,8,20\}$ | 52.80 ms | 0.513 ms | **102.9×** |
| $\{3,15,30\}$ | 108.5 ms | 1.240 ms | 87.5× |

The first call at a given $(p,M,n)$ builds the memo tables; the figures are warm, which is the
regime any real use is in.

**`$RandomFacetLabeledParallelThreshold = 2000`.** It was 100, which suited a draw costing 3–60 ms.
What the parallel branch pays is fixed rather than per-sample — `DistributedContexts` ships the
whole of `ECGrav`Private` , memo tables included, on every call — so cutting the per-draw cost by
20–100× moved the crossover out with it, and at 100 samples parallel had become a 3–16×
*pessimisation*. Remeasured on 11 subkernels with kernels up and tables warm, the crossover runs
200 at $\{2,3,4\}$, 300 at $\{2,4,6\}$, 800 at $\{3,6,12\}$, 1700 at $\{4,6,15\}$, 2800 at
$\{3,10,20\}$ and past 3200 at $\{4,8,20\}$. No single constant is right everywhere and the
penalties are lopsided — too low costs up to 16×, too high costs at most the parallel speedup,
only 2–3× anywhere near the crossover — so the constant sits at the top of the range. At 2000 the
serial/parallel ratios run 0.87–2.34; by 5000, 1.49–3.48.

### The alternative that was measured and rejected

A rejection-free **dynamic program over row-type multiplicities** carrying the current
column-equality partition $\pi$ as part of its state — conditioning on separation instead of
testing for it — is correct (its counts agree with $s_F$ on all 234 triples with $p \le 4$,
$M \le 6$) but is exponential where this sampler is polynomial. Its state space is
$2^M (n+1) (p+1)^M B_M$, and since $n \le pM$ the only free axes are $p$ and $M$, both exponential.
Building the table costs 9.7 s / 164 MB at $\{3,6,12\}$, 118 s / 1.8 GB at $\{3,7,14\}$, and did
not finish in 180 s / 4 GB at $\{3,8,16\}$ — where this sampler draws in 0.27 ms. The reason it
loses as a *counter* is structural: the $\big(N\big)^{(M)}$ falling factorial in
`NumFacetLabeledPureComplexes` already performs the Möbius inversion over the partition lattice
[11,12] in closed form, which is exactly what the $\pi$ state re-derives by brute force.

---

## 12. Tests

In `Tests/PureComplexes.wlt`:

| TestID | what it pins |
| --- | --- |
| `RandomUniformFacetLabeled-burnside-weights` | **(I1)**, exactly, on 10 parameter sets |
| `RandFLPCCandTable-keys-track-uncovered` | the §7 invariants, over every cycle type with $p \le 4$, $n \le 11$ |
| `RandomUniformFacetLabeled-uniform-over-classes` | chi-square over the 7 classes of $\{3,3,5\}$, pinned to the serial branch |
| `RandomUniformFacetLabeled-vertex-count-marginal` | the `{p,M}` form draws $n$ from $s_F$ (§8) |
| `RandomUniformSamplers-structure`, `-parallel-branch`, `-empty-space` | shape, the parallel path, and the guards of §8 |

The counter has its own seven, including a brute-force enumeration of the defining description and
the Stirling identity $B(p,n,M) = \sum_k S(M,k)\, F(p,n,k)$, which exercises the derivation rather
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
   Ch. 1–2. — The falling factorial and the Stirling identity $x^n = \sum_k S(n,k)\,x^{(k)}$ that
   imposes "pairwise distinct facets" (§6), and inclusion–exclusion.
3. M. Jerrum, "Uniform sampling modulo a group of symmetries using Markov chain simulation",
   *DIMACS Ser. Discrete Math. Theoret. Comput. Sci.* **10** (1993), 37–47. — The Burnside
   process: the same fixed-pair set $P$, sampled by Markov chain where §3 samples it directly.
4. L. A. Goldberg and M. Jerrum, "The 'Burnside process' converges slowly", *Combin. Probab.
   Comput.* **11** (2002), 21–34. — Why the chain is not a shortcut, and direct sampling from $P$
   is worth the effort when $P$ is tractable.
5. G. Pólya, "Kombinatorische Anzahlbestimmungen für Gruppen, Graphen und chemische
   Verbindungen", *Acta Math.* **68** (1937), 145–254; English translation in G. Pólya and
   R. C. Read, *Combinatorial Enumeration of Groups, Graphs, and Chemical Compounds*,
   Springer, 1987. See also Stanley, *Enumerative Combinatorics*, Vol. 2, §7.24, for cycle
   indicators in modern notation. — Cycle-index enumeration; the generating function
   $\prod_k (1+x^k)^{m_k}$ of §4 is a truncated instance.
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
    lattice, with $\mu(\hat 0, \pi) = \prod_B (-1)^{|B|-1}(|B|-1)!$; the closed form that the
    falling-factorial weight encodes (§11).
12. R. P. Stanley, *Enumerative Combinatorics*, Vol. 1, 2nd ed., Ch. 3. — The partition lattice
    $\Pi_M$, its Möbius function, and the Bell-number growth that sinks the DP alternative.

Companion specifications in this repository: `Theory.md` (the model),
`UnlabeledCount.md` and `UnlabeledSampler.md` (the fully unlabeled case, where the §6.1 lemma
fails and a different scheme is needed).
