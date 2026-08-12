(* ::Package:: *)

(* Regression tests for the PureComplexes-domain public functions
   (simplicial-complex combinatorics, graph observables, dimensions,
   geometric predicates). Expected values were verified against the current
   source; many are independent oracles (known surface topology, notebook
   author's stated outputs). Run with the Wolfram TestReport tool. *)

With[{local = If[StringQ[$InputFileName] && $InputFileName =!= "",
        FileNameJoin[{DirectoryName[$InputFileName], "TestPrelude.wl"}], ""]},
    Get[If[local =!= "" && FileExistsQ[local], local,
        "/Users/012759760/Desktop/Research/ECGravMathematicaPackage/Paclet/ECGrav/Tests/TestPrelude.wl"]]];

(* ---------- Basic constructions ---------- *)

VerificationTest[
    ECGrav`FacetIncidenceMatrix[fst1],
    {{1, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 1, 0}, {0, 1, 0, 1, 0, 0, 1}},
    TestID -> "FacetIncidenceMatrix-fst1"
];

VerificationTest[
    ECGrav`FacetAdjacencyMatrix[fst1],
    {{3, 0, 1}, {0, 3, 1}, {1, 1, 3}},
    TestID -> "FacetAdjacencyMatrix-fst1"
];

VerificationTest[
    ECGrav`CliquesFromFacetIncidence[ECGrav`FacetIncidenceMatrix[fst1]],
    fst1,
    TestID -> "CliquesFromFacetIncidence-roundtrip"
];

VerificationTest[
    ECGrav`FacetLabeledVertexListFromComplex[fst1],
    {{1, 3}, {2, 3}, {1}, {1}, {2}, {2}, {3}},
    TestID -> "FacetLabeledVertexListFromComplex-fst1"
];

VerificationTest[ECGrav`PureComplexQ[fst1], True,  TestID -> "PureComplexQ-pure"];
VerificationTest[ECGrav`PureComplexQ[{{1, 2, 3}, {4, 5}}], False, TestID -> "PureComplexQ-mixed"];
VerificationTest[ECGrav`FacetOrder[fst1], 3, TestID -> "FacetOrder-fst1"];
VerificationTest[ECGrav`CliqueOrder[ECGrav`GraphFromCliques[fst1]], 3, TestID -> "CliqueOrder-fst1"];

(* ---------- Euler characteristic (known topology) ---------- *)

VerificationTest[ECGrav`EulerChi[tetrahedron], 2, TestID -> "EulerChi-tetrahedron-S2"];
VerificationTest[ECGrav`EulerChi[octahedron], 2, TestID -> "EulerChi-octahedron-S2"];
VerificationTest[ECGrav`EulerChi[torus], 0, TestID -> "EulerChi-torus"];
VerificationTest[ECGrav`EulerChi[kleinbottle], 0, TestID -> "EulerChi-kleinbottle"];
VerificationTest[ECGrav`EulerChi[Table[0, {5}, {5}]], 5, TestID -> "EulerChi-5-isolated-vertices"];

(* ---------- Counting / combinatorics ---------- *)

(* NumVertexLabeledPureComplexes[p,q,n] counts q-element sets of p-subsets of [n] whose union is all of
   [n]. The values below were checked by direct enumeration of exactly that
   (Select[Subsets[Subsets[Range[n],{p}],{q}], Union@@#===Range[n]&]), so they are oracles
   independent of the recursion. NumVertexLabeledPureComplexes[4,3,12]=5775 is independently 12!/(4!^3 3!),
   the number of ways to partition 12 labelled vertices into 3 blocks of 4. *)

VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 3], 2649, TestID -> "NumVertexLabeledPureComplexes-3-3"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 3, 7], 945, TestID -> "NumVertexLabeledPureComplexes-3-3-7"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 3, 9], 280, TestID -> "NumVertexLabeledPureComplexes-3-3-9"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[2, 4, 6], 330, TestID -> "NumVertexLabeledPureComplexes-2-4-6"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 4, 8], 72380, TestID -> "NumVertexLabeledPureComplexes-3-4-8"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[2, 4], 900, TestID -> "NumVertexLabeledPureComplexes-2-4"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 4], 441061, TestID -> "NumVertexLabeledPureComplexes-3-4"];

(* n outside p <= n <= p q is answered without building a row; both ends must give 0. *)
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 2, 2], 0, TestID -> "NumVertexLabeledPureComplexes-n-below-p"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 2, 7], 0, TestID -> "NumVertexLabeledPureComplexes-n-above-pq"];

(* q<=0: the empty complex covers no vertices, so N(p,0,n) is 1 at n==0 and 0 otherwise, and
   summing it over n gives 1. Negative q counts nothing. These returned Indeterminate (3-arg)
   and an unevaluated Total[Table[..]] (2-arg) before. *)
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 0, 0], 1, TestID -> "NumVertexLabeledPureComplexes-q0-n0"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 0, 5], 0, TestID -> "NumVertexLabeledPureComplexes-q0-n-nonzero"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, -1, 5], 0, TestID -> "NumVertexLabeledPureComplexes-q-negative"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[3, 0], 1, TestID -> "NumVertexLabeledPureComplexes-2arg-q0"];
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[2, -1], 0, TestID -> "NumVertexLabeledPureComplexes-2arg-q-negative"];

(* Degenerate purity: no n in p..p q carries enough distinct p-subsets, so the sum is empty. *)
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[0, 3], 0, TestID -> "NumVertexLabeledPureComplexes-2arg-p0"];

(* The 2-arg overload must agree with summing the 3-arg form over every vertex count. *)
VerificationTest[
    ECGrav`NumVertexLabeledPureComplexes[3, 5],
    Total[Table[ECGrav`NumVertexLabeledPureComplexes[3, 5, n], {n, 0, 15}]],
    TestID -> "NumVertexLabeledPureComplexes-2arg-equals-row-total"
];

(* Rows are cached per (p,q) behind a high-water mark, so a descending sweep of q must return
   the same values as an ascending one. *)
VerificationTest[
    Table[ECGrav`NumVertexLabeledPureComplexes[4, q, 12], {q, 9, 1, -1}],
    Reverse[{0, 0, 5775, 57345750, 29179710945, 5551975972620, 627962951769615,
        51166525553400015, 3312982055194610645}],
    TestID -> "NumVertexLabeledPureComplexes-descending-q"
];

(* Rows are filled by an explicit bottom-up loop, so recursion depth no longer grows with q;
   the earlier recursive form died on TerminatedEvaluation[RecursionLimit] before q=700.
   p=1 keeps the entries tiny: q singleton facets cover exactly q vertices. *)
VerificationTest[ECGrav`NumVertexLabeledPureComplexes[1, 800, 800], 1, TestID -> "NumVertexLabeledPureComplexes-deep-q"];

(* The cache cap and its reset are private by design -- they tune an implementation detail, not
   the mathematics -- so the tests below have to name them in ECGrav`Private`. That is
   deliberate, hence the suppression; it is scoped to this block so a genuine private-context
   slip elsewhere in the file still shows up. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::PrivateContextSymbol:: *)

(* Rows are cached only up to ECGrav`Private`$NumPCMaxCachedQ; past it the tail is advanced
   transiently and discarded. Where that boundary falls must not change any value, so pin the
   same row against literals with the cap below, straddling, and far above the q range. The
   cache is deliberately NOT cleared between these, so a stale high-water mark left by one cap
   would surface in the next. *)
With[{expected = {0, 0, 0, 0, 8408400, 12234151930, 3513295545760, 465021334918140,
        38280314838325560, 2271754355448413836, 105805316724776028360,
        4080228574475590322786}},
    VerificationTest[
        Block[{ECGrav`Private`$NumPCMaxCachedQ = 4},
            Table[ECGrav`NumVertexLabeledPureComplexes[3, q, 14], {q, 1, 12}]],
        expected,
        TestID -> "NumVertexLabeledPureComplexes-cap-below-range"
    ];
    VerificationTest[
        Block[{ECGrav`Private`$NumPCMaxCachedQ = 7},
            Table[ECGrav`NumVertexLabeledPureComplexes[3, q, 14], {q, 12, 1, -1}]],
        Reverse[expected],
        TestID -> "NumVertexLabeledPureComplexes-cap-straddling-descending"
    ];
    VerificationTest[
        Block[{ECGrav`Private`$NumPCMaxCachedQ = 10000},
            Table[ECGrav`NumVertexLabeledPureComplexes[3, q, 14], {q, 1, 12}]],
        expected,
        TestID -> "NumVertexLabeledPureComplexes-cap-above-range"
    ];
];

(* A cap of 0 switches caching off altogether and must still give the same answer. *)
VerificationTest[
    Block[{ECGrav`Private`$NumPCMaxCachedQ = 0}, ECGrav`NumVertexLabeledPureComplexes[3, 8, 14]],
    465021334918140,
    TestID -> "NumVertexLabeledPureComplexes-caching-disabled"
];

(* NumPCClearCache[] reclaims the rows and returns the bytes freed; nothing else holds a row,
   so the next call must rebuild to the same values. *)
VerificationTest[
    (ECGrav`NumVertexLabeledPureComplexes[4, 60, 120]; ECGrav`Private`NumPCClearCache[] > 0),
    True,
    TestID -> "NumVertexLabeledPureComplexes-clear-cache-reclaims"
];
VerificationTest[
    (ECGrav`Private`NumPCClearCache[]; ECGrav`Private`NumPCClearCache[]),
    0,
    TestID -> "NumVertexLabeledPureComplexes-clear-cache-idempotent"
];
VerificationTest[
    (ECGrav`Private`NumPCClearCache[];
     {ECGrav`NumVertexLabeledPureComplexes[3, 4, 8], ECGrav`NumVertexLabeledPureComplexes[3, 3]}),
    {72380, 2649},
    TestID -> "NumVertexLabeledPureComplexes-correct-after-clear"
];

(* :!CodeAnalysis::EndBlock:: *)

(* ---- NumPureComplexes: former name, kept as an alias ---- *)

(* The rename to NumVertexLabeledPureComplexes did not touch the argument order, so the alias
   must agree with the new name everywhere, not merely on well-formed input. *)
VerificationTest[
    Module[{grid3, grid2},
        grid3 = Flatten[Table[{p, M, n}, {p, -1, 4}, {M, -1, 6}, {n, -1, 14}], 2];
        grid2 = Flatten[Table[{p, M}, {p, -1, 4}, {M, -1, 6}], 1];
        {Length[grid3], Length[grid2],
         AllTrue[grid3, Quiet[ECGrav`NumPureComplexes @@ #] === Quiet[ECGrav`NumVertexLabeledPureComplexes @@ #] &],
         AllTrue[grid2, Quiet[ECGrav`NumPureComplexes @@ #] === Quiet[ECGrav`NumVertexLabeledPureComplexes @@ #] &],
         ECGrav`NumPureComplexes[3, 3, 7]}],
    {768, 48, True, True, 945},
    TestID -> "NumPureComplexes-alias-agrees"
];

(* ---- NumFacetLabeledPureComplexes ---- *)

(* M labelled, pairwise distinct p-subsets of an n-set covering it, counted up to permutation of
   the unlabelled vertices. Values below are the separating-tableaux counts from the author's
   Facet-labeled-count.nb, independently reproduced here by direct enumeration of exactly that
   description (all M-tuples of distinct p-subsets covering [n], canonicalised over the n!
   relabellings). Note the argument order: (purity, facet order, vertex count), matching
   NumVertexLabeledPureComplexes, so these are the notebook's sF[p,n,M] with the last two swapped. *)
VerificationTest[
    ECGrav`NumFacetLabeledPureComplexes @@@
        {{2, 3, 4}, {2, 3, 5}, {3, 3, 5}, {3, 3, 6}, {2, 4, 4}, {2, 4, 5},
         {3, 4, 6}, {3, 4, 7}, {2, 5, 5}, {3, 5, 5}, {3, 5, 6}},
    {4, 3, 7, 10, 15, 29, 154, 207, 222, 252, 2472},
    TestID -> "NumFacetLabeledPureComplexes-known-values"
];

(* In-suite independent oracle: enumerate and canonicalise, rather than trusting the recursion. *)
VerificationTest[
    Module[{brute},
        brute = Function[{p, M, n},
            Module[{subs = Subsets[Range[n], {p}], tuples, canon},
                tuples = Select[Catenate[Permutations[#, {M}] & /@ Subsets[subs, {M}]],
                    Union @@ # === Range[n] &];
                canon = Function[t, First[Sort[Table[Sort /@ (t /. Thread[Range[n] -> perm]),
                    {perm, Permutations[Range[n]]}]]]];
                Length[Union[canon /@ tuples]]]];
        (brute @@@ #) === (ECGrav`NumFacetLabeledPureComplexes @@@ #) &@
            {{2, 3, 4}, {3, 3, 5}, {2, 4, 4}, {3, 2, 6}, {1, 3, 3}}],
    True,
    TestID -> "NumFacetLabeledPureComplexes-brute-force"
];

(* Merging equal columns of an incidence tableau gives B(p,n,M) = Sum_k StirlingS2[M,k] F(p,n,k),
   tying the separating count to the count of ALL (p,n,M) tableaux. Checking that identity against
   independently brute-forced B values exercises the whole derivation, not just the endpoints. *)
VerificationTest[
    Module[{allTab},
        allTab = Function[{p, n, M},
            Module[{rows = Rest[Subsets[Range[M]]]},
                Length[Select[Subsets[Range[Length[rows] + n - 1], {n}],
                    With[{ms = # - Range[0, n - 1]},
                        AllTrue[Range[M], Function[lab,
                            Count[rows[[ms]], r_ /; MemberQ[r, lab]] === p]]] &]]]];
        AllTrue[{{2, 4, 3}, {3, 5, 3}, {2, 4, 4}, {3, 6, 3}}, Function[s,
            With[{p = s[[1]], n = s[[2]], M = s[[3]]},
                allTab[p, n, M] ===
                    Sum[StirlingS2[M, k]*ECGrav`NumFacetLabeledPureComplexes[p, k, n], {k, 1, M}]]]]],
    True,
    TestID -> "NumFacetLabeledPureComplexes-stirling-identity"
];

(* Cost is essentially independent of M, which enters only as an exponent; this 92-digit value
   is the notebook's sF[3,10,50] and returns in milliseconds. *)
VerificationTest[
    ECGrav`NumFacetLabeledPureComplexes[3, 50, 10],
    153895449539108833097171275367657047043654150073334211607720849958357461588026327040000000000,
    TestID -> "NumFacetLabeledPureComplexes-large-M"
];

(* Zero outside p <= n <= p M, zero when there are too few distinct p-subsets to supply M facets,
   and M = 0 is the empty complex, as for NumVertexLabeledPureComplexes. *)
VerificationTest[
    ECGrav`NumFacetLabeledPureComplexes @@@
        {{3, 0, 0}, {3, 0, 5}, {3, -1, 5}, {3, 4, 2}, {3, 4, 100}, {3, 4, 3}, {-1, 3, 4}},
    {1, 0, 0, 0, 0, 0, 0},
    TestID -> "NumFacetLabeledPureComplexes-degenerate"
];

(* The two-argument form sums the three-argument one over the vertex count. *)
VerificationTest[
    Table[{ECGrav`NumFacetLabeledPureComplexes[p, M],
           Total[Table[ECGrav`NumFacetLabeledPureComplexes[p, M, n], {n, 0, p*M}]]},
        {p, 2, 3}, {M, 2, 5}],
    Table[{#, #} &@Total[Table[ECGrav`NumFacetLabeledPureComplexes[p, M, n], {n, 0, p*M}]],
        {p, 2, 3}, {M, 2, 5}],
    TestID -> "NumFacetLabeledPureComplexes-2arg-equals-sum"
];

VerificationTest[
    Quiet[{ECGrav`NumFacetLabeledPureComplexes[1.5, 2],
           ECGrav`NumFacetLabeledPureComplexes[1, 2, 3, 4]}],
    {$Failed, $Failed},
    TestID -> "NumFacetLabeledPureComplexes-argerr"
];

(* ---- NumUnlabeledPureComplexes ---- *)

(* M pairwise distinct p-subsets of an n-set covering it, up to relabelling the vertices AND
   permuting the facets -- the isomorphism classes. Values below are the UNLABELED column of the
   brute-force oracle table, produced by explicitly enumerating the covering M-sets and grouping
   them into S_n orbits; the same enumeration reproduced the vertex- and facet-labeled counts in
   the two tests above, which is what ties the three together. *)
VerificationTest[
    ECGrav`NumUnlabeledPureComplexes @@@
        {{2, 2, 3}, {2, 3, 4}, {2, 3, 5}, {2, 4, 4}, {2, 4, 5}, {2, 4, 6}, {2, 5, 4},
         {2, 5, 5}, {2, 5, 6}, {2, 6, 4}, {2, 6, 5}, {2, 6, 6},
         {3, 2, 6}, {3, 3, 5}, {3, 3, 6}, {3, 4, 5}, {3, 4, 6}, {3, 4, 7},
         {4, 2, 6}, {4, 3, 6}, {1, 3, 3}, {1, 4, 4}},
    {1, 2, 1, 2, 4, 3, 1, 5, 9, 1, 5, 15, 1, 3, 3, 5, 15, 17, 1, 4, 1, 1},
    TestID -> "NumUnlabeledPureComplexes-known-values"
];

(* In-suite independent oracle: enumerate the covering M-sets and canonicalise each over the n!
   relabellings, rather than trusting the cycle-index derivation. *)
VerificationTest[
    Module[{brute},
        brute = Function[{p, M, n},
            Module[{subs = Subsets[Range[n], {p}], sets, perms},
                sets = Select[Subsets[subs, {M}], Union @@ # === Range[n] &];
                perms = Permutations[Range[n]];
                Length[Union[Table[
                    First[Sort[Table[Sort[Map[Sort[pi[[#]]] &, s]], {pi, perms}]]],
                    {s, sets}]]]]];
        (brute @@@ #) === (ECGrav`NumUnlabeledPureComplexes @@@ #) &@
            {{2, 3, 4}, {2, 4, 5}, {3, 3, 5}, {3, 4, 5}, {2, 5, 5}, {3, 2, 6}, {1, 3, 3}}],
    True,
    TestID -> "NumUnlabeledPureComplexes-brute-force"
];

(* The derivation, not just its endpoints. Burnside over S_n says the class count is the average
   over sigma of the number of sigma-invariant covering M-sets; here those fixed sets are counted
   by explicit enumeration instead of through the orbit polynomial Product (1+z^d)^n_d and the
   isolated-vertex differencing, so a wrong orbit-size count would show up even where the totals
   happen to agree. *)
VerificationTest[
    Module[{burnside},
        burnside = Function[{p, M, n},
            Module[{subs = Subsets[Range[n], {p}], sets, perms, act},
                sets = Select[Subsets[subs, {M}], Union @@ # === Range[n] &];
                perms = Permutations[Range[n]];
                act = Function[{s, pi}, Sort[Map[Sort[pi[[#]]] &, s]]];
                Total[Table[Count[sets, s_ /; act[s, pi] === s], {pi, perms}]]/n!]];
        (burnside @@@ #) === (ECGrav`NumUnlabeledPureComplexes @@@ #) &@
            {{2, 3, 4}, {2, 4, 5}, {3, 3, 5}, {3, 4, 5}, {2, 5, 5}}],
    True,
    TestID -> "NumUnlabeledPureComplexes-burnside-identity"
];

(* The Newton recurrence that produces |Fix(sigma)|, pinned at its own level rather than only
   through the totals. For every cycle type of S_n it must agree with a direct count of the
   M-element sets of p-subsets that a permutation of that type leaves invariant. This is the
   level the totals are averaged from, and a wrong l(k) there still yields plausible-looking
   totals, so it is checked separately. The first returned value is the number of (case, cycle
   type) pairs actually compared, so the test cannot pass by comparing nothing. *)
VerificationTest[
    Module[{permOfType, bruteFix, cases, pairs = 0},
        permOfType = Function[lambda,
            Module[{off = 0, res = {}},
                Do[res = Join[res, RotateLeft[Range[off + 1, off + k]]]; off += k, {k, lambda}];
                res]];
        bruteFix = Function[{lambda, p, M},
            Module[{n = Total[lambda], pi, subs},
                pi = permOfType[lambda];
                subs = Subsets[Range[n], {p}];
                Count[Subsets[subs, {M}], s_ /; Sort[Map[Sort[pi[[#]]] &, s]] === s]]];
        cases = Select[Flatten[Table[{p, M, n}, {p, 1, 3}, {M, 0, 4}, {n, 1, 6}], 2],
            Binomial[Binomial[#[[3]], #[[1]]], #[[2]]] <= 6000 &];
        {Length[cases],
         AllTrue[cases, Function[c,
            With[{p = c[[1]], M = c[[2]], n = c[[3]]},
                AllTrue[IntegerPartitions[n], Function[lambda,
                    pairs++;
                    ECGrav`Private`NumULPCFixSets[Tally[lambda], p, M] ===
                        bruteFix[lambda, p, M]]]]]],
         pairs}],
    {90, True, 435},
    TestID -> "NumUnlabeledPureComplexes-fix-recurrence"
];

(* External cross-check at p = 2, where these are the unlabelled graphs on n nodes with no
   isolated vertex. Unlabelled graphs on n nodes number 1, 1, 2, 4, 11, 34, 156, 1044, and every
   such graph is one with no isolated vertex on j <= n nodes, so differencing that sequence gives
   the no-isolated-vertex counts -- which must be the sum over M of the count here. This is an
   oracle from outside the package entirely. *)
VerificationTest[
    Module[{gAll = {1, 1, 2, 4, 11, 34, 156, 1044}},
        Table[Total[Table[ECGrav`NumUnlabeledPureComplexes[2, M, n], {M, 0, Binomial[n, 2]}]],
            {n, 0, 7}] ===
        Table[gAll[[n + 1]] - If[n >= 1, gAll[[n]], 0], {n, 0, 7}]],
    True,
    TestID -> "NumUnlabeledPureComplexes-unlabeled-graphs"
];

(* Each isomorphism class carries at least one vertex-labeled and at least one facet-labeled
   complex, and distinct classes carry disjoint sets of them, so the unlabeled count can never
   exceed either. Cheap, and it catches a count that is right on the table but wrong off it. *)
VerificationTest[
    Module[{grid = Flatten[Table[{p, M, n}, {p, 1, 4}, {M, 0, 6}, {n, 0, 10}], 2]},
        {Count[grid, s_ /; (ECGrav`NumUnlabeledPureComplexes @@ s) > 0],
         AllTrue[grid, ECGrav`NumUnlabeledPureComplexes @@ # <=
             Min[ECGrav`NumVertexLabeledPureComplexes @@ #,
                 ECGrav`NumFacetLabeledPureComplexes @@ #] &]}],
    {93, True},
    TestID -> "NumUnlabeledPureComplexes-dominated-by-labelled"
];

(* Same degenerate conventions as the two labelled counts, which is checked here by comparing
   against them rather than by restating the values. n = 100 is in the list on purpose: the
   guards have to reject it before the cycle-type sum is entered, since reaching
   IntegerPartitions[100] would not return. *)
VerificationTest[
    Module[{deg = {{3, 0, 0}, {3, 0, 5}, {3, -1, 5}, {3, 4, 2}, {3, 4, 100}, {3, 4, 3},
                   {-1, 3, 4}, {0, 1, 0}, {0, 1, 1}, {0, 2, 0}, {0, 0, 0}, {2, 1, 2},
                   {2, 3, -1}, {2, 1, 3}}},
        {ECGrav`NumUnlabeledPureComplexes @@@ deg,
         AllTrue[deg, ECGrav`NumUnlabeledPureComplexes @@ # ===
             ECGrav`NumFacetLabeledPureComplexes @@ # ===
             ECGrav`NumVertexLabeledPureComplexes @@ # &]}],
    {{1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0}, True},
    TestID -> "NumUnlabeledPureComplexes-degenerate"
];

(* The two-argument form telescopes the differencing instead of summing, so this checks the
   telescoping as well as the convention. *)
VerificationTest[
    Table[{ECGrav`NumUnlabeledPureComplexes[p, M],
           Total[Table[ECGrav`NumUnlabeledPureComplexes[p, M, n], {n, 0, p*M}]]},
        {p, 1, 3}, {M, 0, 5}],
    Table[{#, #} &@Total[Table[ECGrav`NumUnlabeledPureComplexes[p, M, n], {n, 0, p*M}]],
        {p, 1, 3}, {M, 0, 5}],
    TestID -> "NumUnlabeledPureComplexes-2arg-equals-sum"
];

(* The memoised cycle-type sums go behind the shared cache, so clearing must reclaim them and
   leave the values unchanged. *)
VerificationTest[
    Module[{before = ECGrav`NumUnlabeledPureComplexes[3, 6, 12], reclaimed},
        reclaimed = ECGrav`Private`NumPCClearCache[];
        {before, reclaimed > 0, ECGrav`NumUnlabeledPureComplexes[3, 6, 12] === before}],
    {241, True, True},
    TestID -> "NumUnlabeledPureComplexes-clear-cache"
];

VerificationTest[
    Quiet[{ECGrav`NumUnlabeledPureComplexes[1.5, 2],
           ECGrav`NumUnlabeledPureComplexes[1, 2, 3, 4]}],
    {$Failed, $Failed},
    TestID -> "NumUnlabeledPureComplexes-argerr"
];

VerificationTest[ECGrav`RankComb[{0, 1, 2}, 5], 0, TestID -> "RankComb-first"];
VerificationTest[
    ECGrav`UnrankComb[ECGrav`RankComb[{0, 1, 2}, 5], 5, 3],
    {0, 1, 2},
    TestID -> "RankComb-UnrankComb-roundtrip"
];

(* ---------- Isomorphism / automorphism ---------- *)

VerificationTest[
    ECGrav`IsomorphicSimplicialComplexQ[{{1, 2, 3}, {3, 5, 6}, {3, 5, 7}}, {{4, 9, 10}}],
    False,
    TestID -> "IsomorphicSimplicialComplexQ-different"
];
VerificationTest[ECGrav`SimplicialComplexAutomorphismGroupOrder[fst1], 8,  TestID -> "AutomorphismOrder-fst1"];
VerificationTest[ECGrav`SimplicialComplexAutomorphismGroupOrder[fst2], 24, TestID -> "AutomorphismOrder-fst2"];

(* ---------- Purity / connectivity observables ---------- *)

VerificationTest[ECGrav`RmsPurity[fst1], 0, TestID -> "RmsPurity-pure-is-0"];
VerificationTest[ECGrav`FractionInLargestComponent[fst1], 1, TestID -> "FractionInLargestComponent-fst1"];
VerificationTest[ECGrav`FacetDeg[fst1, 2], 2, TestID -> "FacetDeg-fst1"];
VerificationTest[ECGrav`HyperDeg[fst1, {2, 4}], 1, TestID -> "HyperDeg-fst1-edge"];
VerificationTest[ECGrav`Branchedness[octaG], 0, TestID -> "Branchedness-octahedron-graph"];

(* Facet-list input must match the clique-graph form and must not leak a Null
   into the result (regression for the AdjacencyGraph-misinterpretation bug). *)
VerificationTest[
    ECGrav`Branchedness[octahedron],
    ECGrav`Branchedness[octaG],
    TestID -> "Branchedness-octahedron-facetlist-matches-graph"
];
VerificationTest[FreeQ[ECGrav`Branchedness[fst1], Null], True, TestID -> "Branchedness-facetlist-no-Null-leak"];
VerificationTest[ECGrav`Branchedness[fst1], 0, TestID -> "Branchedness-fst1-facetlist"];
(* Adjacency-matrix input must still route to the graph computation. *)
VerificationTest[ECGrav`Branchedness[K4], 0, TestID -> "Branchedness-K4-adjacency-matrix"];

(* ---------- Spheres, balls, links, stars ---------- *)

VerificationTest[
    Sort[VertexList[ECGrav`Sph[octaG, 1]]],
    {2, 3, 4, 5},
    TestID -> "Sph-octahedron-vertex1"
];
VerificationTest[
    Sort[Sort /@ ECGrav`Lnk[octahedron, {1}]],
    {{2, 3}, {2, 4}, {3, 5}, {4, 5}},
    TestID -> "Lnk-octahedron-vertex1"
];

(* ---------- Dimensions, holes, f-vector, volume ---------- *)

VerificationTest[ECGrav`AvgKDim[octaG], 2, TestID -> "AvgKDim-octahedron"];
VerificationTest[ECGrav`CountHoles[octaG], {1, 0, 1}, TestID -> "CountHoles-octahedron-Betti-S2"];
VerificationTest[ECGrav`FVector[octahedron], {6, 12, 8}, TestID -> "FVector-octahedron"];
VerificationTest[ECGrav`GVolume[octaG, 1, 1], 4, TestID -> "GVolume-octahedron"];
VerificationTest[ECGrav`Deg[K4, 1], 3, TestID -> "Deg-K4-vertex1"];
VerificationTest[ECGrav`AvgDeg[CompleteGraph[4]], 3, TestID -> "AvgDeg-K4"];

(* ---------- Geometric predicates on surfaces (independent oracles) ---------- *)

VerificationTest[
    ECGrav`OrientableCombinatorialManifoldQ /@ surfaces,
    {False, False, False, False, True, True, True},
    TestID -> "OrientableCombinatorialManifoldQ-surfaces"
];
VerificationTest[
    ECGrav`CombinatorialManifoldQ /@ {tetrahedron, octahedron, torus},
    {True, True, True},
    TestID -> "CombinatorialManifoldQ-surfaces"
];
VerificationTest[ECGrav`CombinatorialSphereQ[tetrahedron], True, TestID -> "CombinatorialSphereQ-tetrahedron"];
VerificationTest[ECGrav`CombinatorialSphereQ[octahedron], True, TestID -> "CombinatorialSphereQ-octahedron"];

(* CombinatorialSphereQ tests genuine sphere-ness: of the seven surfaces, only the two
   2-spheres (tetrahedron, octahedron) qualify. The closed non-spheres (Klein bottle,
   projective plane, torus) and the two bounded Mobius strips are all False. This is the
   regression for the old CombinatorialSphereQ[torus] == True bug (the old code tested
   closed-manifold-ness, now split out into ClosedCombinatorialManifoldQ below). *)
VerificationTest[
    ECGrav`CombinatorialSphereQ /@ surfaces,
    {False, False, False, False, False, True, True},
    TestID -> "CombinatorialSphereQ-surfaces"
];
VerificationTest[ECGrav`CombinatorialSphereQ[torus], False, TestID -> "CombinatorialSphereQ-torus-not-sphere"];

(* ClosedCombinatorialManifoldQ tests combinatorial-manifold-without-boundary-ness (the
   behaviour the old CombinatorialSphereQ actually implemented): the five closed surfaces
   are True; the two Mobius strips, which have boundary, are False. *)
VerificationTest[
    ECGrav`ClosedCombinatorialManifoldQ /@ surfaces,
    {False, False, True, True, True, True, True},
    TestID -> "ClosedCombinatorialManifoldQ-surfaces"
];
VerificationTest[ECGrav`ClosedCombinatorialManifoldQ[torus], True, TestID -> "ClosedCombinatorialManifoldQ-torus-is-closed-manifold"];
VerificationTest[ECGrav`ClosedCombinatorialManifoldQ[{{1, 2, 3}}], False, TestID -> "ClosedCombinatorialManifoldQ-single-simplex-has-boundary"];

VerificationTest[ECGrav`DSphereQ[octaG], True, TestID -> "DSphereQ-octahedron"];

(* DSphereQ now tests genuine sphere-ness of the graph's clique complex: the octahedron graph
   is a 2-sphere (True), but the torus graph's clique complex is a closed manifold that is not
   a sphere, so it is False. ClosedDGraphQ is the closed-combinatorial-manifold test the old
   DSphereQ actually implemented (torus -> True). A K4 is a solid tetrahedron (a 3-ball), which
   has boundary, so ClosedDGraphQ is False for it (whereas DGraphQ, allowing boundary, is True). *)
VerificationTest[ECGrav`DSphereQ[torusG], False, TestID -> "DSphereQ-torus-not-sphere"];
VerificationTest[ECGrav`ClosedDGraphQ[octaG], True, TestID -> "ClosedDGraphQ-octahedron-sphere"];
VerificationTest[ECGrav`ClosedDGraphQ[torusG], True, TestID -> "ClosedDGraphQ-torus-is-closed-manifold"];
VerificationTest[ECGrav`ClosedDGraphQ[tetraG], False, TestID -> "ClosedDGraphQ-K4-has-boundary"];

VerificationTest[ECGrav`DGraphQ[tetraG], True, TestID -> "DGraphQ-tetrahedron"];

(* ---------- Random complex generators (structural invariants) ---------- *)

(* Each generator must return a well-formed pure simplicial complex of the requested shape,
   whatever the random draw. Seeded only so a failure reproduces; the assertions themselves
   are seed-independent structural invariants.
   NB: the three uniform samplers used to be kept out of the suite on the grounds that a first
   call in a fresh kernel cost about 30 s of parallel distribution. That cost was the bug: their
   ParallelTable branches distributed only ECGrav`Private`, so the subkernels could not evaluate
   NumVertexLabeledPureComplexes and the whole draw came back unevaluated. With that fixed, a first parallel
   call is 0.12 s once TestPrelude has launched the kernels and the serial path costs nothing
   special, so all three are covered here. *)

(* ---- RandomVertexLabeledPureSimplicialComplex ---- *)

(* A sample is q distinct sorted p-subsets whose union is exactly [n]. {2,6,4} is included
   because it is the degenerate end: six edges on four vertices is the whole of K4, so every
   facet after the fourth brings no new vertex and goes through the repeat-rejection branch. *)
VerificationTest[
    Module[{specs}, SeedRandom[310];
        specs = {{2, 3, 4}, {3, 4, 8}, {2, 5, 7}, {4, 3, 9}, {1, 5, 5}, {2, 6, 4}};
        AllTrue[specs, Function[sp,
            Module[{p = sp[[1]], q = sp[[2]], n = sp[[3]], s},
                s = ECGrav`RandomVertexLabeledPureSimplicialComplex[sp, 300];
                Length[s] === 300 && AllTrue[s, Function[c,
                    Length[c] === q && Length[DeleteDuplicates[c]] === q &&
                    AllTrue[c, Function[f, Length[f] === p && f === Sort[f]]] &&
                    Union @@ c === Range[n]]]]]]],
    True,
    TestID -> "RandomVertexLabeledPureSimplicialComplex-structure"
];

(* The parallel threshold is private -- it tunes an implementation detail, not the mathematics --
   so the two tests that pin a branch have to name it. Deliberate, hence the suppression, scoped
   so a genuine private-context slip elsewhere still shows up. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::PrivateContextSymbol:: *)

(* Above ECGrav`Private`$RandomComplexParallelThreshold the draw goes to ParallelTable. That
   branch used to distribute only ECGrav`Private`, so the subkernels had no NumVertexLabeledPureComplexes,
   the composition weights never evaluated, and every sample came back as facet lists with
   unevaluated RandomSample[...] expressions inside -- silently, and only above the threshold.
   The threshold is forced low here so this keeps exercising the parallel path whatever the
   shipped default becomes, and $KernelCount is asserted so that a kernel-less environment
   fails loudly instead of quietly running the serial path and passing for the wrong reason. *)
VerificationTest[
    Block[{ECGrav`Private`$RandomComplexParallelThreshold = 100},
        Module[{s}, SeedRandom[311];
            s = ECGrav`RandomVertexLabeledPureSimplicialComplex[{3, 2, 6}, 3000];
            {$KernelCount > 0, Length[s], AllTrue[s, Function[c,
                Length[c] === 2 && Length[DeleteDuplicates[c]] === 2 &&
                AllTrue[c, Function[f, Length[f] === 3 && AllTrue[f, IntegerQ]]] &&
                Union @@ c === Range[6]]]}]],
    {True, 3000, True},
    TestID -> "RandomVertexLabeledPureSimplicialComplex-parallel-branch"
];

(* The generator is meant to be UNIFORM over vertex-labeled pure complexes -- that is what the
   NumVertexLabeledPureComplexes-weighted composition draw buys. Enumerate the 16 complexes at p=2,q=3,n=4
   and chi-square the sample against uniform. BlockRandom keeps the seeding from disturbing the
   rest of the suite; the threshold is loose, so this fails only on a genuinely wrong
   distribution, and Total[counts] guards against a keying slip making it vacuous. Pinned to the
   serial path: 6400 samples is over the parallel threshold, and parallel draws split the random
   stream by kernel count, which would make the outcome differ from machine to machine rather
   than being fixed by the seed. The parallel branch's distribution is checked out of suite. *)
VerificationTest[
    BlockRandom[SeedRandom[312];
      Block[{ECGrav`Private`$RandomComplexParallelThreshold = Infinity},
        Module[{space, sample, counts, chi2},
            space = Sort /@ Select[Subsets[Subsets[Range[4], {2}], {3}],
                Union @@ # === Range[4] &];
            sample = ECGrav`RandomVertexLabeledPureSimplicialComplex[{2, 3, 4}, 6400];
            counts = Lookup[Counts[Sort /@ sample], space, 0];
            chi2 = Total[(counts - 400)^2]/400.;
            {Length[space], Total[counts],
             SurvivalFunction[ChiSquareDistribution[Length[space] - 1], chi2] > 0.001}]]],
    {16, 6400, True},
    TestID -> "RandomVertexLabeledPureSimplicialComplex-uniform"
];

(* :!CodeAnalysis::EndBlock:: *)

(* ---- the two rejection samplers built on top of it ---- *)

(* Both draw from RandomPureComplexFacets and accept with a weight that re-targets the
   distribution, so their output must still be well-formed complexes of the requested shape. *)
VerificationTest[
    Module[{gens, specs}, SeedRandom[320];
        gens = {ECGrav`RandomUniformUnlabeledPureSimplicialComplex,
                ECGrav`RandomUniformFacetLabeledPureSimplicialComplex};
        specs = {{2, 3, 4}, {2, 4, 5}, {3, 3, 5}};
        AllTrue[Tuples[{gens, specs}], Function[gs,
            Module[{g = gs[[1]], sp = gs[[2]], p, q, n, s},
                {p, q, n} = sp;
                s = g[sp, 20];
                Length[s] === 20 && AllTrue[s, Function[c,
                    Length[c] === q && Length[DeleteDuplicates[c]] === q &&
                    AllTrue[c, Function[f, Length[f] === p && f === Sort[f]]] &&
                    Union @@ c === Range[n]]]]]]],
    True,
    TestID -> "RandomUniformSamplers-structure"
];

(* Above numSamples = 20 both switch to ParallelTable. Both used to distribute only
   ECGrav`Private`, leaving the subkernels without NumVertexLabeledPureComplexes, so 100% of samples came
   back as facet lists holding unevaluated RandomSample[...] -- and unlike the vertex-labeled
   generator, whose threshold was 10^4, these fired at 21 samples. $KernelCount is asserted so a
   kernel-less environment fails loudly rather than passing on the serial path. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::PrivateContextSymbol:: *)
(* Both thresholds must be lowered, not just one: since 1.6.0 the unlabeled generator reads its
   own $RandomUnlabeledParallelThreshold (500, measured -- its memo tables are rebuilt per
   subkernel, so parallel pays much later than for the facet-labeled one). Forcing only the
   facet-labeled threshold would leave the unlabeled sampler on the serial path at 25 samples and
   this test would quietly stop exercising the branch it names. *)
VerificationTest[
    Block[{ECGrav`Private`$RandomFacetLabeledParallelThreshold = 10,
           ECGrav`Private`$RandomUnlabeledParallelThreshold = 10},
        Module[{u, f}, SeedRandom[321];
            u = ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{2, 3, 4}, 25];
            f = ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{2, 3, 4}, 25];
            {$KernelCount > 0,
             AllTrue[Join[u, f], Function[c,
                Length[c] === 3 && Length[DeleteDuplicates[c]] === 3 &&
                AllTrue[c, Function[e, Length[e] === 2 && AllTrue[e, IntegerQ]]] &&
                Union @@ c === Range[4]]]}]],
    {True, True},
    TestID -> "RandomUniformSamplers-parallel-branch"
];
(* :!CodeAnalysis::EndBlock:: *)

(* The point of the unlabeled sampler is uniformity over ISOMORPHISM classes, not over labeled
   complexes. At p=2,q=3,n=4 there are 16 labeled complexes but only two classes -- the path,
   with 12 labellings, and the star, with 4 -- so a uniform labeled draw would split them 75/25,
   and it takes the automorphism structure to bring them to 50/50. Drawn 20 at a time to stay on
   the serial branch, which keeps the seed reproducible. 400 samples is deliberately modest: the
   failure this guards against lands at 75/25 and gives chi-square about 100 on one degree of
   freedom, so it is caught with room to spare while keeping the draw fast. *)
VerificationTest[
    BlockRandom[SeedRandom[322];
        Module[{canon, sample, counts, chi2},
            canon = Function[c, First[Sort[Table[Sort[Sort /@ (c /. Thread[Range[4] -> perm])],
                {perm, Permutations[Range[4]]}]]]];
            sample = Join @@ Table[
                ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{2, 3, 4}, 20], {20}];
            counts = Lookup[Counts[canon /@ sample],
                {{{1, 2}, {1, 3}, {1, 4}}, {{1, 2}, {1, 3}, {2, 4}}}, 0];
            chi2 = Total[(counts - 200)^2]/200.;
            {Total[counts], SurvivalFunction[ChiSquareDistribution[1], chi2] > 0.001}]],
    {400, True},
    TestID -> "RandomUniformUnlabeledPureSimplicialComplex-uniform-over-classes"
];

(* ---- RandomUniformUnlabeledPureSimplicialComplex: Burnside pair sampling ---- *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::PrivateContextSymbol:: *)

(* The identity the whole scheme rests on: the step-1 cycle-type weights must sum to
   n! NumUnlabeledPureComplexes[p,M,n] EXACTLY. A wrong cycle-type weight still emits
   perfectly plausible complexes -- the distribution is simply not uniform -- so the failure
   would be silent in any test that only inspects output. Check the identity itself. *)
VerificationTest[
    AllTrue[Flatten[Table[{p, M, n}, {p, 1, 3}, {M, 1, 5}, {n, 1, 8}], 2], Function[s,
        Total[ECGrav`Private`RandULPCTypeWeights[s[[1]], s[[2]], s[[3]]][[2]]] ===
            s[[3]]!*ECGrav`NumUnlabeledPureComplexes @@ s]],
    True,
    TestID -> "RandomUniformUnlabeledPureSimplicialComplex-weight-identity"
];

(* Step 3 draws the size profile -- how many orbits of each size the set uses -- from its own
   marginal. Those marginals must sum back to the per-cycle-type covering count, or the profile
   is drawn from the wrong distribution while every sample still looks valid. *)
VerificationTest[
    AllTrue[Flatten[Table[{p, M, n}, {p, 1, 3}, {M, 1, 4}, {n, 1, 6}], 2], Function[s,
        AllTrue[IntegerPartitions[s[[3]]], Function[lam,
            Total[ECGrav`Private`RandULPCCCov[lam, s[[1]], s[[2]], #] & /@
                  ECGrav`Private`RandULPCProfiles[lam, s[[1]], s[[2]]]] ===
            ECGrav`Private`RandULPCFixCov[lam, s[[1]], s[[2]]]]]]],
    True,
    TestID -> "RandomUniformUnlabeledPureSimplicialComplex-profile-marginals"
];

(* At the identity cycle type the covering sets a permutation fixes ARE the vertex-labeled
   complexes, which pins the covering count against a function computed a completely different
   way. *)
VerificationTest[
    AllTrue[Flatten[Table[{p, M, n}, {p, 1, 3}, {M, 1, 5}, {n, 1, 7}], 2], Function[s,
        ECGrav`Private`RandULPCFixCov[ConstantArray[1, s[[3]]], s[[1]], s[[2]]] ===
            ECGrav`NumVertexLabeledPureComplexes @@ s]],
    True,
    TestID -> "RandomUniformUnlabeledPureSimplicialComplex-fixcov-at-identity"
];
(* :!CodeAnalysis::EndBlock:: *)

(* Every sample must be M pairwise distinct p-subsets covering [n] -- including at n = p M,
   where the facets must nearly partition the vertex set. That regime is the reason covering is
   imposed during the draw rather than by rejecting afterwards: rejection's expected trial count
   is A(p,M,n)/U(p,M,n), which is 11 at {2,4,8} and 66 at {3,4,12} and keeps climbing. *)
VerificationTest[
    BlockRandom[SeedRandom[9091];
        Module[{ok, sets},
            ok = Function[{p, M, n, s},
                Length[s] === M && Length[Union[s]] === M && Union @@ s === Range[n] &&
                AllTrue[s, Length[#] === p &]];
            sets = {{2, 3, 4}, {3, 4, 6}, {2, 6, 6}, {4, 3, 6}, {2, 5, 8}, {3, 4, 10},
                    {2, 4, 8}, {3, 4, 12}, {2, 6, 10}};
            {Length[sets],
             AllTrue[sets, Function[s,
                AllTrue[ECGrav`RandomUniformUnlabeledPureSimplicialComplex[s, 30],
                    ok[s[[1]], s[[2]], s[[3]], #] &]]]}]],
    {9, True},
    TestID -> "RandomUniformUnlabeledPureSimplicialComplex-valid-samples"
];

(* Uniformity over the enumerated isomorphism classes, at a size where a non-uniform draw would
   show. Drawn in one call on the serial branch so the seed is reproducible. *)
VerificationTest[
    BlockRandom[SeedRandom[5150];
        Module[{p = 3, M = 4, n = 6, canon, cls, counts, chi2, k = 900},
            canon = Function[c, First[Sort[Table[Sort[Sort /@ (c /. Thread[Range[n] -> perm])],
                {perm, Permutations[Range[n]]}]]]];
            cls = Union[canon /@ Select[Subsets[Subsets[Range[n], {p}], {M}],
                Union @@ # === Range[n] &]];
            counts = Lookup[Counts[canon /@
                ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{p, M, n}, k]], cls, 0];
            chi2 = Total[(counts - k/Length[cls])^2]/(k/Length[cls]);
            {Length[cls], Total[counts], FreeQ[counts, 0],
             SurvivalFunction[ChiSquareDistribution[Length[cls] - 1], chi2] > 0.001}]],
    {15, 900, True, True},
    TestID -> "RandomUniformUnlabeledPureSimplicialComplex-uniform-15-classes"
];

(* With the vertex count free it must be drawn from NumUnlabeledPureComplexes, not from the
   vertex-labeled counts -- that was the old body's bug, papered over by rejection. *)
VerificationTest[
    BlockRandom[SeedRandom[6161];
        Module[{p = 2, M = 4, k = 1500, ns, want, got},
            ns = Length[Union @@ #] & /@
                ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{p, M}, k];
            want = N[#/Total[#]] &@
                Table[ECGrav`NumUnlabeledPureComplexes[p, M, n], {n, p, p*M}];
            got = N[Table[Count[ns, n], {n, p, p*M}]/k];
            Max[Abs[want - got]] < 0.05]],
    True,
    TestID -> "RandomUniformUnlabeledPureSimplicialComplex-vertex-count-marginal"
];

(* RandomUniformFacetLabeledPureSimplicialComplex no longer rejects. It samples a fixed pair
   (sigma, x) of the S_n action instead: every orbit contributes exactly n! such pairs whatever
   its size, so dropping sigma leaves the orbit uniform. The whole scheme rests on the Burnside
   weights summing to n! NumFacetLabeledPureComplexes[p,M,n] -- if the cycle-type weights were
   wrong the draw would still look plausible, so check the identity itself, exactly, rather than
   only its consequences. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::PrivateContextSymbol:: *)
VerificationTest[
    AllTrue[{{2, 3, 4}, {2, 3, 5}, {3, 3, 5}, {2, 4, 4}, {2, 4, 5}, {3, 3, 6},
             {2, 4, 6}, {3, 4, 7}, {3, 4, 9}, {2, 5, 8}}, Function[s,
        With[{p = s[[1]], M = s[[2]], n = s[[3]]},
            Total[ECGrav`Private`RandFLPCTypeWeights[p, M, n][[3]]] ===
                n!*ECGrav`NumFacetLabeledPureComplexes[p, M, n]]]],
    True,
    TestID -> "RandomUniformFacetLabeled-burnside-weights"
];

(* Uniformity over the facet-labeled classes, the property the sampler exists to have. Enumerate
   the classes by canonicalising over the n! relabellings and chi-square the draw. Pinned to the
   serial branch so the seed fixes the outcome; the parallel branch splits the stream by kernel
   count, which would vary by machine. *)
VerificationTest[
    BlockRandom[SeedRandom[330];
      Block[{ECGrav`Private`$RandomFacetLabeledParallelThreshold = Infinity},
        Module[{p = 3, M = 3, n = 5, numClasses, canon, smp, counts, chi2},
            numClasses = ECGrav`NumFacetLabeledPureComplexes[p, M, n];
            canon = Function[c, First[Sort[Table[Sort /@ (c /. Thread[Range[n] -> perm]),
                {perm, Permutations[Range[n]]}]]]];
            smp = ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{p, M, n}, 1400];
            counts = Lookup[Counts[canon /@ smp], Union[canon /@ smp], 0];
            chi2 = Total[(counts - 1400/numClasses)^2]/(1400/numClasses);
            {numClasses, Length[counts], Total[counts],
             SurvivalFunction[ChiSquareDistribution[numClasses - 1], N[chi2]] > 0.001}]]],
    {7, 7, 1400, True},
    TestID -> "RandomUniformFacetLabeled-uniform-over-classes"
];

(* With the vertex count free, it is now drawn from NumFacetLabeledPureComplexes directly rather
   than from the vertex-labeled counts with rejection repairing the difference, so the marginal
   must follow those counts. *)
VerificationTest[
    BlockRandom[SeedRandom[331];
      Block[{ECGrav`Private`$RandomFacetLabeledParallelThreshold = Infinity},
        Module[{p = 2, M = 4, smp, ns, support, expected, obs, chi2},
            support = Select[Range[p, p*M], ECGrav`NumFacetLabeledPureComplexes[p, M, #] > 0 &];
            smp = ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{p, M}, 4000];
            ns = Max[Max /@ #] & /@ smp;
            obs = Lookup[Counts[ns], support, 0];
            expected = N[4000*#/Total[#]] &@ (ECGrav`NumFacetLabeledPureComplexes[p, M, #] & /@ support);
            chi2 = Total[(obs - expected)^2/expected];
            {Total[obs],
             SurvivalFunction[ChiSquareDistribution[Length[support] - 1], N[chi2]] > 0.001}]]],
    {4000, True},
    TestID -> "RandomUniformFacetLabeled-vertex-count-marginal"
];
(* :!CodeAnalysis::EndBlock:: *)

(* Same empty-sample-space guards as the vertex-labeled generator. *)
VerificationTest[
    Quiet[Join[
        ECGrav`RandomUniformUnlabeledPureSimplicialComplex[#, 2] & /@
            {{3, 4, 100}, {2, 3, 7}, {3, 4, 3}, {0, 3, 4}, {0, 3}},
        ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[#, 2] & /@
            {{3, 4, 100}, {2, 3, 7}, {3, 4, 3}, {0, 3, 4}, {3, 0}},
        {ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{2, 3, 4}, -1],
         ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{2, 3, 4}, -1]}]],
    ConstantArray[$Failed, 12],
    TestID -> "RandomUniformSamplers-empty-space"
];

(* An empty sample space must fail loudly. Every composition weight carries a
   NumVertexLabeledPureComplexes factor, so on an empty space they are all zero and the draw cannot proceed;
   the old code let that unevaluated result flow on and returned malformed "complexes". *)
VerificationTest[
    Quiet[ECGrav`RandomVertexLabeledPureSimplicialComplex[#, 2] & /@
        {{3, 4, 100}, {2, 3, 7}, {3, 4, 3}, {3, 4, 2}, {3, 4, 0}, {3, 4, -3},
         {0, 3, 4}, {3, 0, 5}, {3, -2, 5}, {0, 3}, {3, 0}}],
    ConstantArray[$Failed, 11],
    TestID -> "RandomVertexLabeledPureSimplicialComplex-empty-space"
];

(* The one-complex overloads must propagate that $Failed rather than First[] into it. *)
VerificationTest[
    Quiet[{ECGrav`RandomVertexLabeledPureSimplicialComplex[{2, 4, 3}],
           ECGrav`RandomVertexLabeledPureSimplicialComplex[{0, 3}],
           ECGrav`RandomVertexLabeledPureSimplicialComplex[{2, 3, 4}, -5],
           ECGrav`RandomVertexLabeledPureSimplicialComplex[{2, 3, 4}, 0]}],
    {$Failed, $Failed, $Failed, {}},
    TestID -> "RandomVertexLabeledPureSimplicialComplex-degenerate-args"
];

(* Compositions are always feasible: each prefix of j facets covering v_j vertices must have
   Binomial[v_j,p] >= j distinct p-subsets to draw from. The generator used to re-draw the whole
   composition while this failed; that loop was dead, because a k is only ever drawn with
   positive weight and a positive weight already implies the condition. This pins that down. *)
VerificationTest[
    Module[{s, comps}, SeedRandom[313];
        s = ECGrav`RandomVertexLabeledPureSimplicialComplex[{2, 6, 8}, 500];
        comps = Function[c, Module[{seen = {}},
            Table[With[{k = Length[Complement[f, seen]]}, seen = Union[seen, f]; k], {f, c}]]] /@ s;
        AllTrue[comps, Function[cp,
            Total[cp] === 8 && AllTrue[Range[2, 6], Binomial[Total[Take[cp, #]], 2] >= # &]]]],
    True,
    TestID -> "RandomVertexLabeledPureSimplicialComplex-composition-feasible"
];

VerificationTest[
    Module[{c}, SeedRandom[301];
        c = ECGrav`RandomFacetLabeledPureSimplicialComplex[{3, 4}];
        {ECGrav`PureComplexQ[c], Length[c], AllTrue[c, Length[#] == 3 &]}],
    {True, 4, True},
    TestID -> "RandomFacetLabeledPureSimplicialComplex-shape"
];

(* RandomUnlabeledPseudoManifold returns {complex, pastingSites}; its complex (part 1) is pure. *)
VerificationTest[
    Module[{r}, SeedRandom[304];
        r = ECGrav`RandomUnlabeledPseudoManifold[{3, 4}];
        {ECGrav`PureComplexQ[r[[1]]], Length[r[[1]]]}],
    {True, 4},
    TestID -> "RandomUnlabeledPseudoManifold-shape"
];

(* One MCMC sweep returns an association describing the current complex / weight / energy. *)
VerificationTest[
    Module[{r}, SeedRandom[305];
        r = ECGrav`RandomPureSimplicialComplexMCMCSweep[
            ECGrav`RandomFacetLabeledPureSimplicialComplex[{3, 4}], 0, 5];
        {Head[r], KeyExistsQ[r, "complex"], ECGrav`PureComplexQ[r["complex"]]}],
    {Association, True, True},
    TestID -> "RandomPureSimplicialComplexMCMCSweep-smoke"
];

(* ---- the MCMC driver's two stages: equilibriation and correlation time ---- *)

(* Unlike the graph MC drivers, these two are fully serial -- no ParallelTable anywhere in
   their six definitions -- so two runs under the same SeedRandom agree exactly and values
   can be pinned rather than only checked for internal consistency.

   Both stages Print a diagnostic plot, so the calls are wrapped in Block[{Print = Null &}]
   except where the plot itself is the thing under test. Note that this does NOT suppress
   PrintTemporary, which both use for progress; those vanish on their own. *)

(* Equilibriation reports where it stopped and whether it got there. outWinLength = 400 and
   inWinLength = 30, so a converged run reports eqlT = numsweeps - 370; small systems pass the
   criterion at the very first check, the sweep after the 400-long window fills, which pins
   eqlT to 31. The chain must also have moved off the seed. *)
VerificationTest[
    Module[{r}, SeedRandom[401];
        r = Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed24, 0]];
        {Keys[r], r["converged"], r["eqlT"],
         ECGrav`PureComplexQ[r["state"]["complex"]],
         Sort[Length /@ r["state"]["complex"]],
         SubsetQ[Range[8], Union @@ r["state"]["complex"]],
         r["state"]["vertexCount"] === Length[Union @@ r["state"]["complex"]],
         r["state"]["complex"] =!= Sort[Sort /@ mcmcSeed24]}],
    {{"eqlT", "converged", "state"}, True, 31, True, {2, 2, 2, 2}, True, True, True},
    TestID -> "RandomPureSimplicialComplexMCMCEquilibriate-converges"
];

(* ... and still emits its energy-vs-sweep diagnostic plot. *)
VerificationTest[
    Module[{n, res}, SeedRandom[406];
        {n, res} = emittedGraphicsAndResult[
            ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed24, 1]];
        {n >= 1, res["converged"], res["eqlT"]}],
    {True, True, 31},
    TestID -> "RandomPureSimplicialComplexMCMCEquilibriate-emits-plot"
];

(* Exhausting the sweep budget must be reported, not silently returned as a number that a
   successful run could also have produced. A budget below the 400-sweep comparison window
   never reaches the convergence criterion, so this failure is deterministic rather than a
   race with the chain. The reported eqlT is the honest lower bound
   maxNumSweeps - 400 + 30, clamped to stay positive; 100 sweeps exercises the clamp and
   380 the unclamped branch. ::noconv must carry the budget it was given and the eqlT it is
   about to return. *)
VerificationTest[
    Module[{r, msgs},
        msgs = messageArguments["RandomPureSimplicialComplexMCMCEquilibriate::noconv",
            Quiet@Block[{ECGrav`$ECGravMaxEquilibriationSweeps = 100}, SeedRandom[407];
                r = Block[{Print = Null &},
                    ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed24, 0]]]];
        {r["eqlT"], r["converged"], Length[msgs], msgs[[1, {1, 4}]]}],
    {1, False, 1, {100, 1}},
    TestID -> "RandomPureSimplicialComplexMCMCEquilibriate-budget-exhausted-clamped"
];

VerificationTest[
    Module[{r, msgs},
        msgs = messageArguments["RandomPureSimplicialComplexMCMCEquilibriate::noconv",
            Quiet@Block[{ECGrav`$ECGravMaxEquilibriationSweeps = 380}, SeedRandom[408];
                r = Block[{Print = Null &},
                    ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed24, 0]]]];
        {r["eqlT"], r["converged"], Length[msgs], msgs[[1, {1, 4}]]}],
    {10, False, 1, {380, 10}},
    TestID -> "RandomPureSimplicialComplexMCMCEquilibriate-budget-exhausted-lower-bound"
];

(* Negative control for the two tests above: a converged run must NOT issue ::noconv. Both
   of those run under Quiet, where MUnit's own unexpected-message check cannot see anything,
   so without this they would pass just as well against a function that always reported
   failure -- and against a messageArguments that never returned an empty list. *)
VerificationTest[
    (SeedRandom[401];
     messageArguments["RandomPureSimplicialComplexMCMCEquilibriate::noconv",
        Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed24, 0]]]),
    {},
    TestID -> "RandomPureSimplicialComplexMCMCEquilibriate-converged-run-is-quiet"
];

(* The three-argument overload runs a different chain: it permutes vertices between two
   facets, so the support never moves. The returned complex must therefore still cover
   exactly the seed's vertex set, and the tracked observable is the edge count rather than
   the vertex count. *)
VerificationTest[
    Module[{r}, SeedRandom[409];
        r = Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed34, True, 1]];
        {r["converged"], r["eqlT"], Keys[r["state"]],
         Union @@ r["state"]["complex"],
         ECGrav`PureComplexQ[r["state"]["complex"]],
         Sort[Length /@ r["state"]["complex"]],
         r["state"]["edgeCount"] ===
            Length[Union @@ (Subsets[#, {2}] & /@ r["state"]["complex"])],
         r["state"]["complex"] =!= Sort[Sort /@ mcmcSeed34]}],
    {True, 31, {"complex", "edgeCount", "weight", "energy"}, Range[7], True, {3, 3, 3, 3},
     True, True},
    TestID -> "RandomPureSimplicialComplexMCMCEquilibriate-fixed-vertex-support"
];

(* HoldNumberOfVerticesFixed -> False hands straight over to the two-argument overload, so
   under the same seed the two calls agree exactly, down to the free-vertex chain's
   "vertexCount" key. *)
VerificationTest[
    Module[{a, b}, SeedRandom[403];
        a = Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed24, False, 1]];
        SeedRandom[403];
        b = Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[mcmcSeed24, 1]];
        {a === b, Keys[a["state"]]}],
    {True, {"complex", "vertexCount", "weight", "energy"}},
    TestID -> "RandomPureSimplicialComplexMCMCEquilibriate-free-vertex-delegates"
];

(* The correlation time echoes back the equilibriation time it was handed and reports one
   value per measured observable: the energy, the vertex (or edge) count, and each supplied
   operator. corrT is the maximum of those, and every one of them is floored at 2. *)
VerificationTest[
    Module[{r}, SeedRandom[506];
        r = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[mcmcSeed24, 12,
                {Function[c, Total[Flatten[c]]], Function[c, Max[Flatten[c]]]}, 2]];
        {Keys[r], r["eqlT"], Length[r["corrTValues"]],
         AllTrue[r["corrTValues"], IntegerQ[#] && # >= 2 &],
         r["corrT"] === Max[r["corrTValues"]],
         ECGrav`PureComplexQ[r["state"]["complex"]], Length[r["state"]["complex"]]}],
    {{"eqlT", "corrT", "corrTValues", "state"}, 12, 4, True, True, True, 4},
    TestID -> "RandomPureSimplicialComplexMCMCCorrelationTime-shape"
];

(* Differential test of the reported correlation times against an independent
   implementation of the documented integration rule (mcmcCorrT in the prelude).

   The driver measures its own observables internally, so the test recovers them: one of the
   supplied operators records the complex it is handed, which is exactly the post-sweep
   complex the driver measures in the same iteration, and the energy and vertex count are
   deterministic functions of that complex. The second operator is a sweep counter -- a
   monotone ramp whose autocorrelation never turns over, which is what forces the lag loop
   to run to its end instead of stopping after a handful of terms. Without it every value
   sits on the floor of 2 and the comparison is nearly vacuous, which is the failure mode
   this test was rebuilt to avoid: turnedOver below records that both branches of the
   integral -- terminating on a non-positive term, and running out of lags -- were taken.

   Vacuity guard: the same oracle on shuffled series, whose autocorrelation is destroyed,
   must give a different answer, so the agreement is not an artifact of every observable
   defaulting to 2. *)
VerificationTest[
    Module[{recorded = {}, sweepNo = 0, r, numsweeps, rows, oracle, shuffled},
        SeedRandom[504];
        r = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[mcmcSeed24, 31,
                {Function[c, AppendTo[recorded, c]; Total[Flatten[c]]],
                 Function[c, ++sweepNo]}, 1]];
        numsweeps = Length[recorded];
        (* the driver's four observable series, rebuilt from the recorded complexes; the
           energy is -Log of the facet-labeled weight, written exactly as the sweep writes
           it so the reconstruction is bit-for-bit rather than merely close *)
        rows = {-Log[(1.0/Binomial[8, Length[Union @@ #]])*
                    ECGrav`PureComplexFacetStabilizerGroupOrder[#]*(4!/Length[Union @@ #]!)] & /@ recorded,
                Length[Union @@ #] & /@ recorded,
                Total[Flatten[#]] & /@ recorded,
                Range[numsweeps]};
        oracle = mcmcCorrT[#, numsweeps - 10] & /@ rows;
        SeedRandom[9];
        shuffled = mcmcCorrT[RandomSample[#], numsweeps - 10][[1]] & /@ rows;
        {numsweeps, Last[recorded] === r["state"]["complex"],
         r["corrTValues"] === oracle[[All, 1]],
         r["corrT"] === Max[oracle[[All, 1]]],
         oracle[[All, 2]],                       (* turnover reached on each series? *)
         Max[r["corrTValues"]] > 2,              (* the lag loop went past the floor *)
         shuffled =!= oracle[[All, 1]]}],
    {310, True, True, True, {True, True, True, False}, True, True},
    TestID -> "RandomPureSimplicialComplexMCMCCorrelationTime-integration-rule"
];

(* The same differential test on a short run, which is where the end of the lag range stops
   being a rounding error and starts setting the answer.

   The integral runs to lastLag = numsweeps - 10, and for a ramp the normalized
   autocorrelation at lag t is ((N-t)/N)^2, so over 310 sweeps the last ten lags are worth
   about 0.02 in total and no plausible change to the lag budget can move Ceiling -- the test
   above cannot see one. Over 20 sweeps those same ten lags are half the integral. At this
   length the linear counter reports 7 where a budget of numsweeps - 20 would report 2, and
   the quadratic counter reports 6 where dropping the one-lag-short convention would report 7
   -- which is the only place either convention is visible in the output at all.

   Both counters are synthetic on purpose: a real observable turns over within a handful of
   lags and never reaches the end of the range, so it cannot probe it. *)
VerificationTest[
    Module[{recorded = {}, sweepNo = 0, sweepNoSq = 0, r, numsweeps, rows, oracle},
        SeedRandom[508];
        r = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[mcmcSeed24, 2,
                {Function[c, AppendTo[recorded, c]; Total[Flatten[c]]],
                 Function[c, ++sweepNo], Function[c, (++sweepNoSq)^2]}, 1]];
        numsweeps = Length[recorded];
        rows = {-Log[(1.0/Binomial[8, Length[Union @@ #]])*
                    ECGrav`PureComplexFacetStabilizerGroupOrder[#]*(4!/Length[Union @@ #]!)] & /@ recorded,
                Length[Union @@ #] & /@ recorded,
                Total[Flatten[#]] & /@ recorded,
                Range[numsweeps],
                Range[numsweeps]^2};
        oracle = mcmcCorrT[#, numsweeps - 10] & /@ rows;
        {numsweeps,
         AllTrue[rows, Length[DeleteDuplicates[#]] > 1 &],  (* nothing excluded as stuck *)
         r["corrTValues"],
         r["corrTValues"] === oracle[[All, 1]],
         r["corrT"] === Max[oracle[[All, 1]]],
         oracle[[All, 2]]}],
    {20, True, {2, 2, 2, 7, 6}, True, True, {True, True, True, False, False}},
    TestID -> "RandomPureSimplicialComplexMCMCCorrelationTime-short-run-lag-range"
];

(* Same differential test on the fixed-vertex overload, whose second observable is the edge
   count and whose chain must hold the support at the seed's vertex set for every one of the
   sweeps it measures. p = 3 here because for p = 2 the facets are the edges, so the edge
   count would be pinned at M and the observable would never fluctuate. *)
VerificationTest[
    Module[{recorded = {}, sweepNo = 0, r, numsweeps, rows, oracle, shuffled},
        SeedRandom[505];
        r = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[mcmcSeed34, True, 31,
                {Function[c, AppendTo[recorded, c]; Total[Flatten[c]]],
                 Function[c, ++sweepNo]}, 1]];
        numsweeps = Length[recorded];
        rows = {-Log[ECGrav`PureComplexFacetStabilizerGroupOrder[#]*(1.0*(4!)/(7!))] & /@ recorded,
                Length[Union @@ (Subsets[#, {2}] & /@ #)] & /@ recorded,
                Total[Flatten[#]] & /@ recorded,
                Range[numsweeps]};
        oracle = mcmcCorrT[#, numsweeps - 10] & /@ rows;
        SeedRandom[9];
        shuffled = mcmcCorrT[RandomSample[#], numsweeps - 10][[1]] & /@ rows;
        {AllTrue[recorded, Union @@ # === Range[7] &],
         r["state"]["edgeCount"] ===
            Length[Union @@ (Subsets[#, {2}] & /@ r["state"]["complex"])],
         Length[DeleteDuplicates[rows[[2]]]] > 1,   (* the edge count really did fluctuate *)
         r["corrTValues"] === oracle[[All, 1]],
         r["corrT"] === Max[oracle[[All, 1]]],
         oracle[[All, 2]],
         shuffled =!= oracle[[All, 1]]}],
    {True, True, True, True, True, {True, True, True, False}, True},
    TestID -> "RandomPureSimplicialComplexMCMCCorrelationTime-fixed-vertex-integration-rule"
];

(* As with the equilibriator, False delegates to the shorter overload. *)
VerificationTest[
    Module[{a, b}, SeedRandom[507];
        a = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[mcmcSeed24, False, 12,
                {Function[c, Total[Flatten[c]]]}, 1]];
        SeedRandom[507];
        b = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[mcmcSeed24, 12,
                {Function[c, Total[Flatten[c]]]}, 1]];
        {a === b, Keys[a["state"]]}],
    {True, {"complex", "vertexCount", "weight", "energy"}},
    TestID -> "RandomPureSimplicialComplexMCMCCorrelationTime-free-vertex-delegates"
];

(* An observable that never moves has no autocorrelation time, and the driver says so per
   observable rather than inventing a number. The triangle is the only 2-complex with three
   facets on three vertices, so the fixed-vertex chain cannot leave it: the energy (weight 1
   for labelingChoise 0), the edge count and a constant operator are all frozen, which is
   the one situation that reaches ::alldefault.

   The value each ::stuck message quotes is the check on a fixed bug: for an operator the
   message used to name operators[[i-2]] but quote observablesTable[[i-2, 1]], i.e. the
   energy's value, so the third message here would have read 0. rather than 7. *)
VerificationTest[
    Module[{r, stuck, allDefault},
        stuck = messageArguments["RandomPureSimplicialComplexMCMCCorrelationTime::stuck",
            allDefault = messageArguments[
                "RandomPureSimplicialComplexMCMCCorrelationTime::alldefault",
                Quiet@Block[{Print = Null &}, SeedRandom[503];
                    r = ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[
                        {{1, 2}, {1, 3}, {2, 3}}, True, 5, {Function[c, 7]}, 0]]]];
        {r["corrT"], r["corrTValues"], r["state"]["complex"],
         stuck[[All, 1]], stuck[[All, 2]], Length[allDefault], allDefault[[1, 2]]}],
    {2, {2, 2, 2}, {{1, 2}, {1, 3}, {2, 3}},
     {"the energy", "the magnetization", Function[c, 7]}, {0., 3, 7}, 1, {2, 2, 2}},
    TestID -> "RandomPureSimplicialComplexMCMCCorrelationTime-frozen-observables"
];

(* Negative control for the test above: a chain that does move must report no stuck
   observable at all. *)
VerificationTest[
    messageArguments["RandomPureSimplicialComplexMCMCCorrelationTime::stuck",
        Quiet@Block[{Print = Null &}, SeedRandom[502];
            ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[mcmcSeed24, 5,
                {Function[c, Length[Union @@ c]]}, 0]]],
    {},
    TestID -> "RandomPureSimplicialComplexMCMCCorrelationTime-moving-chain-is-not-stuck"
];

(* ---- RandomPureSimplicialComplexMCMC: the driver that composes the two stages ---- *)

(* Four definitions: from a seed complex, free-vertex and fixed-vertex, and two continuation
   overloads that take a previous run's association and skip straight to sampling. All four
   are serial and reproduce exactly under a seed. Each returns {measurements, data}, where a
   measurement row is {sweep index, energy, vertex or edge count, operator values..., complex}
   and data carries the two stages' output alongside the final state. *)

(* The rows must describe the states they record: recomputing the energy, the vertex count and
   the operators from each row's own complex has to give back that row's numbers. This is what
   a shape check cannot see -- a driver that Sowed a stale state, or whose operator columns
   slipped out of step, still returns a well-formed list of the right length. The corrupted
   copy is the negative control on the predicate itself. *)
VerificationTest[
    Module[{ops, meas, data}, SeedRandom[601];
        ops = {Function[c, Total[Flatten[c]]], Function[c, Max[Flatten[c]]]};
        {meas, data} = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[mcmcSeed24, ops, 1, 5]];
        {Length[meas], Length /@ meas, meas[[All, 1]],
         mcmcMeasurementsConsistentQ[meas, ops, 1],
         mcmcMeasurementsConsistentQ[ReplacePart[meas, {1, 2} -> meas[[1, 2]] + 1.0], ops, 1],
         Length[DeleteDuplicates[meas[[All, -1]]]] > 1,   (* the chain moved between samples *)
         Last[meas][[-1]] === data[["state", "complex"]],
         Keys[data], data["converged"], data["eqlT"], data["corrT"], data["corrTValues"]}],
    {5, {6, 6, 6, 6, 6}, Range[5], True, False, True, True,
     {"eqlT", "converged", "state", "corrT", "corrTValues"}, True, 31, 2, {2, 2, 2, 2}},
    TestID -> "RandomPureSimplicialComplexMCMC-measurements-are-self-consistent"
];

(* Samples drawn from a chain not shown to have reached its stationary distribution would not
   be the uniform draws this function documents, so an equilibriation that exhausts its budget
   aborts the whole call rather than sampling anyway. A budget below the 400-sweep comparison
   window makes that deterministic. The operator is a call counter: it stays at zero, which is
   what shows neither the correlation-time stage nor the sampling loop ran at all -- $Failed
   alone would not, since it is also what a run that sampled and then gave up would return. *)
VerificationTest[
    Module[{calls = 0, out, noneql, noconv},
        noneql = messageArguments["RandomPureSimplicialComplexMCMC::noneql",
            noconv = messageArguments["RandomPureSimplicialComplexMCMCEquilibriate::noconv",
                Quiet@Block[{ECGrav`$ECGravMaxEquilibriationSweeps = 100}, SeedRandom[602];
                    out = Block[{Print = Null &}, ECGrav`RandomPureSimplicialComplexMCMC[
                        mcmcSeed24, {Function[c, ++calls]}, 1, 5]]]]];
        {out, calls, noneql, Length[noconv]}],
    {$Failed, 0, {{100}}, 1},
    TestID -> "RandomPureSimplicialComplexMCMC-aborts-when-not-equilibriated"
];

(* The fixed-vertex overload holds the support at the seed's vertex set for every sample it
   returns, and its third column is the edge count rather than the vertex count. *)
VerificationTest[
    Module[{ops, meas, data}, SeedRandom[605];
        ops = {Function[c, Total[Flatten[c]]]};
        {meas, data} = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[mcmcSeed34, True, ops, 1, 5]];
        {Length /@ meas,
         AllTrue[meas[[All, -1]], Union @@ # === Range[7] &],
         mcmcMeasurementsConsistentQ[meas, ops, True, 1],
         meas[[All, 3]] === (Length[Union @@ (Subsets[#, {2}] & /@ #)] & /@ meas[[All, -1]]),
         Length[DeleteDuplicates[meas[[All, -1]]]] > 1,
         data["converged"], data["eqlT"], Length[data["corrTValues"]]}],
    {{5, 5, 5, 5, 5}, True, True, True, True, True, 31, 3},
    TestID -> "RandomPureSimplicialComplexMCMC-fixed-vertex-support"
];

(* Handed a previous run's association the driver skips both preparatory stages and samples
   from where that run left off. The check that it really skipped them is the diagnostic plot:
   equilibriation is the only stage that Prints one, so the full pipeline emits at least one
   and the continuation must emit none. It carries the earlier run's eqlT and corrT forward
   and advances the state rather than restarting from it. *)
VerificationTest[
    Module[{ops, pipelinePlots, contPlots, meas, data, meas2, data2},
        ops = {Function[c, Total[Flatten[c]]]};
        SeedRandom[603];
        {pipelinePlots, {meas, data}} = emittedGraphicsAndResult[
            Quiet@ECGrav`RandomPureSimplicialComplexMCMC[mcmcSeed24, ops, 1, 4]];
        SeedRandom[604];
        {contPlots, {meas2, data2}} = emittedGraphicsAndResult[
            Quiet@ECGrav`RandomPureSimplicialComplexMCMC[data, ops, 1, 4]];
        {pipelinePlots >= 1, contPlots, Length[meas2], meas2[[All, 1]],
         mcmcMeasurementsConsistentQ[meas2, ops, 1],
         {data2["eqlT"], data2["corrT"]} === {data["eqlT"], data["corrT"]},
         data2[["state", "complex"]] =!= data[["state", "complex"]]}],
    {True, 0, 4, Range[4], True, True, True},
    TestID -> "RandomPureSimplicialComplexMCMC-continues-from-a-prior-run"
];

(* corrT is not just reported, it is the spacing between recorded samples: the driver runs
   that many sweeps between them, which is the whole point of measuring it. Nothing else here
   pins that -- a driver that reported corrT and then swept a constant number of times would
   still return self-consistent rows.

   The continuation overload takes the association as given, so feeding it corrT -> 0 makes
   the sweep run no steps at all and freezes the chain: every sample must then come back
   identical, and identical to the state it started from. corrT -> 6 is the control that the
   same call does move when it is allowed to. Zero is a probe, not a value the two stages can
   produce -- the correlation time is floored at 2. *)
VerificationTest[
    Module[{ops = {Function[c, Total[Flatten[c]]]}, prior, frozen, moving},
        SeedRandom[608];
        prior = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[mcmcSeed24, ops, 1, 2]][[2]];
        frozen = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[Append[prior, "corrT" -> 0], ops, 1, 3]][[1]];
        moving = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[Append[prior, "corrT" -> 6], ops, 1, 3]][[1]];
        {DeleteDuplicates[frozen[[All, -1]]] === {prior[["state", "complex"]]},
         Length[DeleteDuplicates[Rest /@ frozen]],   (* rows differ only in the sweep index *)
         Length[DeleteDuplicates[moving[[All, -1]]]] > 1}],
    {True, 1, True},
    TestID -> "RandomPureSimplicialComplexMCMC-corrT-sets-the-sample-spacing"
];

(* HoldNumberOfVerticesFixed -> False delegates to the shorter overload, for both the
   seed-complex and the prior-run forms. *)
VerificationTest[
    Module[{ops = {Function[c, Total[Flatten[c]]]}, a, b, prior, c1, c2},
        SeedRandom[606];
        a = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[mcmcSeed24, False, ops, 1, 4]];
        SeedRandom[606];
        b = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[mcmcSeed24, ops, 1, 4]];
        prior = a[[2]];
        SeedRandom[607];
        c1 = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[prior, False, ops, 1, 4]];
        SeedRandom[607];
        c2 = Quiet@Block[{Print = Null &},
            ECGrav`RandomPureSimplicialComplexMCMC[prior, ops, 1, 4]];
        {a === b, c1 === c2}],
    {True, True},
    TestID -> "RandomPureSimplicialComplexMCMC-free-vertex-delegates"
];
