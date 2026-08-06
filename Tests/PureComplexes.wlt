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
VerificationTest[
    Block[{ECGrav`Private`$RandomFacetLabeledParallelThreshold = 10},
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
   with 12 labellings, and the star, with 4 -- so a uniform labeled draw would split them 75/25
   and only the rejection weight |Aut|/n! brings them to 50/50. Drawn 20 at a time to stay on
   the serial branch, which keeps the seed reproducible. 400 samples is deliberately modest: the
   failure this guards against is the rejection weight being dropped or inverted, which lands at
   75/25 and gives chi-square about 100 on one degree of freedom, so it is caught with room to
   spare while keeping the draw near two seconds. *)
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
