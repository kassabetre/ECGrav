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

(* NumPureComplexes[p,q,n] counts q-element sets of p-subsets of [n] whose union is all of
   [n]. The values below were checked by direct enumeration of exactly that
   (Select[Subsets[Subsets[Range[n],{p}],{q}], Union@@#===Range[n]&]), so they are oracles
   independent of the recursion. NumPureComplexes[4,3,12]=5775 is independently 12!/(4!^3 3!),
   the number of ways to partition 12 labelled vertices into 3 blocks of 4. *)

VerificationTest[ECGrav`NumPureComplexes[3, 3], 2649, TestID -> "NumPureComplexes-3-3"];
VerificationTest[ECGrav`NumPureComplexes[3, 3, 7], 945, TestID -> "NumPureComplexes-3-3-7"];
VerificationTest[ECGrav`NumPureComplexes[3, 3, 9], 280, TestID -> "NumPureComplexes-3-3-9"];
VerificationTest[ECGrav`NumPureComplexes[2, 4, 6], 330, TestID -> "NumPureComplexes-2-4-6"];
VerificationTest[ECGrav`NumPureComplexes[3, 4, 8], 72380, TestID -> "NumPureComplexes-3-4-8"];
VerificationTest[ECGrav`NumPureComplexes[2, 4], 900, TestID -> "NumPureComplexes-2-4"];
VerificationTest[ECGrav`NumPureComplexes[3, 4], 441061, TestID -> "NumPureComplexes-3-4"];

(* n outside p <= n <= p q is answered without building a row; both ends must give 0. *)
VerificationTest[ECGrav`NumPureComplexes[3, 2, 2], 0, TestID -> "NumPureComplexes-n-below-p"];
VerificationTest[ECGrav`NumPureComplexes[3, 2, 7], 0, TestID -> "NumPureComplexes-n-above-pq"];

(* q<=0: the empty complex covers no vertices, so N(p,0,n) is 1 at n==0 and 0 otherwise, and
   summing it over n gives 1. Negative q counts nothing. These returned Indeterminate (3-arg)
   and an unevaluated Total[Table[..]] (2-arg) before. *)
VerificationTest[ECGrav`NumPureComplexes[3, 0, 0], 1, TestID -> "NumPureComplexes-q0-n0"];
VerificationTest[ECGrav`NumPureComplexes[3, 0, 5], 0, TestID -> "NumPureComplexes-q0-n-nonzero"];
VerificationTest[ECGrav`NumPureComplexes[3, -1, 5], 0, TestID -> "NumPureComplexes-q-negative"];
VerificationTest[ECGrav`NumPureComplexes[3, 0], 1, TestID -> "NumPureComplexes-2arg-q0"];
VerificationTest[ECGrav`NumPureComplexes[2, -1], 0, TestID -> "NumPureComplexes-2arg-q-negative"];

(* Degenerate purity: no n in p..p q carries enough distinct p-subsets, so the sum is empty. *)
VerificationTest[ECGrav`NumPureComplexes[0, 3], 0, TestID -> "NumPureComplexes-2arg-p0"];

(* The 2-arg overload must agree with summing the 3-arg form over every vertex count. *)
VerificationTest[
    ECGrav`NumPureComplexes[3, 5],
    Total[Table[ECGrav`NumPureComplexes[3, 5, n], {n, 0, 15}]],
    TestID -> "NumPureComplexes-2arg-equals-row-total"
];

(* Rows are cached per (p,q) behind a high-water mark, so a descending sweep of q must return
   the same values as an ascending one. *)
VerificationTest[
    Table[ECGrav`NumPureComplexes[4, q, 12], {q, 9, 1, -1}],
    Reverse[{0, 0, 5775, 57345750, 29179710945, 5551975972620, 627962951769615,
        51166525553400015, 3312982055194610645}],
    TestID -> "NumPureComplexes-descending-q"
];

(* Rows are filled by an explicit bottom-up loop, so recursion depth no longer grows with q;
   the earlier recursive form died on TerminatedEvaluation[RecursionLimit] before q=700.
   p=1 keeps the entries tiny: q singleton facets cover exactly q vertices. *)
VerificationTest[ECGrav`NumPureComplexes[1, 800, 800], 1, TestID -> "NumPureComplexes-deep-q"];

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
            Table[ECGrav`NumPureComplexes[3, q, 14], {q, 1, 12}]],
        expected,
        TestID -> "NumPureComplexes-cap-below-range"
    ];
    VerificationTest[
        Block[{ECGrav`Private`$NumPCMaxCachedQ = 7},
            Table[ECGrav`NumPureComplexes[3, q, 14], {q, 12, 1, -1}]],
        Reverse[expected],
        TestID -> "NumPureComplexes-cap-straddling-descending"
    ];
    VerificationTest[
        Block[{ECGrav`Private`$NumPCMaxCachedQ = 10000},
            Table[ECGrav`NumPureComplexes[3, q, 14], {q, 1, 12}]],
        expected,
        TestID -> "NumPureComplexes-cap-above-range"
    ];
];

(* A cap of 0 switches caching off altogether and must still give the same answer. *)
VerificationTest[
    Block[{ECGrav`Private`$NumPCMaxCachedQ = 0}, ECGrav`NumPureComplexes[3, 8, 14]],
    465021334918140,
    TestID -> "NumPureComplexes-caching-disabled"
];

(* NumPCClearCache[] reclaims the rows and returns the bytes freed; nothing else holds a row,
   so the next call must rebuild to the same values. *)
VerificationTest[
    (ECGrav`NumPureComplexes[4, 60, 120]; ECGrav`Private`NumPCClearCache[] > 0),
    True,
    TestID -> "NumPureComplexes-clear-cache-reclaims"
];
VerificationTest[
    (ECGrav`Private`NumPCClearCache[]; ECGrav`Private`NumPCClearCache[]),
    0,
    TestID -> "NumPureComplexes-clear-cache-idempotent"
];
VerificationTest[
    (ECGrav`Private`NumPCClearCache[];
     {ECGrav`NumPureComplexes[3, 4, 8], ECGrav`NumPureComplexes[3, 3]}),
    {72380, 2649},
    TestID -> "NumPureComplexes-correct-after-clear"
];

(* :!CodeAnalysis::EndBlock:: *)
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
   NB: RandomUniformUnlabeledPureSimplicialComplex and RandomUniformFacetLabeledPureSimplicialComplex
   are verified interactively but deliberately kept OUT of the suite: their first call in a
   fresh kernel pays a large one-time parallel-distribution cost (~30 s), not worth a
   structural smoke test here. RandomVertexLabeledPureSimplicialComplex IS covered below --
   its serial path pays no such cost (7 ms in a fresh kernel), and its parallel branch costs
   about 0.1 s once TestPrelude has launched the kernels. *)

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
   branch used to distribute only ECGrav`Private`, so the subkernels had no NumPureComplexes,
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
   NumPureComplexes-weighted composition draw buys. Enumerate the 16 complexes at p=2,q=3,n=4
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

(* An empty sample space must fail loudly. Every composition weight carries a
   NumPureComplexes factor, so on an empty space they are all zero and the draw cannot proceed;
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
