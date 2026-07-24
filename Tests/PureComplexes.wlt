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

VerificationTest[ECGrav`NumPureComplexes[3, 3], 2649, TestID -> "NumPureComplexes-3-3"];
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
