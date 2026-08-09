(* ::Package:: *)

(* ::Input:: *)
(*(*:Name: PureComplexes - Functions for doing combinatorics on simplical complexes*)*)
(*(*:Author: Kassahun Betre*)*)
(*(*:Contributors: Khai Luong*)*)
(*(*:Date: 08/11/2025*)*)
(*(*: To install package: *)
(*1. Save as ECGrav.wl, *)
(* 2. Evaluate Notebook, *)
(* 3. File->Install->Choose: Type = Package, Source=From File, Select ECGrav.wl,     Instal Name = ECGrav  *)*)
(*(*: To load package evaluate Needs["ECGrav`"] in a new notebook*)*)
(**)


(* ::Title:: *)
(*Begin PureComplexes Package*)


(* ::Title:: *)
(*PureComplexes Public Functions*)


(* ::Subtitle:: *)
(*Main PureComplexes Public Symbols*)


(* ::Title:: *)
(*PureComplexes Private*)


(* Private helper messages ------------------------------------------------------
   Error messages for the internal connected-complex automorphism/stabilizer
   helpers (SimplicialComplexAutomorphismGroupOrderConn, PureComplexAutomorphismGroupOrderConn,
   PureComplexFacetAutomorphismGroupOrderConn, CliqueFacetStabilizerGroupOrderConn,
   CliqueFacetAutomorphismGroupOrderConn). Being undeclared (no usage message), these
   resolve into the package-private ECGrav`Private` context and are NOT public: the
   public non-Conn functions call them on each already-connected component. The messages
   are defensive diagnostics only. *)

SimplicialComplexAutomorphismGroupOrderConn::argerr="Internal helper: input should be a connected complex of the form SimplicialComplexAutomorphismGroupOrderConn[facetsLst_List].";
PureComplexAutomorphismGroupOrderConn::argerr="Internal helper: input should be a connected pure complex of the form PureComplexAutomorphismGroupOrderConn[facetsLst_List].";
PureComplexFacetAutomorphismGroupOrderConn::argerr="Internal helper: input should be a connected pure complex of the form PureComplexFacetAutomorphismGroupOrderConn[facetsLst_List].";
CliqueFacetStabilizerGroupOrderConn::argerr="Internal helper: input should be a connected clique complex of the form CliqueFacetStabilizerGroupOrderConn[facetsLst_List].";
CliqueFacetAutomorphismGroupOrderConn::argerr="Internal helper: input should be a connected clique complex of the form CliqueFacetAutomorphismGroupOrderConn[facetsLst_List].";


(* ::Chapter:: *)
(*Helper Functions*)


(* ::Section::Closed:: *)
(*General*)


(* ::Subsection::Closed:: *)
(*Basic Simplicial Complex Constructions*)


(* ::Item::Closed:: *)
(*FacetIncidenceMatrix*)


(* Primary Pattern *)
FacetIncidenceMatrix[facetsLst_List]:=
With[{facetOrder=Length[facetsLst],vlist=DeleteDuplicates[Flatten[facetsLst]]},
	Table[Table[If[MemberQ[i,j],1,0],{j,vlist}],{i,facetsLst}]
];

(* Overload Pattern *)
FacetIncidenceMatrix[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},
	FacetIncidenceMatrix[clqs]
];

(* Catch-all Pattern *)
FacetIncidenceMatrix[args___]:=(Message[FacetIncidenceMatrix::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetAdjacencyMatrix*)


(* Primary Pattern *)
FacetAdjacencyMatrix[facetsLst_List]:=Module[{clqOrder=Length[facetsLst],facetAdjMat=DiagonalMatrix[Length/@facetsLst]},
Do[
Do[
facetAdjMat[[i,j]]=facetAdjMat[[j,i]]=Length[Intersection[facetsLst[[i]],facetsLst[[j]]]]
,{j,i+1,clqOrder}]
,{i,1,clqOrder-1}];
facetAdjMat
];

(* Overload Pattern *)
FacetAdjacencyMatrix[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},
FacetAdjacencyMatrix[clqs]
];

(* Catch-all Pattern *)
FacetAdjacencyMatrix[args___]:=(Message[FacetAdjacencyMatrix::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*GraphFromCliques*)


(* Primary Pattern *)
GraphFromCliques[clqs_List]:=With[{vertexlist=DeleteDuplicates[Flatten[clqs]],edgelist=UndirectedEdge@@@(DeleteDuplicates[Flatten[Subsets[Sort[#],{2}]&/@clqs,1]])},
Graph[vertexlist,edgelist]
];

(* Catch-all Pattern *)
GraphFromCliques[args___]:=(Message[GraphFromCliques::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*GraphFromFacetIncidence*)


(* Primary Pattern *)
GraphFromFacetIncidence[incidenceMat_List/;MatrixQ[incidenceMat]]:=
With[{clqsLst=Table[Pick[Range[Length[incidenceMat[[1]]]],i,1],{i,incidenceMat}]},
GraphFromCliques[clqsLst]
];

(* Catch-all Pattern *)
GraphFromFacetIncidence[args___]:=(Message[GraphFromFacetIncidence::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliquesFromFacetIncidence*)


(* Primary Pattern *)
CliquesFromFacetIncidence[incidenceMat_List/;MatrixQ[incidenceMat]]:=
With[{clqsLst=Table[Pick[Range[Length[incidenceMat[[1]]]],i,1],{i,incidenceMat}]},
clqsLst
];

(* Catch-all Pattern *)
CliquesFromFacetIncidence[args___]:=(Message[CliquesFromFacetIncidence::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ComplexFromFacetLabeledVertexList*)


(* Primary Pattern *)
ComplexFromFacetLabeledVertexList[facetLabeledVertices_List]:=
With[{facetLabels=DeleteDuplicates[Flatten[facetLabeledVertices]],cmlxAsn=<|Table[i->facetLabeledVertices[[i]],{i,1,Length[facetLabeledVertices]}]|>},
Sort[
	Sort/@Table[
			With[{k=Select[facetLabeledVertices,MemberQ[#,q]&]},
				Keys[Select[cmlxAsn,MemberQ[k,#]&]]
			],{q,facetLabels}
		]
	]
];

(* Catch-all Pattern *)
ComplexFromFacetLabeledVertexList[args___]:=(Message[ComplexFromFacetLabeledVertexList::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetLabeledVertexListFromComplex*)


(* Primary Pattern *)
FacetLabeledVertexListFromComplex[facetsLst_List]:=
With[{vlist=DeleteDuplicates[Flatten[facetsLst]],
cmlxAsn=<|Table[i->facetsLst[[i]],{i,1,Length[facetsLst]}]|>},
(*ReverseSort[Table[Keys[Select[cmlxAsn,MemberQ[#,v]&]],{v,vlist}]]*)
Join@@(Map[LexicographicSort[#]&,
GatherBy[ReverseSortBy[Table[Keys[Select[cmlxAsn,MemberQ[#,v]&]],{v,vlist}],Length],Length],{1}])
];

(* Catch-all Pattern *)
FacetLabeledVertexListFromComplex[args___]:=(Message[FacetLabeledVertexListFromComplex::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexQ*)


(* Primary Pattern *)
PureComplexQ[facets_List]:=If[Length[DeleteDuplicates[Length/@facets]]==1,True,False];

(* Catch-all Pattern *)
PureComplexQ[args___]:=(Message[PureComplexQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureGraphQ*)


(* Primary Pattern *)
PureGraphQ[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},
If[Length[DeleteDuplicates[Map[Length,clqs]]]==1,True,False]
];

(* Catch-all Pattern *)
PureGraphQ[args___]:=(Message[PureGraphQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueComplexQ*)


(* Primary Pattern *)
CliqueComplexQ[facetsLst_List]:=
(*Given a list of maximal cliques, it checks whether or not it passes the clique condition test*)
Module[{threeSubsets=Subsets[facetsLst,{3}],fIsInAClique},
Catch[If[Length[facetsLst]<3,Throw[True, "ECGravReturn$1"]];
fIsInAClique[candClq_List]:=AnyTrue[facetsLst,SubsetQ[#,candClq]&];
AllTrue[threeSubsets,fIsInAClique[Union[Intersection[#[[1]],#[[2]]],Intersection[#[[1]],#[[3]]],Intersection[#[[2]],#[[3]]]]]&], "ECGravReturn$1"]
];

(* Catch-all Pattern *)
CliqueComplexQ[args___]:=(Message[CliqueComplexQ::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Spheres and Balls, Links and Stars*)


(* ::Item::Closed:: *)
(*Sph*)


(* Primary Pattern *)
Sph[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer]:=
With[{size = Length[Amat],rowcolsToKeep=Flatten[Position[Amat[[i]],1]]},
	Catch[If[i>size,Message[Sph::vtxnotfound, i];
	Throw[{}, "ECGravReturn$2"]];
	If[rowcolsToKeep=={},Throw[{}, "ECGravReturn$2"]];
	If[size==0,Throw[{}, "ECGravReturn$2"]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep], "ECGravReturn$2"]
];

(* Overload Pattern *)
Sph[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer,r_Integer]:=
Module[{rowcolsToKeep,dm},
	Catch[If[Length[Amat]==0,Throw[{}, "ECGravReturn$3"]];
	dm=GraphDistanceMatrix[AdjacencyGraph[Amat]];
	rowcolsToKeep=Flatten[Position[dm[[i]],r]];
	If[rowcolsToKeep=={},Throw[{}, "ECGravReturn$3"]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep], "ECGravReturn$3"]
];

(* Overload Pattern *)
Sph[g_Graph,i_Integer]:=
Module[{},
	Catch[If[VertexCount[g]==0,Throw[{}, "ECGravReturn$4"]];
	If[MemberQ[VertexList[g],i]==False,Message[Sph::vtxnotfound, i];{}
	];

	Subgraph[g,AdjacencyList[g,i]], "ECGravReturn$4"]
];

(* Overload Pattern *)
Sph[g_Graph,i_Integer,r_Integer]:=
With[{sphVertices=Select[VertexList[g],GraphDistance[g,i,#]==r&]},
	
	Catch[If[VertexCount[g]==0,Throw[{}, "ECGravReturn$5"]];
	If[MemberQ[VertexList[g],i]==False,Message[Sph::vtxnotfound, i];{}
	];

	Subgraph[g,sphVertices], "ECGravReturn$5"]
];

(* Catch-all Pattern *)
Sph[args___]:=(Message[Sph::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Bll*)


(* Primary Pattern *)

Bll[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer]:=
With[{size = Length[Amat],rowcolsToKeep=Flatten[Position[Amat[[i]],1]]},
	Catch[If[i>size,Message[Bll::vtxnotfound, i];
	Throw[{}, "ECGravReturn$6"]];
	If[size==0,Throw[{}, "ECGravReturn$6"]];
	If[rowcolsToKeep=={},Throw[{}, "ECGravReturn$6"]];
	Join[
		Transpose[
			Join[
				Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep],
				Table[1,{Length[rowcolsToKeep]}]
			]
		],
		PadRight[Table[1,{Length[rowcolsToKeep]}],Length[rowcolsToKeep]+1]
	], "ECGravReturn$6"]
];

(* Overload Pattern *)
Bll[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer,r_Integer]:=
Module[{rowcolsToKeep,dm},
	Catch[If[Length[Amat]==0,Throw[{}, "ECGravReturn$7"]];
	dm=GraphDistanceMatrix[AdjacencyGraph[Amat]];
	rowcolsToKeep=Flatten[Position[dm[[1]],_?(#<=r&)]];
	If[rowcolsToKeep=={},Throw[{}, "ECGravReturn$7"]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep], "ECGravReturn$7"]
];

(* Overload Pattern *)
Bll[g_Graph,i_Integer]:=
(*Unit ball in a graph g at vertex i*)
With[{unitBallVertices=Union[{i},AdjacencyList[g,i]]},
	Catch[If[VertexCount[g]==0,Throw[{}, "ECGravReturn$8"]];
	If[MemberQ[VertexList[g],i]==False,Message[Bll::vtxnotfound, i];{}
	];

Subgraph[g,unitBallVertices], "ECGravReturn$8"]

];

(* Overload Pattern *)
Bll[g_Graph,i_Integer,r_Integer]:=
(*Ball of radius r in a graph g at vertex i*)
With[{ballVertices=Select[VertexList[g],GraphDistance[g,i,#]<=r&]},
	Catch[If[VertexCount[g]==0,Throw[{}, "ECGravReturn$9"]];
	If[MemberQ[VertexList[g],i]==False,Message[Bll::vtxnotfound, i];{}
	];

	Subgraph[g,Union[{i},ballVertices]], "ECGravReturn$9"]

];

(* Catch-all Pattern *)
Bll[args___]:=(Message[Bll::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Lnk*)


Lnk[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1}),q_List/;Depth[q]==2]:=
(*Returns the link of the face q in the complex given by the list of facets facetsLst*)
With[{link=Select[facetsLst,SubsetQ[#,q]&],rules=Table[i->Nothing,{i,q}]},
	link/.rules
];

(* Catch-all Pattern *)
Lnk[args___]:=(Message[Lnk::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Str*)


Str[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1}),q_List/;Depth[q]==2]:=
(*Returns the link of the face q in the complex given by the list of facets facetsLst*)
Select[facetsLst,SubsetQ[#,q]&];

(* Catch-all Pattern *)
Str[args___]:=(Message[Str::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Basic Graph Observables*)


(* ::Item::Closed:: *)
(*Deg*)


(* Primary Pattern *)
Deg[Amat_List,i_Integer]:=(*Given an adjacency matrix of a graph Amat, and a vertex i, it gives the degree of the vertex*)
Total[Amat[[i]]];

(* Catch-all Pattern *)
Deg[args___]:=(Message[Deg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*AvgDeg*)


(* Primary Pattern *)
AvgDeg[Amat_List/;!GraphQ[Amat]]:=
(*The average degree of a graph given as an adjacency matrix Amat*)
(1/Length[Amat])*Total[Amat,2];

(* Overload Pattern *)
AvgDeg[g_Graph]:=
(*The average degree of a graph given as a graph object g*)
(1/VertexCount[g])*Total[VertexDegree[g]];

(* Catch-all Pattern *)
AvgDeg[args___]:=(Message[AvgDeg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetDeg*)


(* Primary Pattern *)
FacetDeg[facetsLst_List,v_Integer]:=
(*Given a simplicial complex as a list facetsLst of facets and a vertex v, 
it computes the facet degree (number of facets containing the vertex) *)
Length[Select[facetsLst,MemberQ[#,v]&]];

(* Overload Pattern *)
FacetDeg[facetsLst_List]:=
(*Given a simplicial complex, it lists the facet degrees of all vertices  *)
With[{vertices=Sort[DeleteDuplicates[Flatten[facetsLst]]]},
FacetDeg[facetsLst,#]&/@vertices
];

(* Catch-all Pattern *)
FacetDeg[args___]:=(Message[FacetDeg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HyperDeg*)


(* Primary Pattern *)
HyperDeg[g_Graph,clq_List]:=
(*Computes the degree of the clique clq given as a list of vertices. 
Checks whether or not the input is a clique in the graph*)
With[{spheres=Table[AdjacencyList[g,v],{v,clq}]},
Catch[If[CompleteGraphQ[Subgraph[g,clq]]==False,Throw[Null, "ECGravReturn$10"]];
Length[Intersection@@spheres], "ECGravReturn$10"]
];

(* Overload Pattern *)
HyperDeg[Amat_List/;(SymmetricMatrixQ[Amat]&&Sort[DeleteDuplicates[Flatten[Amat]]]=={0,1}),clq_List]:=
(*Computes the degree of the clique clq given as a list of vertices.
Checks whether or not the input is a clique in the graph.

Pure adjacency-matrix implementation: no Graph object is built, so this no longer
pays for AdjacencyGraph + Subgraph + CompleteGraphQ (twice - the old version ran the
CompleteGraphQ test here and then again inside the Graph overload it delegated to).

clqV is a clique iff its induced submatrix Amat[[clqV,clqV]] has all Length[clqV]
(Length[clqV]-1) off-diagonal entries equal to 1, i.e. iff those entries sum to that
value (the diagonal of an adjacency matrix is 0, so it contributes nothing).

The common neighbours of clqV are the columns where every row Amat[[v]], v in clqV,
carries a 1, so their count is Total[Times@@Amat[[clqV]]]. Members of clqV are excluded
automatically: for v in clqV the entry Amat[[v,v]] is 0, which zeroes that column of the
product. For a single vertex this correctly reduces to its degree.*)
With[{clqV=DeleteDuplicates[clq]},
Catch[If[clqV=={},Throw[Null, "ECGravReturn$11"]];
If[Total[Amat[[clqV,clqV]],2]!=Length[clqV](Length[clqV]-1),Throw[Null, "ECGravReturn$11"]];
Total[Times@@Amat[[clqV]]], "ECGravReturn$11"]
];

(* Overload Pattern *)
HyperDeg[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&
	Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1}),clq_List]:=
(*Computes the hyperdegree of the face clq given as a list of vertices. 
Checks whether or not the input is a face in the graph*)
With[{lnk=Lnk[facetsLst,clq]},
	Catch[If[NoneTrue[facetsLst,SubsetQ[#,clq]&],Throw[Null, "ECGravReturn$12"]];
	Total[Length/@lnk], "ECGravReturn$12"]
];

(* Catch-all Pattern *)
HyperDeg[args___]:=(Message[HyperDeg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FVector*)


(* Primary Pattern *)
FVector[g_Graph]:=
(*Computes the f-vector of the graph g which is the vector {f0,f1,...} of 
the number of vertices, number of edges, etc. *)
With[{n0=VertexCount[g],n1=EdgeCount[g],maxClqs=FindClique[g,\[Infinity],All]},
	Catch[If[n1==0,Throw[{n0}, "ECGravReturn$13"]];
	Join[{n0,n1},
		Table[Length[Union@@(Subsets[#,{q}]&/@maxClqs)],{q,3,Max[Length/@maxClqs]}]
	], "ECGravReturn$13"]
];

(* Overload Pattern *)
FVector[Amat_List/;(SymmetricMatrixQ[Amat]&&Sort[DeleteDuplicates[Flatten[Amat]]]=={0,1})]:=
(*Computes the f-vector of the graph given by its adjacency matrix Amat which is the 
vector {f0,f1,...} of the number of vertices, number of edges, etc. *)
With[{n0=Length[Amat],n1=Total[Amat,2]/2,maxClqs=FindClique[AdjacencyGraph[Amat],\[Infinity],All]},
Catch[If[n1==0,Throw[{n0}, "ECGravReturn$14"]];
Join[{n0,n1},Table[Length[Union@@(Subsets[#,{q}]&/@maxClqs)],{q,3,Max[Length/@maxClqs]}]], "ECGravReturn$14"]
];

(* Overload Pattern *)
FVector[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&
	Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1})]:=
(*Computes the hyperdegree of the face clq given as a list of vertices. 
Checks whether or not the input is a face in the graph*)
With[{n0=Length[DeleteDuplicates[Flatten[facetsLst]]],maxFacetSize=Max[Length/@facetsLst]},
Catch[If[n0==0,Throw[{}, "ECGravReturn$15"]];
Join[{n0},Table[Length[Union@@(Subsets[#,{q}]&/@facetsLst)],{q,2,maxFacetSize}]], "ECGravReturn$15"]
];

(* Catch-all Pattern *)
FVector[args___]:=(Message[FVector::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*GVolume*)


(* Primary Pattern *)
GVolume[Amat_List,Dmat_List,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, where the distance 
between nodes is given by Dmat, centered at node i and having radius r*)
With[{n=Length[Amat],inball=Table[If[Dmat[[i,j]]<=r,1,0],{j,n}]},
Total[inball]-1
];

(* Overload Pattern *)
GVolume[g_Graph,Dmat_List,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, where the distance 
between nodes is given by Dmat, centered at node i and having radius r*)
With[{n=VertexCount[g],inball=Table[If[Dmat[[i,j]]<=r,1,0],{j,n}]},
Total[inball]-1
];

(* Overload Pattern *)
GVolume[Amat_List,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, centered at node i 
and having radius r*)
With[{n=Length[Amat],dmat=GraphDistanceMatrix[AdjacencyGraph[Amat]]},
Total[Table[If[dmat[[i,j]]<=r,1,0],{j,n}]]-1
];

(* Overload Pattern *)
GVolume[g_Graph,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, centered at node i 
and having radius r*)
With[{n=VertexCount[g],dmat=GraphDistanceMatrix[g]},
Total[Table[If[dmat[[i,j]]<=r,1,0],{j,n}]]-1
];

(* Catch-all Pattern *)
GVolume[args___]:=(Message[GVolume::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*KFaceDistance*)


(* Primary Pattern *)
KFaceDistance[facets_List/;Depth[facets]==3,clq1_List/;Depth[clq1]==2,clq2_List/;Depth[clq2]==2]:=

(*Given a list of maximal cliques (facets) of a clique complex or a graph and two equal
 cardinality faces, clq1 and clq2 of cardinality d in the graph,  
 this method computes the minimum clique distance between them where 
 the path steps on only (d+1)-dimensional faces joined along d-dimensional subfaces.  
 It  uses the built in Mathematica GraphDistance function on a graph constructed 
 such that the vertices are all d-faces and edges between them are 1 if two faces 
 are contained in a bigger clique and 0 otherwise. *)

Module[{d=Length[clq1],dests,pos1,pos2,clqGraphAmat},
Catch[If[Complement[clq1,clq2]=={},Throw[0, "ECGravReturn$16"]];
(*If the two cliques are the same sets, then the distance between them is zero.*)

dests=Union@@(Subsets[#,{d}]&/@facets);

pos1=Position[dests,Sort[clq1]][[1,1]];
pos2=Position[dests,Sort[clq2]][[1,1]];

clqGraphAmat = 
	Table[
		Table[
			If[Complement[i,j]=={},0,
				If[(Length[Union[i,j]]==d+1)&&AnyTrue[facets,SubsetQ[#,Union[i,j]]&],1,0]
				],{i,dests}
			],
{j,dests}];

GraphDistance[AdjacencyGraph[clqGraphAmat],pos1,pos2], "ECGravReturn$16"]

];

(* Catch-all Pattern *)
KFaceDistance[args___]:=(Message[KFaceDistance::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*KFaceDistanceMatrix*)


(* Primary Pattern *)
KFaceDistanceMatrix[facets_List, facedim_Integer]:=
(*Given a list of maximal cliques (facets) of a clique complex or a graph and 
dimensionality clqdim, this method computes the clique distance matrix between 
every pair of faces of cardinality clqdim, where each path steps on only 
(facedim+1)-dimensional faces.  It  uses the built in Mathematica 
GraphDistanceMatrix function on a graph constructed such that the vertices 
are all clqdim-faces and edges between them are 1 if two faces are contained in a 
bigger clique and 0 otherwise. *)

Module[{dests,clqGraphAmat},

Catch[If[facedim==1,Throw[GraphDistanceMatrix[GraphFromCliques[facets]], "ECGravReturn$17"]];

dests=Union@@(Subsets[#,{facedim}]&/@facets);

clqGraphAmat = 
	Table[
		Table[
			If[Complement[i,j]=={},0,
				If[(Length[Union[i,j]]==facedim+1)&&AnyTrue[facets,SubsetQ[#,Union[i,j]]&],1,0]
			],{i,dests}
		],
{j,dests}];

GraphDistanceMatrix[AdjacencyGraph[clqGraphAmat]], "ECGravReturn$17"]

];

(* Catch-all Pattern *)
KFaceDistanceMatrix[args___]:=(Message[KFaceDistanceMatrix::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*KpathConnectedComponents*)


(* Primary Pattern *)
KpathConnectedComponents[facets_List,k_Integer]:=

(*(*****************************)
(* Last Updated: 01/17/2026  *)
(*****************************)*)
(*Given a complex as a list of facets or a graph and dimensionality k, this method 
computes the vertices contained in (k+1)-path connected components where a 
(k+1)-path connected component of the complex contains all k-faces that form a 
component connected by a path of (k+1)-dimensional faces.  It outputs a list of list 
of vertices in each (k+1)-path connected component.
It  uses the built in Mathematica ConnectedComponents function on a graph 
constructed such that the vertices are all k-faces and two such vertices are
connected if they are contained in a bigger face. *)

Module[{tooSmallFacets,isolated,destinations,clqGraphAmat},

Catch[If[k==1,Throw[ConnectedComponents[GraphFromCliques[facets]], "ECGravReturn$18"]];

tooSmallFacets=Select[facets,Length[#]<=k&];
isolated={#}&/@
	Complement[DeleteDuplicates[Flatten[tooSmallFacets]],
		DeleteDuplicates[Flatten[Complement[facets,tooSmallFacets]]
		]
	];
destinations=Union@@(Subsets[#,{k}]&/@facets);

clqGraphAmat = 
	Table[
		Table[
			If[Complement[i,j]=={},0,
				If[(Length[Union[i,j]]==k+1)&&AnyTrue[facets,SubsetQ[#,Union[i,j]]&],1,0]
			],{i,destinations}
		],
{j,destinations}];

Join[
	(Union@@Part[destinations,#])&/@ConnectedComponents[AdjacencyGraph[clqGraphAmat]],
	isolated
	], "ECGravReturn$18"]

];

(* Overload Pattern *)
KpathConnectedComponents[g_Graph,k_Integer]:=
With[{clqs=FindClique[g,\[Infinity],All]},
	KpathConnectedComponents[clqs,k]
];

(* Catch-all Pattern *)
KpathConnectedComponents[args___]:=(Message[KpathConnectedComponents::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ConnectedComplexComponents*)


(* Primary Pattern *)
ConnectedComplexComponents[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it gives a list of simplicial
 complexes which are connected. 
 E.g.ECGrav`ConnectedComplexComponents[{{1,2},{3,4}}] = {{{1,2}},{{3,4}}} *)
With[{connectedvertices=ConnectedComponents[GraphFromCliques[facelst]]},
Table[Select[facelst,SubsetQ[i,#]&],{i,connectedvertices}]
];

(* Catch-all Pattern *)
ConnectedComplexComponents[args___]:=(Message[ConnectedComplexComponents::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FractionInLargestComponent*)


(* Primary Pattern *)
FractionInLargestComponent[g_Graph]:=
(*
(*Notes: *)
(*Given a graph, it computes the ratio of the number of vertices in the largest 
connected component to th etotal numbe rof vertices in the graph.  *)*)
With[{numV=VertexCount[g],largestComponentLength=Max[Length/@(ConnectedComponents[g])]},
largestComponentLength/numV
];

(* Overload Pattern *)
FractionInLargestComponent[amat_List]:=FractionInLargestComponent[AdjacencyGraph[amat]];

(* Catch-all Pattern *)
FractionInLargestComponent[args___]:=(Message[FractionInLargestComponent::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FractionInLargestKPathComponent*)


(* Primary Pattern *)
FractionInLargestKPathComponent[g_Graph,k_Integer]:=
(*
(*Given a graph, it computes the ratio of the number of vertices in the largest 
connected component to the total number of vertices in the graph.  *)*)
With[{kpcomps=KpathConnectedComponents[g,k]},
Max[Length/@kpcomps]/VertexCount[g]
];

(* Overload Pattern *)
FractionInLargestKPathComponent[g_Graph]:=With[{k=Length@@FindClique[g]},
Catch[If[k==1,Throw[1/VertexCount[g], "ECGravReturn$19"]];(*graph is fully isolated*)
FractionInLargestKPathComponent[g,k-1], "ECGravReturn$19"]
];

(* Overload Pattern *)
FractionInLargestKPathComponent[amat_List,k_Integer]:=FractionInLargestKPathComponent[AdjacencyGraph[amat],k];

(* Overload Pattern *)
FractionInLargestKPathComponent[amat_List]:=FractionInLargestKPathComponent[AdjacencyGraph[amat]];

(* Catch-all Pattern *)
FractionInLargestKPathComponent[args___]:=(Message[FractionInLargestKPathComponent::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueOrder*)


(* Primary Pattern *)
CliqueOrder[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},Length[clqs]];

(* Catch-all Pattern *)
CliqueOrder[args___]:=(Message[CliqueOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetOrder*)


(* Primary Pattern *)
FacetOrder[facetsLst_List]:=Length[facetsLst];

(* Catch-all Pattern *)
FacetOrder[args___]:=(Message[FacetOrder::argerr, args];
$Failed);


(* ::Item:: *)
(*EulerChi*)


(* Primary Pattern *)
EulerChi[facetsLst_List/;(Depth[facetsLst]==3&&Max[facetsLst]>1)]:=
(*(*****************************)
(* Last Updated: 04/02/2026  *)
(*****************************)*)

(*Given the list of facets of a simplicial complex (not necessarily 
clique complex), it finds the Euler Characteristic by counting all simplices.*)
Module[{connectedComponents=ConnectedComplexComponents[facetsLst],sortedfacetsLst,maxDim,altSignVector,fVector},

	Catch[If[Length[facetsLst]==1,Throw[1, "ECGravReturn$20"]];
	If[Length[connectedComponents]==1,
		sortedfacetsLst=Sort[#]&/@facetsLst;
		maxDim=Max[Length/@facetsLst];
		altSignVector=Table[(-1)^(i+1),{i,1,maxDim}];
		fVector=Table[Length[DeleteDuplicates[Join@@(Subsets[#,{i}]&/@sortedfacetsLst)]],{i,1,maxDim}];
		Throw[altSignVector . fVector, "ECGravReturn$20"],
		Throw[Total[EulerChi[#]&/@connectedComponents], "ECGravReturn$20"]
	], "ECGravReturn$20"]
];

(* Private helper *)
EulerChiAM[Amat_List]:=
(*Euler characteristic of the clique complex of the graph with adjacency matrix Amat, i.e.
chi = f1 - f2 + f3 - ... where fk is the number of k-cliques.

Rather than enumerate the cliques, the alternating sum is accumulated directly. Writing
	g[C] = Sum over cliques T contained in C, the empty one included, of (-1)^Length[T],
and splitting on whether a chosen vertex v of C is used, gives the recursion
	g[C] = g[C without v] - g[C intersect N[v]],      g[{}] = 1,
while chi = Sum over NONEMPTY cliques of (-1)^(Length[T]+1) = 1 - g[V]. So one recursion
over vertex subsets yields chi without ever listing a clique, and FindClique -- together
with the Graph object it requires -- is avoided entirely.

Vertex sets are integer bitmasks, and g is memoised on the Module-local symbol so the memo
is discarded when the call returns. The chosen vertex is cleared from the candidate set
before intersecting, so a self-loop on the diagonal cannot make it recur on itself; a
1x1 matrix therefore gives 1 whether or not its single entry is a self-loop.*)
Module[{n=Length[Amat],nbrMask,g},
	If[n==0,
		0,
		nbrMask=Table[FromDigits[Reverse[Amat[[v]]],2],{v,n}];
		g[0]=1;
		g[c_]:=g[c]=With[{p=BitLength[c]-1},
			g[BitClear[c,p]]-g[BitAnd[BitClear[c,p],nbrMask[[p+1]]]]];
		1-g[2^n-1]
	]
];

(* Overload Pattern *)
EulerChi[Amat_List/;(Depth[Amat]==3&&Max[Amat]<=1)]:=
(*Given adjacency matrix of a graph, it computes its Euler Characteristic by counting all
simplices. Read off the matrix alone -- no Graph object and no FindClique; see the
EulerChiAM helper above.*)
EulerChiAM[Amat];

(* Overload Pattern *)
EulerChi[g_Graph]:=
(*Given a graph it computes its Euler Characteristic by counting all simplices.*)
With[{maxCliques=FindClique[g,\[Infinity],All]},
EulerChi[maxCliques]
];

(* Catch-all Pattern *)
EulerChi[args___]:=(Message[EulerChi::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*RmsPurity*)


(* Primary Pattern *)
RmsPurity[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&
	Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1})]:=
(*Given a complex, it computes its mean purity = 1/|F| sum_{fi,fj}(dim(f_i)-dim(fj))^2, 
where fi, fj run over all facets. If there's just one facet it 
outputs 0. *)
If[Length[facetsLst]==1,0,StandardDeviation[Length/@facetsLst]]

(* Overload Pattern *)
RmsPurity[g_Graph]:=
(*Given a graph, it computes its mean purity = 1/|F| sum_{fi,fj}(dim(f_i)-dim(fj))^2, 
where fi, fj run over all maximal cliques. If there's just one maximal clique it 
outputs 0. *)

With[{clqsLst=FindClique[g,\[Infinity],All]},
	Catch[If[Length[clqsLst]==1,Throw[0, "ECGravReturn$22"]];StandardDeviation[Length/@clqsLst], "ECGravReturn$22"]
];

(* Overload Pattern *)
RmsPurity[Amat_List/;(SymmetricMatrixQ[Amat]&&Sort[DeleteDuplicates[Flatten[Amat]]]=={0,1})]:=
RmsPurity[AdjacencyGraph[Amat]];

(* Catch-all Pattern *)
RmsPurity[args___]:=(Message[RmsPurity::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Branchedness*)


(* Primary Pattern *)
Branchedness[g_Graph,n_Integer]:=
(*
(*Given a graph and a dimension n, it computes the rms value of the degree of 
the n facets, minus 3/2. This quantity will be zero for pseudo manifolds 
if n is the purity of the graph minus one since all codimension 1 
faces have degree 1 or 2. ratio of the number of vertices in the largest 
connected component to th etotal numbe rof vertices in the graph.  *)*)
Module[{clqsLst=FindClique[g,\[Infinity],All],maxDim,nfaces},
Catch[maxDim=Max[Length/@clqsLst];
nfaces=Union@@(Subsets[#,{n}]&/@clqsLst);
If[maxDim<=1,Throw[2, "ECGravReturn$23"]];(*a bunch of isolated vertices have degree 0 each so the output would be 2*)
Mean[((HyperDeg[g,#]&/@nfaces-3/2)^2-1/4)], "ECGravReturn$23"]
];

(* Overload Pattern *)
Branchedness[g_Graph]:=
With[{dMad=Length@@FindClique[g]},
	Branchedness[g,dMad-1]
];

(* Overload Pattern: adjacency-matrix input *)
Branchedness[amat_List/;(SymmetricMatrixQ[amat]&&SubsetQ[{0,1},DeleteDuplicates[Flatten[amat]]]),n_Integer]:=
	Branchedness[AdjacencyGraph[amat],n];

(* Overload Pattern: adjacency-matrix input *)
Branchedness[amat_List/;(SymmetricMatrixQ[amat]&&SubsetQ[{0,1},DeleteDuplicates[Flatten[amat]]])]:=
	Branchedness[AdjacencyGraph[amat]];

(* Overload Pattern: facet-list input, interpreted as a clique complex *)
Branchedness[facetsLst_List,n_Integer]:=Branchedness[GraphFromCliques[facetsLst],n];

(* Overload Pattern: facet-list input, interpreted as a clique complex *)
Branchedness[facetsLst_List]:=Branchedness[GraphFromCliques[facetsLst]];

(* Catch-all Pattern *)
Branchedness[args___]:=(Message[Branchedness::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Graph Dimensions*)


(* ::Item::Closed:: *)
(*AvgKDim*)


(* Primary Pattern *)
AvgKDim[g_Graph]:=
(*Computes the inductive (Knill) dimension of a graph as defined by Oliver Knill.*)
Block[{n,ne,isolatedVertices,ni,components,joinComponents,joinCount,cliques,temp,clqAssoc,numClqSizes,w,gm,degassoc,iassoc,rassoc,resultassoc,ComputeDeg,ComputeI,ComputeR,curClqs,curClqFaces,curI,curR,curval,curdeg,dimg},
Catch[n=VertexCount[g];
ne=EdgeCount[g];
If[n==0,Throw[-1, "ECGravReturn$24"]];
If[n==1||ne==0,Throw[0, "ECGravReturn$24"]];
If[ne==n (n-1)/2,Throw[n-1, "ECGravReturn$24"]];
If[ne==n (n-1)/2-1, Throw[n-2, "ECGravReturn$24"]];

isolatedVertices=Select[VertexList[g],VertexDegree[g,#]==0&];
ni=Length[isolatedVertices];
If[ni>0,
Throw[AvgKDim[Subgraph[g,Complement[VertexList[g],isolatedVertices]]]*(n-ni)/n, "ECGravReturn$24"]
];
If[TreeGraphQ[g],Throw[1, "ECGravReturn$24"]];

components=ConnectedGraphComponents[g];
If[Length[components]>1,Throw[Sum[VertexCount[components[[i]]]*AvgKDim[components[[i]]],{i,Length[components]}]/n, "ECGravReturn$24"]];

(************************************************************
Check if the graph is a join of subgraphs.
*******************************************************)
joinComponents=ConnectedGraphComponents[GraphComplement[g]];
joinCount=Length[joinComponents];

If[joinCount>1,Throw[joinCount-1+Sum[AvgKDim[Subgraph[g,VertexList[i]]],{i,joinComponents}], "ECGravReturn$24"]
];


cliques=FindClique[g,\[Infinity],All];
temp=Reap[Do[Sow[i,Length[i]],{i,cliques}]][[2]];
clqAssoc=<|Table[Length[i[[1]]]->i,{i,temp}]|>;
numClqSizes=Length[clqAssoc];
w=Max[Keys[clqAssoc]];


If[numClqSizes==1,Throw[w-1, "ECGravReturn$24"]];

degassoc=<||>;
iassoc=<||>;
rassoc=<||>;
resultassoc=<||>;

ComputeDeg[clq1_List]:=Block[{spheres,deg},
Catch[deg=Lookup[degassoc,Key[clq1],-1];
If[deg>=0,Throw[deg, "ECGravReturn$25"]];
spheres=Table[AdjacencyList[g,v],{v,clq1}];
(*Print["spheres",spheres]*);
deg=Length[Intersection@@spheres];
degassoc[clq1]=deg;
deg, "ECGravReturn$25"]
];

ComputeI[clq2_List]:=Block[{spheres,sph,isolatedNodes,isolatedNodesCt},
Catch[isolatedNodesCt=Lookup[iassoc,Key[clq2],-1];
If[isolatedNodesCt>=0,Throw[isolatedNodesCt, "ECGravReturn$26"]];
spheres=Table[AdjacencyList[g,v],{v,clq2}];
sph=Subgraph[g,Intersection@@spheres];
isolatedNodes=Select[VertexList[sph],VertexDegree[sph,#]==0&];
isolatedNodesCt=Length[isolatedNodes];
iassoc[clq2]=isolatedNodesCt;
isolatedNodesCt, "ECGravReturn$26"]
];

ComputeR[clq3_List]:=Block[{k,rval,deg,faces,faceRvals},
Catch[rval=Lookup[rassoc,Key[clq3],-1];
If[rval>=0,Throw[rval, "ECGravReturn$27"]];

k=Length[clq3];
If[k==1,rval=1/VertexDegree[g,clq3[[1]]];rassoc[clq3]=rval;Throw[rval, "ECGravReturn$27"]
];
deg=Lookup[degassoc,Key[clq3],ComputeDeg[clq3]];
faces=Subsets[clq3,{k-1}];
faceRvals=Table[Lookup[rassoc,Key[i],ComputeR[i]],{i,faces}];
rval=Total[faceRvals]/deg;
rassoc[clq3]=rval;
rval, "ECGravReturn$27"]
];

Do[
curClqs=clqAssoc[i];
curClqFaces=Union@@Table[Subsets[j,{i-1}],{j,curClqs}];
Do[
curI=ComputeI[k];
curR=ComputeR[k];
curval=(w+1-i)*curI*curR;
resultassoc[k]=curval;
,{k,curClqFaces}]
,{i,Keys[clqAssoc]}
];

dimg=w-(1/n)*Total[resultassoc];
dimg, "ECGravReturn$24"]
];

(* Overload Pattern *)
AvgKDim[Amat_List]:=With[{g=AdjacencyGraph[Amat]},AvgKDim[g]];

(* Catch-all Pattern *)
AvgKDim[args___]:=(Message[AvgKDim::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*SpectralDim*)


(* Primary Pattern *)
SpectralDim[Amat_List/;ArrayDepth[Amat]>1,s_Integer/;s>1]:=
(* Spectral dimension for all nodes at step s *)
Module[{n=Length[Amat],AmatPows= MatrixPower[Amat,s],totTable,Ps,avgP},

totTable=Table[Total[AmatPows[[i]]],{i,n}];

Ps=Table[If[totTable[[i]]==0,0,(AmatPows[[i,i]]/totTable[[i]])],{i,n}];
avgP=Total[Ps]/n;

-2.0*Log[avgP*1.0]/Log[s*1.0]
];

(* Overload Pattern *)
SpectralDim[g_Graph,s_Integer/;s>1]:=
(* Spectral dimension for all nodes at step s *)
With[{Amat=Normal[AdjacencyMatrix[g]]},
	SpectralDim[Amat,s]
];

(* Catch-all Pattern *)
SpectralDim[args___]:=(Message[SpectralDim::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HausdorffDim*)


(* Primary Pattern *)
HausdorffDim[Amat_List/;ArrayDepth[Amat]>1,Dmat_List, s_Integer/;s>1]:=
(* Hausdorff dimension for the first node at step s when the distance between nodes is 
given by the matrix Dmat*)
Module[{n=Length[Amat],vol,result},
vol=Sum[GVolume[Amat,Dmat,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,Dmat,s,1];(*Vertex 1 picked at random*)*)
result=Log[vol*1.0]/Log[s*1.0];
result
];

(* Overload Pattern *)
HausdorffDim[g_Graph,Dmat_List, s_Integer/;s>1]:=
(* Hausdorff dimension for the first node at step s when the distance between nodes is 
given by the matrix Dmat*)
Module[{n=VertexCount[g],vol,result},
vol=Sum[GVolume[g,Dmat,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,Dmat,s,1];(*Vertex 1 picked at random*)*)
result=Log[vol*1.0]/Log[s*1.0];
result
];

(* Overload Pattern *)
HausdorffDim[Amat_List/;ArrayDepth[Amat]>1,s_Integer/;s>1]:=
(* Hausdorff dimension for all nodes at step s *)
Module[{n=Length[Amat],vol},
n=Length[Amat];
vol=Sum[GVolume[Amat,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,s,1]; (*Vertex 1 picked at random*) *)
Log[vol*1.0]/Log[s*1.0]
];

(* Overload Pattern *)
HausdorffDim[g_Graph,s_Integer/;s>1]:=
(* Hausdorff dimension for all nodes at step s *)
Module[{n=VertexCount[g],vol},
vol=Sum[GVolume[g,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,s,1]; (*Vertex 1 picked at random*) *)
Log[vol*1.0]/Log[s*1.0]
];

(* Catch-all Pattern *)
HausdorffDim[args___]:=(Message[HausdorffDim::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Pure and Geometric Complexes*)


(* ::Item::Closed:: *)
(*ClosedDGraphQ*)


(* Primary Pattern *)
ClosedDGraphQ[g_Graph]:=
(*Outputs True if the clique complex of the graph g is a combinatorial manifold without
boundary (a closed d-graph): the unit sphere of every vertex is recursively a closed d-graph,
bottoming out at a single cycle (the 1-sphere) and at two isolated vertices (the 0-sphere).
This is the topological manifold condition only and does NOT test sphere-ness -- a graph whose
clique complex is a torus returns True. Use DSphereQ to test whether it is a sphere.*)
Module[{n=VertexCount[g],vlist=VertexList[g],clqs=FindClique[g,\[Infinity],All],numClqSizes,dim,spheres},
Catch[If[EdgeCount[g]==0,If[n==2,Throw[True, "ECGravReturn$28"],Throw[False, "ECGravReturn$28"]]];
numClqSizes=DeleteDuplicates[Length/@clqs];
If[Length[numClqSizes]>1,Throw[False, "ECGravReturn$28"]];
dim=numClqSizes[[1]];
If[dim==2,Throw[PathGraphQ[g]&&!TreeGraphQ[g], "ECGravReturn$28"]];
spheres=Table[Sph[g,i],{i,vlist}];
AllTrue[spheres,ClosedDGraphQ[#]&], "ECGravReturn$28"]
];

(* Catch-all Pattern *)
ClosedDGraphQ[args___]:=(Message[ClosedDGraphQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*DSphereQ*)


(* Primary Pattern *)
DSphereQ[g_Graph]:=
(*Outputs True if the clique complex of the graph g is a combinatorial manifold homeomorphic
to a sphere: a closed d-graph (see ClosedDGraphQ) whose clique complex additionally has the
Euler characteristic of a sphere of the same dimension, chi = 1 + (-1)^d (d = dimension of the
clique complex = one less than the largest clique). For surfaces this is exact -- only the
2-sphere returns True, while a graph whose clique complex is a torus returns False; in
dimension >= 3 the Euler-characteristic condition is necessary but not sufficient.*)
If[!ClosedDGraphQ[g],False,
	EulerChi[g]==1+(-1)^(Length[First[FindClique[g]]]-1)];

(* Catch-all Pattern *)
DSphereQ[args___]:=(Message[DSphereQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*DGraphQ*)


(* Primary Pattern *)
DGraphQ[g_Graph]:=
(*Returns true of the graph is a geometric graph and false if not.*)
Module[{n=VertexCount[g],vlist=VertexList[g],clqs=FindClique[g,\[Infinity],All],numClqSizes,dim,spheres},
Catch[If[EdgeCount[g]==0,If[n<=2,Throw[True, "ECGravReturn$29"],Throw[False, "ECGravReturn$29"]]];
If[CompleteGraphQ[g],Throw[True, "ECGravReturn$29"]];(*Complete graphs are d-graphs*)
numClqSizes=DeleteDuplicates[Length/@clqs];
If[Length[numClqSizes]>1,Throw[False, "ECGravReturn$29"]];
dim=numClqSizes[[1]];
If[dim==2,Throw[PathGraphQ[g], "ECGravReturn$29"]];
spheres=Table[Sph[g,i],{i,vlist}];
AllTrue[spheres,DGraphQ[#]&], "ECGravReturn$29"]
];

(* Catch-all Pattern *)
DGraphQ[args___]:=(Message[DGraphQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*DGraphBoundary*)


(* Primary Pattern *)
DGraphBoundary[g_Graph]:=
(*Outputs the induced subgraph which is the boundary of the graph g. The boudnary of a geometric graph is the induced subgraph over boundary nodes, i.e. those whose unit spheres are contractible, or paths. *)
Module[{vlist=VertexList[g],interiorPoints,bdryPoints},
(*If[!DGraphQ[g],Return["Error! The graph is not geometric."]];*)
Catch[interiorPoints=Select[vlist,ClosedDGraphQ[Sph[g,#]]&];(*interior vertices have a sphere link; ClosedDGraphQ preserves the original behaviour of this test (it was DSphereQ, which used to be the closed-manifold predicate)*)
bdryPoints=Complement[vlist,interiorPoints];
Throw[Subgraph[g,bdryPoints], "ECGravReturn$30"], "ECGravReturn$30"]
];

(* Overload Pattern *)
DGraphBoundary[Amat_List]:=
(*Outputs the induced subgraph which is the boundary of the graph with adjacency matrix 
Amat. The boudnary of a geometric graph is the induced subgraph over boundary nodes, i.e. those whose unit spheres are contractible, or paths. *)
With[{bdryGraph=DGraphBoundary[AdjacencyGraph[Amat]]},
	Catch[If[VertexCount[bdryGraph]==0,Throw[{}, "ECGravReturn$31"]];
	Normal[AdjacencyMatrix[bdryGraph]], "ECGravReturn$31"]
];

(* Catch-all Pattern *)
DGraphBoundary[args___]:=(Message[DGraphBoundary::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CombinatorialManifoldQ*)


(* Primary Pattern *)
CombinatorialManifoldQ[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it checks whether or not it is a combinatorial manifold. 
Returns true if a combinatorial manifold and false if not.*)
Module[{purity,codim1faces,links},
Catch[If[facelst=={},Throw[True, "ECGravReturn$32"]];(*The empty graph is the -1 sphere*)
If[Length[facelst]==1,Throw[True, "ECGravReturn$32"]];
(*an n-simplex is a combinatorial manifold, an n-ball *)
If[Length[DeleteDuplicates[Length/@facelst]]>1,Throw[False, "ECGravReturn$32"]];(*the complex has to be pure*)

purity=Length[facelst[[1]]];
If[purity==1,
If[Length[facelst]<=2,Throw[True, "ECGravReturn$32"],Throw[False, "ECGravReturn$32"]]
];
If[purity==2,
If[PathGraphQ[GraphFromCliques[facelst]],Throw[True, "ECGravReturn$32"],Throw[False, "ECGravReturn$32"]]
];
If[Length[ConnectedComplexComponents[facelst]]>1,Throw[False, "ECGravReturn$32"]];(*the complex has to be connected*)
codim1faces=Union@@(Subsets[#,{purity-1}]&/@facelst);
If[AnyTrue[codim1faces,HyperDeg[facelst,#]>2&],Throw[False, "ECGravReturn$32"]];
links=Lnk[facelst,{#}]&/@DeleteDuplicates[Flatten[facelst]];
AllTrue[links,CombinatorialManifoldQ[#]&], "ECGravReturn$32"]
];

(* Catch-all Pattern *)
CombinatorialManifoldQ[args___]:=(Message[CombinatorialManifoldQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ClosedCombinatorialManifoldQ*)


(* Primary Pattern *)
ClosedCombinatorialManifoldQ[facelst_List]:=
(*Given an input pure simplicial complex as a list of facets, checks whether it is a
combinatorial manifold without boundary (a closed combinatorial manifold): connected,
pure, every codimension-1 face contained in exactly two facets, and every vertex link
recursively a closed combinatorial manifold. Returns True or False. This tests the
topological manifold condition only and does NOT test sphere-ness -- a torus, Klein
bottle, or sphere all return True. Use CombinatorialSphereQ to test for a sphere.*)
Module[{purity,codim1faces,links},
Catch[If[facelst=={},Throw[True, "ECGravReturn$33"]];(*The empty complex is closed (the -1 sphere)*)
If[Length[facelst]==1,Throw[False, "ECGravReturn$33"]];
(*a single n-simplex is an n-ball: it has boundary, so it is not closed *)
If[Length[DeleteDuplicates[Length/@facelst]]>1,Throw[False, "ECGravReturn$33"]];(*the complex has to be pure*)

purity=Length[facelst[[1]]];

If[purity==1,
	If[Length[facelst]==2,Throw[True, "ECGravReturn$33"],Throw[False, "ECGravReturn$33"]]
];(*the only closed 0-manifold is S^0 = two points*)

If[purity==2,Throw[With[{g=GraphFromCliques[facelst]},
	If[PathGraphQ[g]&&!AcyclicGraphQ[g],True,False]], "ECGravReturn$33"]
];(*the only closed 1-manifold is a single cycle*)

If[Length[ConnectedComplexComponents[facelst]]>1,Throw[False, "ECGravReturn$33"]];(*the complex has to be connected*)

codim1faces=Union@@(Subsets[#,{purity-1}]&/@facelst);
If[AnyTrue[codim1faces,HyperDeg[facelst,#]!=2&],Throw[False, "ECGravReturn$33"]];(*closed: every ridge lies in exactly two facets*)

links=Lnk[facelst,{#}]&/@DeleteDuplicates[Flatten[facelst]];
AllTrue[links,ClosedCombinatorialManifoldQ[#]&], "ECGravReturn$33"]
];

(* Catch-all Pattern *)
ClosedCombinatorialManifoldQ[args___]:=(Message[ClosedCombinatorialManifoldQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CombinatorialSphereQ*)


(* Primary Pattern *)
CombinatorialSphereQ[facelst_List]:=
(*Given an input pure simplicial complex as a list of facets, checks whether it is a
combinatorial sphere. A combinatorial d-sphere is a closed combinatorial manifold
(see ClosedCombinatorialManifoldQ) that additionally has the Euler characteristic of the
d-sphere, chi = 1 + (-1)^d. For surfaces (dimension <= 2) this is exact: only the
2-sphere returns True, while a torus, Klein bottle, or projective plane return False.
In dimension >= 3 the Euler-characteristic condition is necessary but not sufficient
(recognizing PL spheres is not algorithmically decidable in general), so there the test
is a homology-level filter. Returns True or False.*)
Module[{purity},
Catch[If[facelst=={},Throw[True, "ECGravReturn$34"]];(*The empty complex is the -1 sphere*)
If[!ClosedCombinatorialManifoldQ[facelst],Throw[False, "ECGravReturn$34"]];(*a sphere is first of all a closed manifold*)
purity=Length[facelst[[1]]];
EulerChi[facelst]==1+(-1)^(purity-1), "ECGravReturn$34"]
];

(* Catch-all Pattern *)
CombinatorialSphereQ[args___]:=(Message[CombinatorialSphereQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*OrientableCombinatorialManifoldQ*)


(* Primary Pattern *)
OrientableCombinatorialManifoldQ[inputfacetsLst_List]:=
(*Checks whether or not a pseudo-combinatorial manifold given as a set of facets is 
orientable. Does not check whether or not the input is a combinatorial pseudomanifold *)
Module[{p=Length[inputfacetsLst[[1]]],q=Length[inputfacetsLst],facetsLst=Sort/@inputfacetsLst,orientations,RelativeOrientation,edgeList,graph,spanningtree,m1,m2,unexploredEdges},
Catch[If[q<=2,Throw[True, "ECGravReturn$35"]];
If[p<=2,Throw[True, "ECGravReturn$35"]];

RelativeOrientation[simplex1_List,simplex2_List]:=
(*gives the relative orientation between neighboring facets*)
Block[{sharedFace=Intersection[simplex1,simplex2],remaining1,pos1,remaining2,pos2},(*Find the vertex not in the shared face for each simplex*)
Catch[If[Length[sharedFace]!=Length[simplex1]-1,Throw[0, "ECGravReturn$36"]];(*relative orientation of non-adjacent is set to 0*)
remaining1=First@Complement[simplex1,sharedFace];
remaining2=First@Complement[simplex2,sharedFace];
pos1=Position[simplex1,remaining1][[1,1]];
pos2=Position[simplex2,remaining2][[1,1]];

((-1)^pos1)*Signature[Delete[simplex1,pos1]]*(-1)^pos2*Signature[Delete[simplex2,pos2]], "ECGravReturn$36"]
];

orientations=<|Table[i->i,{i,facetsLst}]|>;
edgeList=Join@@(Table[Table[If[Length[Intersection[facetsLst[[i]],facetsLst[[j]]]]==p-1,UndirectedEdge[facetsLst[[j]],facetsLst[[i]]],Nothing],{i,j+1,q}],{j,1,q-1}]);

graph=Graph[facetsLst,edgeList];

spanningtree=FindSpanningTree[graph];


(*Find Orientation on the spanning tree*)
DepthFirstScan[spanningtree,{"DiscoverVertex"->(({m1,m2}={#1,#2};

If[RelativeOrientation[orientations[[Key[m1]]],orientations[[Key[m2]]]]==+1,
orientations[[Key[m1]]]=SelectFirst[Permutations[orientations[[Key[m1]]]],RelativeOrientation[orientations[[Key[m2]]],#]==-1&]];
(*Print[" relative orientation after revision ",orientations];*))&)}];

(*Check Compatability on the remaining edges*)
unexploredEdges=Complement[Sort/@edgeList,Sort/@EdgeList[spanningtree],SameTest->(DeleteDuplicates[{Sort[#1], Sort[#2]}]==={Sort[#1]}&)];

(*Print["unexploredEdges ",unexploredEdges];
Print["checks on unexplored edges ",<|Table[i->RelativeOrientation[orientations[[Key[i[[1]]]]],orientations[[Key[i[[2]]]]]],{i,unexploredEdges}]|>];*)

AllTrue[unexploredEdges,RelativeOrientation[orientations[[Key[#[[1]]]]],orientations[[Key[#[[2]]]]]]==-1&], "ECGravReturn$35"]
];


(* Overload Pattern *)
OrientableCombinatorialManifoldQ[g_Graph]:=
(*Checks whether or not the clique complex of the graph g is an orientable. 
Does not check whether or not the input is a combinatorial pseudomanifold *)
Module[{facetsLst=Sort/@FindClique[g,\[Infinity],All],p,q,orientations,RelativeOrientation,edgeList,graph,spanningtree,m1,m2,unexploredEdges},
Catch[p=Length[facetsLst[[1]]];
q=Length[facetsLst];
If[q<=2,Throw[True, "ECGravReturn$37"]];
If[p<=2,Throw[True, "ECGravReturn$37"]];

RelativeOrientation[simplex1_List,simplex2_List]:=
(*gives the relative orientation between neighboring facets*)
Block[{sharedFace=Intersection[simplex1,simplex2],remaining1,pos1,remaining2,pos2},(*Find the vertex not in the shared face for each simplex*)
Catch[If[Length[sharedFace]!=Length[simplex1]-1,Throw[0, "ECGravReturn$38"]];(*relative orientation of non-adjacent is set to 0*)
remaining1=First@Complement[simplex1,sharedFace];
remaining2=First@Complement[simplex2,sharedFace];
pos1=Position[simplex1,remaining1][[1,1]];
pos2=Position[simplex2,remaining2][[1,1]];

((-1)^pos1)*Signature[Delete[simplex1,pos1]]*(-1)^pos2*Signature[Delete[simplex2,pos2]], "ECGravReturn$38"]
];

orientations=<|Table[i->i,{i,facetsLst}]|>;
edgeList=Join@@(Table[Table[If[Length[Intersection[facetsLst[[i]],facetsLst[[j]]]]==p-1,UndirectedEdge[facetsLst[[j]],facetsLst[[i]]],Nothing],{i,j+1,q}],{j,1,q-1}]);

graph=Graph[facetsLst,edgeList];

spanningtree=FindSpanningTree[graph];


(*Find Orientation on the spanning tree*)
DepthFirstScan[spanningtree,{"DiscoverVertex"->(({m1,m2}={#1,#2};

If[RelativeOrientation[orientations[[Key[m1]]],orientations[[Key[m2]]]]==+1,
orientations[[Key[m1]]]=SelectFirst[Permutations[orientations[[Key[m1]]]],RelativeOrientation[orientations[[Key[m2]]],#]==-1&]];
(*Print[" relative orientation after revision ",orientations];*))&)}];

(*Check Compatability on the remaining edges*)
unexploredEdges=Complement[Sort/@edgeList,Sort/@EdgeList[spanningtree],SameTest->(DeleteDuplicates[{Sort[#1], Sort[#2]}]==={Sort[#1]}&)];

(*Print["unexploredEdges ",unexploredEdges];
Print["checks on unexplored edges ",<|Table[i->RelativeOrientation[orientations[[Key[i[[1]]]]],orientations[[Key[i[[2]]]]]],{i,unexploredEdges}]|>];*)

AllTrue[unexploredEdges,RelativeOrientation[orientations[[Key[#[[1]]]]],orientations[[Key[#[[2]]]]]]==-1&], "ECGravReturn$37"]
];


(* Catch-all Pattern *)
OrientableCombinatorialManifoldQ[args___]:=(Message[OrientableCombinatorialManifoldQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CombinatorialBoundary*)


(* Primary Pattern *)
CombinatorialBoundary[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it returns the boundary. 
It first checks whether or not the manifold is a combinatorial manifold and returns 
$Failed if the complex is not a combinatorial manifold with a warning. *)
Module[{purity=Length[facelst[[1]]],codim1faces},

If[CombinatorialManifoldQ[facelst],
codim1faces=Union@@(Subsets[#,{purity-1}]&/@facelst);
Select[codim1faces,HyperDeg[facelst,#]==1&],
Message[CombinatorialBoundary::notmanifold];$Failed]

];

(* Catch-all Pattern *)
CombinatorialBoundary[args___]:=(Message[CombinatorialBoundary::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CountHoles*)


(* Primary Pattern *)
CountHoles[g_Graph,k_Integer]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a graph. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@FindClique[g,Infinity,All]},
ResourceFunction["BettiNumbers"][simplices][[k]]
];

(* Overload Pattern *)
CountHoles[g_Graph]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a graph. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@FindClique[g,Infinity,All]},
ResourceFunction["BettiNumbers"][simplices]
];

(* Overload Pattern *)
CountHoles[maxClqLst_List,k_Integer]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a simplicial complex. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@maxClqLst},
ResourceFunction["BettiNumbers"][simplices][[k]]
];

(* Overload Pattern *)
CountHoles[maxClqLst_List]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a simplicial complex. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@maxClqLst},
ResourceFunction["BettiNumbers"][simplices]
];

(* Catch-all Pattern *)
CountHoles[args___]:=(Message[CountHoles::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Visualize2DComplex*)


(* Primary Pattern *)
Visualize2DComplex[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it returns a 3D embedding. *)
Module[{facets,cmap,coord,scoords},

cmap=FindGraphIsomorphism[GraphFromCliques[facelst],CanonicalGraph[GraphFromCliques[facelst]]];

facets=facelst/.cmap[[1]];

coord=GraphEmbedding[CanonicalGraph[GraphFromCliques[facets]],"SpringEmbedding",3];
(*"SpringElectricalEmbedding","SpringEmbedding","HighDimensionalEmbedding","CircularEmbedding","SpiralEmbedding","RandomEmbedding","RadialEmbedding","SpectralEmbedding","StarEmbedding"}*)
(*Graphics3D[{Opacity[1.7],Yellow,GraphicsComplex[coord,Polygon[facets]]}];*)
scoords=Map[Part[coord,#]&,facets];
Graphics3D[{Opacity[0.9],Table[Simplex[i],{i,scoords}]}]

];

(* Catch-all Pattern *)
Visualize2DComplex[args___]:=(Message[Visualize2DComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Counting Pure Complexes*)


(* ::Item::Closed:: *)
(*RankComb*)


(* Primary Pattern *)

RankComb[set_List,numLabels_Integer]:=
(***************************************
Given a sorted set which is a subset of the set {0,1,,...,numLabels-1}, it assigns a 
unique integer to the set between 0 and numLabels choose length(set).
E.g., Rank[{0,1},3] = 0, and  Rank[{1,2},3]=3.
***************************************)
Module[{setSize=Length[set],result},
Catch[If[AnyTrue[set,#>=numLabels&],Throw[Null, "ECGravReturn$39"]];
result=Sum[Sum[Binomial[numLabels-k-j,setSize-k],{j,1,set[[k]]-k+1}],{k,1,setSize}];
result, "ECGravReturn$39"]
];

(* Catch-all Pattern *)
RankComb[args___]:=(Message[RankComb::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*UnrankComb*)


(* Primary Pattern *)
UnrankComb[l_Integer,numLabels_Integer,setSize_Integer]:=
(***************************************
Given an integer l between 0 and (numLabels choose setSize)-1, it assigns a unique sorted 
set which is a sorted setSize-subset of {0,1,,...,numLabels-1}.
E.g., UnrankComb[0,3,2] = {0,1}, and  UnrankComb[2,3,2] = {1,2}.
***************************************)
Module[{result={},lCur,iCur,iNext,lowVal,highVal,stopnum},
Catch[If[Binomial[numLabels,setSize]<=l,
Throw[Null, "ECGravReturn$40"];
];
lCur=l;
iCur=0;
iNext=1;

Do[

stopnum=0;
lowVal=Binomial[numLabels-iCur-1,setSize-k];
highVal=lowVal+Binomial[numLabels-iNext-1,setSize-k];

While[highVal<=lCur&&stopnum<10^10,
stopnum++;

iCur=iNext;
iNext++;

lowVal=highVal;
highVal+=Binomial[numLabels-iNext-1,setSize-k];

];


If[lowVal<=lCur,

AppendTo[result,iNext];
iCur=iNext+1;
iNext=iCur+1;
lCur=lCur-lowVal,

AppendTo[result,iCur];
iCur=iNext;
iNext=iCur+1;

];

,{k,1,setSize}
];

result, "ECGravReturn$40"]
];

(* Catch-all Pattern *)
UnrankComb[args___]:=(Message[UnrankComb::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*NumVertexLabeledPureComplexes*)


(* Private helper *)
NumPCKernel[p_Integer]:=NumPCKernel[p]=Reverse[Table[Binomial[p,k],{k,1,p}]];

(* Private helper *)
NumPCAdvance[p_Integer,q_Integer,prevRow_List]:=
(*One row of the table: {N(p,q,n) : n=p,...,p q}, built from the q-1 row. Writing C(n) for
Binomial[n,p], the recursion of the primary pattern below reads, term by term,
	q N(p,q,n) = (C(n)-(q-1)) N(p,q-1,n) + C(n) Sum[Binomial[p,k] N(p,q-1,n-k),{k,1,p}],
so a whole row is one vector expression. The inner sum is a fixed-kernel correlation of the
previous row against {Binomial[p,p],...,Binomial[p,1]}, hence the p leading zeros. Every N is
a count, so the right-hand side is divisible by q entrywise and the arithmetic stays in exact
integers -- dividing once at the end avoids building a Rational, and its GCD, per term.
Takes the previous row as an argument rather than reading the cache, so NumPCRow can use it
both to fill the cache and to run past the cap without storing anything.*)
Module[{len=p*(q-1)+1,prev,cvec,conv},
	prev=PadRight[prevRow,len];
	cvec=Table[Binomial[n,p],{n,p,p*q}];
	conv=Take[ListCorrelate[NumPCKernel[p],Join[ConstantArray[0,p],prev]],len];
	((cvec-(q-1))*prev+cvec*conv)/q
];

(* Private helper *)
NumPCStep[p_Integer,1]:=NumPCStep[p,1]={1};

(* Private helper *)
NumPCStep[p_Integer,q_Integer]:=NumPCStep[p,q]=NumPCAdvance[p,q,NumPCStep[p,q-1]];

(* Private helper *)
NumPCTop[p_Integer]:=0;

(* Private helper: highest facet order whose rows are kept in the cache. -----------------
   A row for (p,q) is p(q-1)+1 entries wide and the entries themselves grow with q, so the
   cache costs on the order of p q^4 and has to stop somewhere. Measured cumulative size at
   this cap: 2.5 MB at p=2, 6.3 MB at p=3, 11.8 MB at p=4, 28 MB at p=6, 51 MB at p=8. That
   is far above any facet order the samplers in this package actually use.
   Raise it, or set it to 0 to switch caching off, with
	ECGrav`Private`$NumPCMaxCachedQ = 400;
   and reclaim the memory with ECGrav`Private`NumPCClearCache[]. *)
$NumPCMaxCachedQ=150;

(* Private helper *)
NumPCRow[p_Integer,q_Integer]:=
(*The row for (p,q), cached up to $NumPCMaxCachedQ. Rows are filled bottom-up by an explicit
loop rather than by letting NumPCStep recurse on itself, so the recursion depth stays at 1
instead of q: the older recursive form ran at q=300 but died on
TerminatedEvaluation[RecursionLimit] by q=700.

NumPCTop[p] is the highest row already built for this p, so a later call pays only for the
rows it actually adds, in any order of q. The rows persist for the session -- that is the
point, since the callers ask for many (q,n) at one fixed p.

Past the cap the tail is advanced transiently and thrown away, which bounds the cache without
losing the vectorised row recursion. It still seeds from the deepest cached row at or below q,
so the cost past the cap is the tail alone, not a rebuild from q=1, and it is paid again on
every call. Seeding is capped at q rather than at NumPCTop[p]: lowering $NumPCMaxCachedQ in a
session that has already cached past the new value would otherwise leave NumPCTop[p] above q
and return that higher row outright, with the advancing loop running empty.*)
Module[{cap=$NumPCMaxCachedQ,target=Min[q,$NumPCMaxCachedQ],top=NumPCTop[p],row,q0},
	If[target>top,
		Do[NumPCStep[p,k],{k,top+1,target}];
		NumPCTop[p]=target];
	If[q<=cap,
		NumPCStep[p,q],
		q0=Min[NumPCTop[p],q];
		row=If[q0>=1,NumPCStep[p,q0],{1}];
		q0=Max[q0,1];
		Do[row=NumPCAdvance[p,k,row],{k,q0+1,q}];
		row
	]
];

(* Private helper *)
NumPCClearCache[]:=
(*Drops everything memoised behind NumVertexLabeledPureComplexes and NumUnlabeledPureComplexes, and
behind the composition weights that
RandomVertexLabeledPureSimplicialComplex draws from, keeping only the definitions themselves --
the memoised entries are the DownValues whose left-hand side carries no Pattern. Returns the
number of bytes reclaimed. Nothing else in the package holds on to either, so calling this is
always safe; it only costs the rebuilding.*)
Module[{before,after,patterned=!FreeQ[First[#],Pattern]&,
		cached={NumPCStep,NumPCTop,NumPCKernel,RandomPureComplexCDF,
			RandFLPCWeightCounts,RandFLPCCompletions,RandFLPCTypeWeights,NumULPCA,
			RandULPCOrbits,RandULPCNd,RandULPCSubMulti,RandULPCFixCov,RandULPCTypeWeights,
			RandULPCProfiles,RandULPCCCov,RandULPCComps}},
	before=Total[ByteCount/@Map[DownValues,cached]];
	Scan[(DownValues[#]=Select[DownValues[#],patterned])&,cached];
	after=Total[ByteCount/@Map[DownValues,cached]];
	before-after
];

(* Private helper *)
NumPCRec[pp_Integer,qq_Integer,nn_Integer]:=
(*The original per-call recursion, kept only for degenerate purity (p<1) with q>=1, where the
row picture does not apply; it reproduces those boundary values exactly as before. Memoised on
the Block-local symbol, so the memo is discarded when the call returns.*)
Block[{NumPCRec},

NumPCRec[p_Integer,q_Integer,n_Integer]:=NumPCRec[p,q,n]=
Which[q<0,0,q==1&&n==p,1,q==1&&n!=p,0,True,
(((Binomial[n,p]-(q-1))/q)*NumPCRec[p,q-1,n]
	+(Binomial[n,p]/q)*Sum[Binomial[p,k]*NumPCRec[p,q-1,n-k],{k,1,p}])];

NumPCRec[pp,qq,nn]

];

(* Primary Pattern *)
NumVertexLabeledPureComplexes[pp_Integer,MM_Integer,nn_Integer]:=
(*Gives the number of vertex labeled pure simplicial complexes of purity p, facet order M,
and number of vertices n. Computes bottom-up over M and memoizes a whole row of n at a time;
see NumPCRow above. N(p,M,n) vanishes unless p<=n<=p M -- fewer than p vertices cannot carry
a facet, and M facets of p vertices cannot cover more than p M of them -- so those n are
answered without building any row.

M==0 is the empty complex: it covers no vertices, so N(p,0,n) is 1 at n==0 and 0 otherwise, for
every p. Feeding that into the recursion reproduces N(p,1,n)=KroneckerDelta[n,p], so the M==0
row is the consistent seed one step below the M==1 row the table actually starts from. (Earlier
versions had no M==0 branch at all and divided by M, returning Indeterminate.)*)
Which[
	MM<0,0,
	MM==0,If[nn==0,1,0],
	pp<1,NumPCRec[pp,MM,nn],
	nn<pp||nn>pp*MM,0,
	True,NumPCRow[pp,MM][[nn-pp+1]]
];

(* Overload Pattern *)
NumVertexLabeledPureComplexes[p_Integer,M_Integer]:=
(*The number of pure simplicial complexes of purity p and facet order M, summed over the vertex
count n. n0 is the first n carrying enough distinct p-subsets to supply M facets; below it every
term is zero. When no n in p..p M qualifies -- which happens only for degenerate p or M -- there
is nothing to sum and the answer is 0; earlier versions used the resulting n0=Null as a Table
iterator bound and returned the expression unevaluated.*)
Which[
	M<0,0,
	M==0,1,
	True,
	With[{n0=Catch[Do[If[Binomial[nn,p]>=M,Throw[nn]],{nn,p,p*M}]]},
		If[IntegerQ[n0],Total[Table[NumVertexLabeledPureComplexes[p,M,n],{n,n0,p*M}]],0]
	]
];

(* Catch-all Pattern *)
NumVertexLabeledPureComplexes[args___]:=(Message[NumVertexLabeledPureComplexes::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*NumPureComplexes*)


(* Alias Pattern *)
NumPureComplexes[args___]:=
(*Former name of NumVertexLabeledPureComplexes. The rename did not touch the argument order --
it is still (purity, facet order, vertex count), the facet order simply written M rather than q
-- so forwarding is a straight hand-off and old calls return exactly what they always did. Kept
deliberately rather than removed, so existing notebooks keep working.*)
NumVertexLabeledPureComplexes[args];


(* ::Item::Closed:: *)
(*NumFacetLabeledPureComplexes*)


(* Counting facet-labeled complexes. -------------------------------------------------------
   A facet-labeled pure complex is M labelled, pairwise distinct p-subsets of an n-element
   vertex set covering it, with the vertices NOT labelled. Transposing the incidence matrix
   turns that into an incidence tableau: a multiset of n nonempty subsets of [M] -- one row per
   vertex, listing the facets it lies in -- in which every label occurs in exactly p rows. The
   rows are a multiset because the vertices are unlabelled; distinct facets means the M columns
   are pairwise distinct, which is the "separating" condition.

   Counting them directly over n rows is hopeless, so the count goes through Burnside on the n
   row slots. Allowing empty rows first, let A(n) be the number of multisets of exactly n
   subsets of [M] with every label in exactly p of them; padding with empty rows gives
   A(n) = Sum_{k<=n} F(k), hence F(n) = A(n) - A(n-1). A tuple fixed by a permutation sigma is
   constant on each cycle, so choosing one amounts to choosing a subset per cycle, and the
   condition "label i lies in exactly p rows" becomes "the cycles containing i have lengths
   summing to p" -- separately for each label. The labels therefore decouple and

	|Fix(sigma)| = N(lambda,p)^M,   N(lambda,p) = [x^p] Product_k (1+x^k)^m_k,

   with lambda the cycle type of sigma and m_k its number of parts of length k. Cycles longer
   than p can never be used, so N depends only on m_1..m_p, leaving a sum over those few
   exponents times the number of ways to arrange the remaining long cycles.

   The separating condition is one substitution away: N^M counts arbitrary labellings, while
   labellings with pairwise DISTINCT columns are counted by the falling factorial N^(M). Since
   N^M = Sum_k StirlingS2[M,k] N^(k), swapping the weight IS the Stirling inversion of
   B(p,n,M) = Sum_k StirlingS2[M,k] F(p,n,k), applied cycle type by cycle type -- and unlike
   that inversion it has no cancellation, since every summand stays non-negative.

   Ported from the author's Facet-labeled-count.nb (sFDirect), whose values were verified there
   against a brute-force recursion. Cost is essentially independent of M, which enters only as
   an exponent, and grows like n^p. *)

(* Private helper *)
NumFLPCLongCycleTable[p_Integer,smax_Integer]:=
(*Entry s+1 is the number of permutations of [s] whose every cycle is longer than p. Read off
the cycle containing the point 1: it has some length k>p and can be laid down in (s-1)!/(s-k)!
ways, leaving the same problem on the other s-k points.*)
Module[{gg=ConstantArray[0,smax+1],tot,fall},
	gg[[1]]=1;
	Do[
		tot=0;
		fall=1;(*fall = (s-1)!/(s-k)!*)
		Do[
			If[k>p,tot+=fall*gg[[s-k+1]]];
			fall*=(s-k)
		,{k,1,s}];
		gg[[s+1]]=tot
	,{s,1,smax}];
	gg
];

(* Private helper *)
NumFLPCNCoeff[mvec_List,jvecs_List]:=
(*N(lambda,p) = [x^p] Product_k (1+x^k)^m_k, the number of sub-multisets of the cycle lengths
summing to p. Expanded over the partitions of p, jvecs holding each partition's part counts.*)
Sum[Times@@Binomial[mvec,jv],{jv,jvecs}];

(* Private helper *)
NumFLPCTSum[p_Integer,n_Integer,wf_,gg_List,fac_List,jvecs_List]:=
(*n! A(n), summed over the cycle types of the n row slots. The recursion walks the multiplicities
m_1..m_p of the short cycles; s is what is left for cycles longer than p, arranged gg[[s+1]]
ways, and the multinomial fac[[n+1]]/(fac[[s+1]] denom) counts the ways to lay the short cycles
down. wf turns the per-cycle-type count N into its weight: #^M would count all tableaux, while
FactorialPower[#,M] counts only those with distinct columns, which is what separating means.*)
If[n<0,
	0,
	Module[{total=0,rec},
		rec[k_,used_,denom_,mvec_]:=
			If[k>p,
				Module[{s=n-used,nv=NumFLPCNCoeff[mvec,jvecs],w},
					If[nv=!=0,
						w=wf[nv];
						If[w=!=0,
							total+=w*(fac[[n+1]]/(fac[[s+1]]*denom))*gg[[s+1]]]]],
				Do[
					rec[k+1,used+k*mk,denom*k^mk*mk!,Append[mvec,mk]]
				,{mk,0,Quotient[n-used,k]}]];
		rec[1,0,1,{}];
		total
	]
];

(* Private helper *)
NumFLPCCount[p_Integer,MM_Integer,n_Integer]:=
(*F(p,M,n) = (T(n) - n T(n-1))/n!, the empty-row padding differenced away.*)
Module[{gg,fac,jvecs,wf,tn,tn1},
	gg=NumFLPCLongCycleTable[p,n];
	fac=Table[k!,{k,0,n}];
	jvecs=Table[Count[nu,k],{nu,IntegerPartitions[p]},{k,1,p}];
	wf=FactorialPower[#,MM]&;
	tn=NumFLPCTSum[p,n,wf,gg,fac,jvecs];
	tn1=If[n>=1,NumFLPCTSum[p,n-1,wf,gg,fac,jvecs],0];
	(tn-n*tn1)/fac[[n+1]]
];

(* Primary Pattern *)
NumFacetLabeledPureComplexes[p_Integer,MM_Integer,n_Integer]:=
(*Gives the number of facet labeled pure simplicial complexes of purity p, facet order M and
number of vertices n. The argument order matches NumVertexLabeledPureComplexes: purity, facet
order, vertex count. Zero outside p<=n<=p M -- fewer than p vertices cannot carry a facet, and M
facets of p vertices cannot cover more than p M of them -- and zero when Binomial[n,p]<M, since
the M facets must be distinct p-subsets. M==0 is the empty complex, counted at n==0 as for
NumVertexLabeledPureComplexes.*)
Which[
	p<0||MM<0,0,
	MM==0,If[n==0,1,0],
	n<0||n<p||n>p*MM,0,
	Binomial[n,p]<MM,0,
	True,NumFLPCCount[p,MM,n]
];

(* Overload Pattern *)
NumFacetLabeledPureComplexes[p_Integer,MM_Integer]:=
(*The number of facet labeled pure complexes of purity p and facet order M, summed over the
vertex count n, mirroring the two-argument NumVertexLabeledPureComplexes.*)
Which[
	MM<0,0,
	MM==0,1,
	True,
	With[{n0=Catch[Do[If[Binomial[nn,p]>=MM,Throw[nn]],{nn,p,p*MM}]]},
		If[IntegerQ[n0],Total[Table[NumFacetLabeledPureComplexes[p,MM,n],{n,n0,p*MM}]],0]
	]
];

(* Catch-all Pattern *)
NumFacetLabeledPureComplexes[args___]:=(Message[NumFacetLabeledPureComplexes::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*NumUnlabeledPureComplexes*)


(* Counting fully unlabeled complexes. ------------------------------------------------------
   Neither the vertices nor the facets carry labels, so the objects are the isomorphism classes
   of pure complexes: M pairwise distinct p-subsets of an n-element vertex set covering it, up to
   relabelling the vertices AND permuting the facets.

   Acting on the ordered M-tuples that the facet-labeled count uses, the group is the product
   S_n x S_M, and (sigma,tau) fixes (F_1,...,F_M) iff sigma(F_i) = F_tau(i). The S_M factor
   collapses, though, and it is worth saying why rather than carrying it through a double
   Burnside: forgetting the order maps the tuples onto the SETS of M facets, the map commutes
   with S_n, and because the facets in a tuple are pairwise distinct S_M acts simply transitively
   on each fibre. So the S_n x S_M orbits on tuples ARE the S_n orbits on sets, and one Burnside
   over S_n alone suffices, with

	|Fix(sigma)| = #{sigma-invariant M-element sets of distinct p-subsets covering [n]}.

   That is why the tau-cycle bookkeeping this looks like it needs -- a tau-cycle of length l
   pinning one facet F with sigma^l(F) = F, cycle types of sigma^l, gcd(k,l) -- never appears.

   A sigma-invariant set is a union of <sigma>-orbits of p-subsets, so it is chosen by choosing
   which orbits to include. If sigma has n_d orbits of size exactly d on the p-subsets, then
   dropping the covering condition,

	#{sigma-invariant M-element sets} = [z^M] Product_d (1 + z^d)^n_d.

   Covering comes back by the padding trick the facet-labeled count also uses. Let A(n) count
   those sets up to relabelling with covering NOT required; an unlabeled object whose facets span
   j vertices is an unlabeled covering object on j vertices plus n-j indistinguishable isolated
   vertices, so A(n) = Sum_{j<=n} U(j) and U(n) = A(n) - A(n-1).

   The n_d themselves are never formed. They would come from Moebius inversion over the divisor
   lattice -- f(e), the number of p-subsets fixed by sigma^e, counts those of period dividing e,
   so n_d = (1/d) Sum_{e|d} mu(d/e) f(e) -- but taking the logarithmic derivative of the product
   above turns the whole thing into a recurrence whose coefficients are the f directly:

	m c(m) = Sum_{k=1..m} l(k) c(m-k),   c(0) = 1,   c(m) = [z^m] Product_d (1 + z^d)^n_d,
	l(k) = f(k) for k odd,   f(k) - 2 f(k/2) for k even.

   That is Newton's identity: log Product_d (1+z^d)^n_d = Sum_m (z^m/m) Sum_{d|m} d n_d
   (-1)^(m/d-1), and Sum_{d|m} d n_d is f(m), since the orbits of size dividing m are exactly the
   points fixed by sigma^m; splitting the sign by the parity of m/d leaves l. So the Moebius
   inversion, and with it the binomial expansion of (1+z^d)^n_d over enormous n_d, both drop out.
   f(k) is needed for every k <= M rather than only the divisors, but f(k) = f(gcd(k,LCM)) since
   sigma^k and sigma^gcd(k,LCM) fix the same subsets, so the evaluations are the same handful.

   sigma^e splits each k-cycle of sigma into gcd(k,e) cycles of length k/gcd(k,e), a fixed subset
   being a union of those:

	f(e) = [x^p] Product_k (1 + x^(k/gcd(k,e)))^(gcd(k,e) m_k).

   Unlike the facet-labeled count there is no restriction to partitions with parts <= p -- with
   the facets permutable a long cycle CAN be covered, by a facet whose sigma-orbit walks around
   it, as the four edges of a square are covered by a 4-cycle of vertices. Every cycle type of
   S_n therefore contributes and the cost grows like PartitionsP[n], which is the limiting factor
   here rather than n or M individually. Only divisors d <= M of LCM matter to f, since a factor
   (1 + z^d)^n_d with d > M contributes its constant term and nothing else below z^(M+1). *)

(* Private helper *)
NumULPCPolyMul[a_List,b_List,deg_Integer]:=
(*Product of two coefficient lists, truncated at degree deg.*)
Take[PadRight[ListConvolve[a,b,{1,-1},0],deg+1],deg+1];

(* Private helper *)
NumULPCPowPoly[r_Integer,c_Integer,deg_Integer]:=
(*(1+z^r)^c as a coefficient list truncated at degree deg. Only the multiples of r are stocked,
so this is the sparse form of the binomial expansion rather than a repeated multiplication.*)
Module[{res=ConstantArray[0,deg+1]},
	Do[res[[r*j+1]]=Binomial[c,j],{j,0,Quotient[deg,r]}];
	res
];

(* Private helper *)
NumULPCFixedSubsets[parts_List,e_Integer,p_Integer]:=
(*f(e): the number of p-subsets of [n] fixed by sigma^e, with parts the {length,multiplicity}
pairs of sigma's cycle type. A k-cycle becomes gcd(k,e) cycles of length k/gcd(k,e) under
sigma^e, and a fixed subset is a union of them. Cycles of sigma^e longer than p can never be
used, so those factors are 1 up to degree p and are skipped outright.*)
Module[{poly=PadRight[{1},p+1],g,r},
	Do[
		g=GCD[part[[1]],e];
		r=part[[1]]/g;
		If[r<=p,poly=NumULPCPolyMul[poly,NumULPCPowPoly[r,g*part[[2]],p],p]]
	,{part,parts}];
	poly[[p+1]]
];

(* Private helper *)
NumULPCFixSets[parts_List,p_Integer,MM_Integer]:=
(*|Fix(sigma)| on the M-element sets of distinct p-subsets, covering NOT imposed, by the Newton
recurrence of the header. f is evaluated only at the divisors of the cycle-length LCM that are at
most MM -- tested by divisibility over Range[MM] rather than by Divisors[LCM], since the LCM can
have far more divisors than that -- and every f(k) the recurrence asks for is one of them, because
sigma^k fixes the same subsets as sigma^gcd(k,LCM).

The recurrence carries the whole row c(0..MM) whether or not the caller wants it. That is not
waste: it is the same O(M^2) either way, and it is what a caller wanting all facet orders at one
vertex count would need.*)
Module[{lcmParts,f,ell,c},
	lcmParts=If[parts==={},1,LCM@@parts[[All,1]]];
	f=Association[Table[e->NumULPCFixedSubsets[parts,e,p],
		{e,Select[Range[MM],Divisible[lcmParts,#]&]}]];
	ell=Table[
		If[OddQ[k],
			f[GCD[k,lcmParts]],
			f[GCD[k,lcmParts]]-2*f[GCD[k/2,lcmParts]]]
	,{k,1,MM}];
	c=ConstantArray[0,MM+1];
	c[[1]]=1;
	Do[c[[m+1]]=(Take[ell,m] . Reverse[Take[c,m]])/m,{m,1,MM}];
	c[[MM+1]]
];

(* Private helper *)
NumULPCA[p_Integer,MM_Integer,n_Integer]:=NumULPCA[p,MM,n]=
(*A(p,M,n): unlabeled M-element sets of distinct p-subsets of [n] with covering NOT required, by
Burnside over S_n. n!/z_lambda is the number of permutations of cycle type lambda. Memoised for
the session and released by the shared NumPCClearCache[], so a scan over n pays for each row
once rather than twice through the differencing.*)
If[n<0,
	0,
	Total[Table[
		With[{parts=Tally[lambda]},
			(n!/(Times@@(#[[1]]^#[[2]]*#[[2]]!&/@parts)))*NumULPCFixSets[parts,p,MM]]
	,{lambda,IntegerPartitions[n]}]]/n!
];

(* Private helper *)
NumULPCCount[p_Integer,MM_Integer,n_Integer]:=
(*U(p,M,n) = A(n) - A(n-1), the isolated-vertex padding differenced away.*)
NumULPCA[p,MM,n]-NumULPCA[p,MM,n-1];

(* Primary Pattern *)
NumUnlabeledPureComplexes[p_Integer,MM_Integer,n_Integer]:=
(*Gives the number of fully unlabeled pure simplicial complexes of purity p, facet order M and
number of vertices n -- the isomorphism classes, with neither the vertices nor the facets
labelled. The argument order matches NumVertexLabeledPureComplexes and
NumFacetLabeledPureComplexes: purity, facet order, vertex count. The guards are theirs as well,
and are not merely an optimisation here: the count sums over the cycle types of S_n, so falling
through to it on an out-of-range n would call IntegerPartitions on it.*)
Which[
	p<0||MM<0,0,
	MM==0,If[n==0,1,0],
	n<0||n<p||n>p*MM,0,
	Binomial[n,p]<MM,0,
	True,NumULPCCount[p,MM,n]
];

(* Overload Pattern *)
NumUnlabeledPureComplexes[p_Integer,MM_Integer]:=
(*The number of unlabeled pure complexes of purity p and facet order M, summed over the vertex
count n, mirroring the two-argument NumVertexLabeledPureComplexes. The sum telescopes: U(n) is
zero past n = p M, since M facets of p vertices cover no more than that, so summing A(n)-A(n-1)
over all n leaves A(p,M,p M) alone. Note that this is the whole cost of the three-argument form
at its largest n, so it is only practical while PartitionsP[p M] is.*)
Which[
	p<0||MM<0,0,
	MM==0,1,
	True,NumULPCA[p,MM,p*MM]
];

(* Catch-all Pattern *)
NumUnlabeledPureComplexes[args___]:=(Message[NumUnlabeledPureComplexes::argerr, args];
$Failed);


(* ::Section::Closed:: *)
(*Choosing Isomorphism Classes*)


(* ::Subsection::Closed:: *)
(*Graph Isomorphism Classes*)


(* ::Item::Closed:: *)
(*GraphIsContained*)


(* :Code Section: *)

(* Primary Pattern *)
GraphIsContained[glist_List/;GraphQ[glist[[1]]],g_Graph/;GraphQ[g]]:=
(*
(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)
(*Notes: *)
(*Given a list of graphs glist and a new graph g, it returns True if there is a graph 
isomorphic to g in glist, False otherwise. It stops checking at the first occurence 
of isimorphic graph.*)
If[VertexCount[g]==0,True,
AnyTrue[glist,IsomorphicGraphQ[g,#]&]
];

(* Overload Pattern *)
GraphIsContained[Amats_List/;SquareMatrixQ[Amats[[1]]],Am_List/;SquareMatrixQ[Am]]:=(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given a list of adjacency matrices  Amats and a new adjacency matrix Am, it returns True if there is a graph isomorphic to Am in glist, False otherwise. It stops checking at the first occurence of isimorphic graph.*)If[Length[Am[[1]]]==0,True,AnyTrue[Amats,IsomorphicGraphQ[AdjacencyGraph[Am],AdjacencyGraph[#]]&]
];

(* Catch-all Pattern *)
GraphIsContained[args___]:=(Message[GraphIsContained::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicGraphs*)


(* Primary Pattern *)
SetAttributes[ChooseNonIsomorphicGraphs,Flat];
ChooseNonIsomorphicGraphs[li__List/;Length[{li}]>1&&GraphQ[{li}[[1,1]]]]:=
(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given an arbitrary number of lists of graphs, it merges them into one list of graphs where no two graphs belonging to different input lists are isomorphic. This method does not check if the graphs WITHIN a given list are non-isomorphic.*)
Block[{result={li}[[1]],length=Length[{li}]},
Do[
result=Join[result,Select[{li}[[n]],!GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Overload Pattern *)
ChooseNonIsomorphicGraphs[li_List/;GraphQ[li[[1]]]]:=
(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*given a single list of graphs, it generates a new list composed of non-isomorphic graphs.*)Block[{result={},remaining=li},
Reap[
While[Length[remaining]>0,
Sow[First[remaining]];
remaining=Select[remaining,IsomorphicGraphQ[First[remaining],#]==False&];
];
][[2,1]]
];

(* Overload Pattern *)
ChooseNonIsomorphicGraphs[li_List/;GraphQ[li[[1,1]]]]:=
(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given a list of lists of graphs, it merges them into one list where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)Block[{result=li[[1]],length=Length[li]},
Do[
result=Join[result,Select[li[[n]],!GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Overload Pattern *)
ChooseNonIsomorphicGraphs[li__List/;Length[{li}]>1&&SquareMatrixQ[{li}[[1,1]]]]:=(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given an arbitrary number of lists of adjacency matrices of graphs, it merges them into one list of graphs where no two graphs belonging to different input lists are isomorphic. This method does not check if the graphs WITHIN a given list are non-isomorphic.*)Block[{result={li}[[1]],length=Length[{li}]},
Do[
result=Join[result,Select[{li}[[n]],!GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Overload Pattern *)
ChooseNonIsomorphicGraphs[li_List/;MatrixQ[li[[1]]]]:=

(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*given a single list of adjacency matrices of graphs, it generates a new list of adjacency matrices of the largest set of mutually non-isomorphic adjacency matrices.*)Block[{result={},remaining=li},
Reap[
While[Length[remaining]>0,
Sow[First[remaining]];
remaining=Select[remaining,IsomorphicGraphQ[AdjacencyGraph[First[remaining]],AdjacencyGraph[#]]==False&];
];
][[2,1]]
];

(* Overload Pattern *)
ChooseNonIsomorphicGraphs[li_List/;SquareMatrixQ[li[[1,1]]]]:=(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given a list of lists of adjacency matrices, it merges them into one list where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)
Block[{result=li[[1]],length=Length[li]},
Do[
result=Join[result,Select[li[[n]],!GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Catch-all Pattern *)
ChooseNonIsomorphicGraphs[args___]:=(Message[ChooseNonIsomorphicGraphs::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Clique Complex Isomorphism Classes*)


(* ::Item::Closed:: *)
(*IsContainedClqComp*)


(* :Code Section: *)

(* Primary Pattern *)
IsContainedClqComp[cmpxesLst_List/;Depth[cmpxesLst]==4,cmpx_List/;(Depth[cmpx]==3||cmpx=={})]:=
(****************************)
(*Last updated: 09/13/2023  *)
(****************************)(*Given a list of clique complexes cmpxesLst (each given as a list of maximal cliques)  and a single other clique complex cmpx, it returns True if there is a graph isomorphic to clq in clqsLst, False otherwise. It stops checking at the first occurence of isimorphic graph.*)
If[Length[cmpx]==0||cmpx=={{}},True,AnyTrue[cmpxesLst,IsomorphicGraphQ[GraphFromCliques[cmpx],GraphFromCliques[#]]&]
];

(* Catch-all Pattern *)
IsContainedClqComp[args___]:=(Message[IsContainedClqComp::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicClqComplexes*)


(* Primary Pattern *)
SetAttributes[ChooseNonIsomorphicClqComplexes,Flat];
ChooseNonIsomorphicClqComplexes[li__List/;(Depth[{li}]==5&&Length[{li}]>1)]:=
(*
(*****************************)
(* Last Updated: 09/13/2023  *)
(*****************************)*)

(*Given an arbitrary sequence of lists of clique complexes (each clique complex given
  as lists of maximal cliques), it merges them into one list of clique complexes 
  where no two complexes belonging to different input lists are isomorphic. 
  This method does not check if the graphs WITHIN a given list are non-isomorphic.*)

Block[{result={li}[[1]],length=Length[{li}]},


Do[
result=Join[result,Select[{li}[[n]],!IsContainedClqComp[result,#]&]];
,{n,2,length}];

result
];

(* Overload Pattern *)
ChooseNonIsomorphicClqComplexes[li_List/;(Depth[li]==4&&Length[{li}]==1)]:=
(*(*****************************)
(* Last Updated: 09/13/2023  *)
(*****************************)*)

(*given a single list of clique complexes (each given as a list of maximal cliques),
 it generates a new list composed of non-isomorphic clique complexes.*)
Block[{result={},remaining=li},


result=Reap[
While[Length[remaining]>0,
Sow[First[remaining]];
remaining=Select[remaining,IsomorphicGraphQ[GraphFromCliques[First[remaining]],GraphFromCliques[#]]==False&];
];
][[2,1]];

result
];

(* Overload Pattern *)
ChooseNonIsomorphicClqComplexes[li_List/;(Depth[{li}]==6&&Length[li]>1)]:=
(*(*****************************)
(* Last Updated: 09/30/2025  *)
(*****************************)*)

(*Given a list of lists of graph maximal cliques, it merges them into one list 
where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)Block[{result=li[[1]],length=Length[li]},


Catch[If[length==1,Throw[result, "ECGravReturn$41"]];

Do[
result=Join[result,Select[li[[n]],!IsContainedClqComp[result,#]&]];
,{n,2,length}];

result, "ECGravReturn$41"]
];

(* Catch-all Pattern *)
ChooseNonIsomorphicClqComplexes[args___]:=(Message[ChooseNonIsomorphicClqComplexes::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Simplicial Complex Isomorphism Classes*)


(* ::Item::Closed:: *)
(*IsomorphicSimplicialComplexQ*)


(* :Code Section: *)

(* Primary Pattern *)
IsomorphicSimplicialComplexQ[c1_List,c2_List]:=
(*
(****************************)
(*Last updated: 03/05/2024  *)
(****************************)
(*Checks whether two pure simplicial complexes given as lists of facets are isomorphic or not(i.e. whether or not there is a bijection of the vertex sets that is also a bijection of the facets.)
It does so by brute force enumeration of all isomorphism between their respective underlying graphs and testing if any of the isomorphisms are bijections of the facets*)
*)
With[{isomorphisms=Normal[FindGraphIsomorphism[GraphFromCliques[c1],GraphFromCliques[c2],All]]},
Catch[If[c1=={}&&c2=={},Throw[True, "ECGravReturn$42"]];
If[Length[isomorphisms]==0,Throw[False, "ECGravReturn$42"]];

If[Length[Select[isomorphisms, Complement[Sort/@(c1/.#),Sort/@c2]=={}&
]]>0,Throw[True, "ECGravReturn$42"]];

False, "ECGravReturn$42"]
];

(* Catch-all Pattern *)
IsomorphicSimplicialComplexQ[args___]:=(Message[IsomorphicSimplicialComplexQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*IsContainedSimplicialComp*)


(* Primary Pattern *)
IsContainedSimplicialComp[cmpxesLst_List/;Depth[cmpxesLst]==4,cmpx_List/;(Depth[cmpx]==3||cmpx=={})]:=
(*
(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)
(*Given a list of pure simplicial complexes cmpxesLst (each given as a list of facets of equal size)  and a single other pure simplicial complex cmpx, it returns True if there is a complex isomorphic to cmpx in cmpxesLst, False otherwise. It stops checking at the first occurence of isomorphic complex.*)
If[Length[cmpx]==0||cmpx=={{}},True,AnyTrue[cmpxesLst,IsomorphicSimplicialComplexQ[cmpx,#]&]
];

(* Catch-all Pattern *)
IsContainedSimplicialComp[args___]:=(Message[IsContainedSimplicialComp::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicSimplicialComplexes*)


(* Primary Pattern *)
SetAttributes[ChooseNonIsomorphicSimplicialComplexes, Flat];
ChooseNonIsomorphicSimplicialComplexes[li__List/;(Depth[{li}]==5&&Length[{li}]>1)]:=
(*
(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)
(*Given an arbitrary number of lists of facets of simplicial complexes, it merges them into one list of graphs where no two simplexes belonging to different input lists are isomorphic. This method does not check if the complexes WITHIN a given list are non-isomorphic.*)Block[{result={li}[[1]],length=Length[{li}]},


Do[
result=Join[result,Select[{li}[[n]],!IsContainedSimplicialComp[result,#]&]];
,{n,2,length}];

result
];

(* Overload Pattern *)
ChooseNonIsomorphicSimplicialComplexes[li_List/;(Depth[li]==4&&Length[{li}]==1)]:=(*(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)

(*given a single list of simplicial complex facet lists, it generates a new list composed of non-isomorphic facet lists.*)Block[{result={},remaining=li},

result=Reap[
While[Length[remaining]>0,
Sow[First[remaining]];
remaining=Select[remaining[[2;;-1]],IsomorphicSimplicialComplexQ[First[remaining],#]==False&];
];
][[2,1]];

result
];

(* Overload Pattern *)
ChooseNonIsomorphicSimplicialComplexes[li_List/;(Depth[li]==5)]:=
(*(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)

(*Given a list of lists of pure simplicial complexes, it merges them into one list where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)
Block[{result=li[[1]],length=Length[li]},


Catch[If[length==1,Throw[result, "ECGravReturn$43"]];

Do[
result=Join[result,Select[li[[n]],!IsContainedSimplicialComp[result,#]&]];
,{n,2,length}];

result, "ECGravReturn$43"]
];

(* Catch-all Pattern *)
ChooseNonIsomorphicSimplicialComplexes[args___]:=(Message[ChooseNonIsomorphicSimplicialComplexes::argerr, args];
$Failed);


(* ::Section::Closed:: *)
(*Automorphism Groups and Orders*)


(* ::Subsection::Closed:: *)
(*Simplicial Complex Automorphism and Facet Automorphism*)


(* ::Item::Closed:: *)
(*SimplicialComplexAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
SimplicialComplexAutomorphismGroupOrderConn[facetsLst_List]:=

(*
(****************************************)
(*   (* Last updated 2/22/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  Much more improved efficiency by first reducing the leaves of the graph,  *) *)
(****************************************)
(* computes the automorphism group order of a connected simplicial complex (not necessarily a clique complex and not necessarily pure) whose facets are given by facetsLst. It assumes the complex is connected. It returns Null if it is not connected.*)
(*Note: It only works by first relabeling the complex into canonically labeled complex*)
*)

Module[{purity,g,VertexFacetDegree,leafsLst,relabelingRule,reducedFacetsParentsAsn,relabeledreducedFacetsParentsAsn,reducedGraph, reducedAutG,canonicalReducedFacets, thereAreMissingCanonicalReducedFacets,isPermutationOfFacets,relabeledfacetsLst},

Catch[If[Length[facetsLst]==1,Throw[Length[facetsLst[[1]]]!, "ECGravReturn$44"]];


g=GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Throw[Null, "ECGravReturn$44"]
];

VertexFacetDegree[v_Integer]:=Length[Select[facetsLst,MemberQ[#,v]&]];


leafsLst=Table[Select[i,VertexFacetDegree[#]==1&],
{i,facetsLst}];
 

reducedFacetsParentsAsn=<|Table[With[{leafs=leafsLst[[i]]},
If[Length[leafs]>=1,Sort[Join[{leafs[[1]]},Complement[facetsLst[[i]],leafs]]],facetsLst[[i]]]->facetsLst[[i]]],{i,1,Length[facetsLst]}]|>;


reducedGraph=CanonicalGraph[GraphFromCliques[Keys[reducedFacetsParentsAsn]]];


relabelingRule=Normal[FindGraphIsomorphism[GraphFromCliques[Keys[reducedFacetsParentsAsn]],reducedGraph][[1]]];


relabeledreducedFacetsParentsAsn=<|Table[Sort[(i/.relabelingRule)]->reducedFacetsParentsAsn[[Key[i]]],{i,Keys[reducedFacetsParentsAsn]}]|>;


canonicalReducedFacets=Keys[relabeledreducedFacetsParentsAsn];


reducedAutG=GraphAutomorphismGroup[reducedGraph];


thereAreMissingCanonicalReducedFacets=Complement[FindClique[reducedGraph,\[Infinity],All],canonicalReducedFacets,SameTest->(Complement[#1,#2]=={}&)];


isPermutationOfFacets[x_]:=AllTrue[canonicalReducedFacets,With[{permval=Sort[PermutationReplace[#,x]]},
((MemberQ[canonicalReducedFacets,permval])&&(Length[Lookup[relabeledreducedFacetsParentsAsn,Key[#],{}]]==Length[Lookup[relabeledreducedFacetsParentsAsn,Key[permval],{}]]))]&];


If[thereAreMissingCanonicalReducedFacets!={},reducedAutG=PermutationGroup[Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]]
];


GroupOrder[reducedAutG]*Product[Length[i]!,{i,leafsLst}], "ECGravReturn$44"]

];

(* Catch-all Pattern *)
SimplicialComplexAutomorphismGroupOrderConn[args___]:=(Message[SimplicialComplexAutomorphismGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*SimplicialComplexAutomorphismGroupOrder*)


(* Primary Pattern *)
SimplicialComplexAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 2/22/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Given a general simplicial complex (not necessarily clique complex and not necessarily pure) that may or may not be connected, it computes the order of the automorphism group. It uses SimplicialComplexAutomorphismGroupOrderConn to compute the automorphism group orders of each component and then combines them. *)
(*Note: It only works by first relabeling the complex into canonically labeled complex*)
*)
With[{components=Tally[ConnectedComplexComponents[facetsLst],IsomorphicSimplicialComplexQ[#1,#2]&]},


Product[(SimplicialComplexAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
SimplicialComplexAutomorphismGroupOrder[args___]:=(Message[SimplicialComplexAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroup*)


(* :Code Section *)

(* Primary Pattern *)
PureComplexAutomorphismGroup[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 1/19/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Computes the automorphism group of a pure simplicial complex (not necessarily a clique complex) given as a list of its facets. It does so in a brute force and inefficient way without breaking the complex up into its connected components. *)
(*Note: It only works by first relabeling the complex into canonically labeled complex, so the output group will have different labeleing convention than the input unless the input is canonically labeled. The output is a list with two elements, the first is the automorphism group, the second is the reverse relabeling rule.*)
*)
Module[{facets,g,relabelingRule,autGroup, missingfacets,doesntMixupMissingAndNonMissing},

Catch[g=CanonicalGraph[GraphFromCliques[facetsLst]];


(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[GraphFromCliques[facetsLst],g][[1]]];


If[Length[facetsLst]==1,Throw[{SymmetricGroup[Length[facetsLst[[1]]]],Reverse/@relabelingRule}, "ECGravReturn$45"]];

facets=Sort/@(facetsLst/.relabelingRule);


autGroup=GraphAutomorphismGroup[g];


doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingfacets,MemberQ[facets,Sort[PermutationReplace[#,x]]]&];

missingfacets=Complement[Join@@(Subsets[#,{Length[facets[[1]]]}]&/@FindClique[g,\[Infinity],All]),facets];


If[missingfacets!={},autGroup=PermutationGroup[Select[GroupElements[autGroup],doesntMixupMissingAndNonMissing[#]&]
];
];


{autGroup,(Reverse/@relabelingRule)}, "ECGravReturn$45"]

];

(* Catch-all Pattern *)
PureComplexAutomorphismGroup[args___]:=(Message[PureComplexAutomorphismGroup::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
PureComplexAutomorphismGroupOrderConn[facetsLst_List]:=

(*
(****************************************)
(*   (* Last updated 2/21/2024. *) *)
(*   (* Version: 2.0 *) *) 
(*   (* Note:  Much more improved efficiency by first reducing the leaves of the graph, computing the *) *)
(****************************************)
(* computes the automorphism group order of a connected pure simplicial complex 
  (not necessarily a clique complex) whose facets are given by facetsLst. It assumes 
   the complex is connected and pure. It returns Null if it is not connected or if it 
   is not pure.*)
(* Note: It only works by first relabeling the complex into canonically labeled complex, so the output group will have different labeleing convention than the input unless the input is canonically labeled*)
*)

Module[{purity,g,VertexFacetDegree,leafsLst,reducedFacets,reducedGraph, reducedAutG,canonicalReducedFacets, thereAreMissingCanonicalReducedFacets,isPermutationOfFacets,reducedAutGorder},


Catch[If[Length[facetsLst]==1,Throw[Length[facetsLst[[1]]]!, "ECGravReturn$46"]];

If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Throw[Null, "ECGravReturn$46"]];

g=GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Throw[Null, "ECGravReturn$46"]
];

VertexFacetDegree[v_Integer]:=Length[Select[facetsLst,MemberQ[#,v]&]];


leafsLst=Table[Select[i,VertexFacetDegree[#]==1&],
{i,facetsLst}];
 

reducedFacets=Table[With[{leafs=Select[i,VertexFacetDegree[#]==1&]},
If[Length[leafs]>=1,Sort[Join[{leafs[[1]]},Complement[i,leafs]]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];


thereAreMissingCanonicalReducedFacets=Complement[FindClique[reducedGraph,\[Infinity],All],canonicalReducedFacets,SameTest->(Complement[#1,#2]=={}&)];


isPermutationOfFacets[x_]:=AllTrue[canonicalReducedFacets,MemberQ[canonicalReducedFacets,Sort[PermutationReplace[#,x]]]&];


If[thereAreMissingCanonicalReducedFacets!={},reducedAutGorder=Length[Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]],reducedAutGorder=GroupOrder[reducedAutG]
];


reducedAutGorder*Product[Length[i]!,{i,leafsLst}], "ECGravReturn$46"]

];

(* Catch-all Pattern *)
PureComplexAutomorphismGroupOrderConn[args___]:=(Message[PureComplexAutomorphismGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroupOrder*)


(* Primary Pattern *)
PureComplexAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 1/19/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Given a pure simplicial complex (not necessarily clique complex) that may or may not be connected, it computes the automorphism group order. It uses PureComplexAutomorphismGroupOrderConn to compute the automorphism group orders of each component and then combines them. *)
(*Note: It only works by first relabeling the complex into canonically labeled complex, so the output group will have different labeleing convention than the input unless the input is canonically labeled*)
*)
With[{components=Tally[ConnectedComplexComponents[facetsLst],IsomorphicSimplicialComplexQ[#1,#2]&]},


Product[(PureComplexAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
PureComplexAutomorphismGroupOrder[args___]:=(Message[PureComplexAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetStabilizerGroupOrder*)


(* Primary Pattern *)
PureComplexFacetStabilizerGroupOrder[facetsLst_List]:=
(*Given a simplicial complex as a list of facets, it computes the order of the 
	facet stabilizer group. It does not check whether the complex is a clique complex or not, 
	it assumes it is. *)
(*
(****************************************)
(*   (* Last updated 05/01/2026. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  Computes the facet stabilizer group order by first turning the complex
				into facet-representation where each vertex is identified only by the 
				label of the facets it participates in. Then, the facet stabilizer
				group size is simply the product of the factorials of the degeneracies of the 
				facet-labeled vertices. Example, if the original complex is 
				{{1,2,3},{2,3,4},{4,5,6}}, its facet-labeled vertices are 
				{{1},{1,2},{1,2},{2,3},{3},{3}}; i.e,
				1->{1},2->{1,2},3->{1,2},4->{2,3},5->{3},6->{3}. 
				So multiplicities of facet-labeled vertices are 
				{1}->1, {1,2}->2, {2,3}->1, {3}->2. the facet stabilizer group permutes these 
				vertices amongst themselves only, so has size 1!*2!*1!*2! = 4.*) *)
(****************************************)
*)
With[{
	multiplicitiesOfFacetLabeledVertices=
		Counts[(First/@Position[facetsLst,#])&/@(DeleteDuplicates[Flatten[facetsLst]])]},

	Product[i!,{i,multiplicitiesOfFacetLabeledVertices}]
];


(* Catch-all Pattern *)
PureComplexFacetStabilizerGroupOrder[args___]:=(Message[PureComplexFacetStabilizerGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
PureComplexFacetAutomorphismGroupOrderConn[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 2.0 *) *) 
(*   (* Note:  *) *)
(****************************************)(*computes the facet-automorphism group order of a simpli clique complex whose maximal cliques are given by facetsLst. It does not check whether the complex is a clique complex or not, it assumes it is. It returns Null if it is not connected.*)
*)
Module[{purity,g,reducedFacets,reducedGraph, reducedAutG,canonicalReducedFacets, missingCanonicalReducedFacets,doesntMixupMissingAndNonMissing,reducedFacetStabilizerGroup},

Catch[If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Throw[Null, "ECGravReturn$47"]];

If[Length[facetsLst]==1,Throw[1, "ECGravReturn$47"]];

g=GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Throw[Null, "ECGravReturn$47"]
];

purity=Length[facetsLst[[1]]];

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];


missingCanonicalReducedFacets=Complement[Join@@(Subsets[#,{purity}]&/@FindClique[reducedGraph,\[Infinity],All]),canonicalReducedFacets];


doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingCanonicalReducedFacets,MemberQ[canonicalReducedFacets,Sort[PermutationReplace[#,x]]]&];


If[missingCanonicalReducedFacets!={},reducedAutG=PermutationGroup[Select[GroupElements[reducedAutG],doesntMixupMissingAndNonMissing[#]&]];
];


reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];


GroupOrder[reducedAutG]/GroupOrder[reducedFacetStabilizerGroup], "ECGravReturn$47"]

];

(* Catch-all Pattern *)
PureComplexFacetAutomorphismGroupOrderConn[args___]:=(Message[PureComplexFacetAutomorphismGroupOrderConn::argerr, args];$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroupOrder*)


(* Primary Pattern *)
PureComplexFacetAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  *) *)
(****************************************)
(*Given a pure simplicial complex (not necessarily a clique complex) as a list of facets, it computes the order of the facet automorphism group. It does not check whether the complex is pure or not, it assumes it is. *)
*)
With[{components=Tally[ConnectedComplexComponents[facetsLst],IsomorphicSimplicialComplexQ[#1,#2]&]},


Product[(PureComplexFacetAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
PureComplexFacetAutomorphismGroupOrder[args___]:=(Message[PureComplexFacetAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroup*)


(* :Code Section *)

(* Primary Pattern *)
PureComplexFacetAutomorphismGroup[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/09/2024. *) *)
(*   (* Version: 2.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Returns the facet automorphism group of a pure simplicial complex given as a list 
of facets *)
*)
Module[{facets,g,relabelingRule,autGroup, missingfacets,facetStabilizerGroup,doesntMixupMissingAndNonMissing},

Catch[If[Length[facetsLst]==1,Throw[SymmetricGroup[1], "ECGravReturn$48"]];

g=CanonicalGraph[GraphFromCliques[facetsLst]];


(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[GraphFromCliques[facetsLst],g][[1]]];


facets=Sort/@(facetsLst/.relabelingRule);


autGroup=GraphAutomorphismGroup[g];


doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingfacets,MemberQ[facets,Sort[PermutationReplace[#,x]]]&];

missingfacets=Complement[Join@@(Subsets[#,{Length[facets[[1]]]}]&/@FindClique[g,\[Infinity],All]),facets];


If[missingfacets!={},
	autGroup=PermutationGroup[Select[GroupElements[autGroup],doesntMixupMissingAndNonMissing[#]&]];
];


facetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[autGroup,#]]&/@facets)];


(*If[GroupOrder[facetStabilizerGroup]==1,Return[autGroup/.(Reverse/@relabelingRule)]];*)

If[GroupOrder[facetStabilizerGroup]==1,Throw[autGroup, "ECGravReturn$48"]];

(PermutationGroup[DeleteDuplicates[RightCosetRepresentative[facetStabilizerGroup,#]&/@GroupElements[autGroup]]])/.(Reverse/@relabelingRule), "ECGravReturn$48"]

];

(* Catch-all Pattern *)
PureComplexFacetAutomorphismGroup[args___]:=(Message[PureComplexFacetAutomorphismGroup::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Clique Complex Automorphism and Facet Automorphism*)


(* ::Item::Closed:: *)
(*CliqueFacetStabilixerGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
CliqueFacetStabilizerGroupOrderConn[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/09/2024. *) *)
(*   (* Version: 2.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*computes the facet stabilizer group order of a connected clique complex whose 
maximal cliques are given by facetsLst. It does not check whether the complex is a 
clique complex or not, it assumes it is. It returns Null if it is not connected.*)
*)
Module[{purity,g,leafsLst,reducedFacets,reducedGraph, reducedAutG,canonicalReducedFacets, reducedFacetStabilizerGroup},


Catch[If[Length[facetsLst]==1,Throw[Length[facetsLst[[1]]]!, "ECGravReturn$49"]];

g=GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Throw[Null, "ECGravReturn$49"]
];

purity=Length[facetsLst[[1]]];

leafsLst=Table[Select[i,VertexDegree[g,#]==purity-1&],
{i,facetsLst}];
 

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];


GroupOrder[reducedFacetStabilizerGroup]*Product[Length[i]!,{i,leafsLst}], "ECGravReturn$49"]

];

(* Catch-all Pattern *)
CliqueFacetStabilizerGroupOrderConn[args___]:=(Message[CliqueFacetStabilizerGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetStabilixerGroupOrder*)


(* Primary Pattern *)
CliqueFacetStabilizerGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  *) *)
(****************************************)
(*Given a clique complex as a list of facets, it computes the order of the facet stabilizer group. It does not check whether the complex is a clique complex or not, it assumes it is. *)
*)With[{components=Tally[ConnectedComplexComponents[facetsLst],IsomorphicGraphQ[GraphFromCliques[#1],GraphFromCliques[#2]]&]},


Product[(CliqueFacetStabilizerGroupOrderConn[i[[1]]]^(i[[2]])),{i,components}]

];

(* Catch-all Pattern *)
CliqueFacetStabilizerGroupOrder[args___]:=(Message[CliqueFacetStabilizerGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
CliqueFacetAutomorphismGroupOrderConn[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/09/2024. *) *)
(*   (* Version: 2.0 *) *) 
(*   (* Note:  *) *)
(****************************************)(*computes the facet-automorphism group order of a connected clique complex whose maximal cliques are given by facetsLst. It does not check whether the complex is a clique complex or not, it assumes it is. It returns Null if it is not connected.*)
*)
Module[{purity,g,reducedFacets,reducedGraph, reducedAutG,canonicalReducedFacets, reducedFacetStabilizerGroup},


Catch[If[Length[facetsLst]==1,Throw[1, "ECGravReturn$50"]];

g=GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Throw[Null, "ECGravReturn$50"]
];

purity=Length[facetsLst[[1]]];

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];


GroupOrder[reducedAutG]/GroupOrder[reducedFacetStabilizerGroup], "ECGravReturn$50"]

];

(* Catch-all Pattern *)
CliqueFacetAutomorphismGroupOrderConn[args___]:=(Message[CliqueFacetAutomorphismGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroupOrder*)


(* Primary Pattern *)

CliqueFacetAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  *) *)
(****************************************)
(*Given a clique complex as a list of facets, it computes the order of the facet automorphism group. It does not check whether the complex is a clique complex or not, it assumes it is. *)
*)With[{components=Tally[ConnectedComplexComponents[facetsLst],IsomorphicGraphQ[GraphFromCliques[#1],GraphFromCliques[#2]]&]},


Product[(CliqueFacetAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
CliqueFacetAutomorphismGroupOrder[args___]:=(Message[CliqueFacetAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroup*)


(* :Code Section *)

(* Primary Pattern *)
CliqueFacetAutomorphismGroup[facetsLst_List]:=(*
(****************************************)
(*   (* Last updated 09/30/2025. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Given a clique complex as a list of facets, it computes the facet automorphism group. 
It does not check whether the complex is a clique complex or not, it assumes it is. *)
*)
Module[{facets,g,relabelingRule,autGroup,facetStabilizerGroup},

Catch[If[Length[facetsLst]==1,Throw[SymmetricGroup[1], "ECGravReturn$51"]];

g=CanonicalGraph[GraphFromCliques[facetsLst]];


(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[GraphFromCliques[facetsLst],g][[1]]];


facets=Sort/@(facetsLst/.relabelingRule);


autGroup=GraphAutomorphismGroup[g];


facetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[autGroup,#]]&/@facets)];


If[GroupOrder[facetStabilizerGroup]==1,Throw[autGroup, "ECGravReturn$51"]];


(PermutationGroup[DeleteDuplicates[RightCosetRepresentative[facetStabilizerGroup,#]&/@GroupElements[autGroup]]])/.(Reverse/@relabelingRule), "ECGravReturn$51"]

];

(* Catch-all Pattern *)
CliqueFacetAutomorphismGroup[args___]:=(Message[CliqueFacetAutomorphismGroup::argerr, args];
$Failed);


(* ::Chapter:: *)
(*Generate Pure Simplicial Complexes*)


(* ::Section::Closed:: *)
(*Generate All Pure Simplicial Complexes*)


(* ::Subsection::Closed:: *)
(*Generate all vertex labeled pure simplicial complexes*)


(* ::Item::Closed:: *)
(*GenerateAllVertexLabeledPureSimplicialComplexes[{p,q,n}] or [{p,q}]*)


(* :Code Section: *)

(* Primary Pattern *)
GenerateAllVertexLabeledPureSimplicialComplexes[{purity_Integer,facetOrder_Integer, nV_Integer}]:=
(*
(****************************************)
(*   (* Last updated 08/26/2023. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(* Generates all vertex labeled pure complexes of a given purity = purity, facet order = 
facetOrder, and number of vertices nV. E.g. GenerateAllPureSimplicialComplexes[2,4,4] 
returns a list of all (connected and disconnected) 
2-pure vertex-labeled simplicial complexes with 4 edges and 4 vertices. *)
*)

Module[{po=purity*facetOrder,facets=Subsets[Range[nV]-1,{purity}],zeroSeqLength,
partitions,capacityCheckOk,allowedPartitionPermutations,addOneClique,result={}},


zeroSeqLength[ind_Integer,cmpsn_List]:=(*counts the number of continuous sequences of 
zeros following the position ind in the composition of n given as a list cmpsn *)
Module[{count=0},
	Do[
		If[cmpsn[[j]]==0,count++,Break[]]
	,{j,ind+1,Length[cmpsn]}];
count
];

capacityCheckOk[lst_List,p_Integer,q_Integer]:=
(*For a given composition of n as a sum of q numbers all of which are less than or 
equal to the purity, this subroutine checks the maximum number of sequences of 
continuous zeros that can be accomodated based on the number of vertices available so 
far*)
Module[{capacity,res},
capacity[ind_Integer]:=Binomial[Total[lst[[1;;ind]]],p]-ind;

res=Catch[Do[If[capacity[k]<zeroSeqLength[k],Throw[False]],{k,1,Length[lst]}]];

If[res==Null,res=True];

res
];

partitions=PadRight[#,facetOrder]&/@(Select[IntegerPartitions[nV,facetOrder,Range[purity]],MemberQ[#,purity]&]);


allowedPartitionPermutations=Join@@Table[Select[Permutations[i],(#[[1]]==purity&&capacityCheckOk[#,purity,facetOrder])&],{i,partitions}];


addOneClique[curlst_List,ptn_List,facetsLst_List]:=
(*
(* recursively adds one more facet to a given initial list of facets (curlst) based on 
a list of number of new vertices to be added (or vertex count increment order given by 
ptn). ptn is a partition of the number of vertices n of the simplicial complex into 
the sum of clique order values none of which is greater than the purity of the complex 
obeying some restriction). The facetsLst is the list of all passible facets which is 
the set of purity-subsets of {0,1,...,n-1}. 
E.g. addOneClique[{{1,2}},{2,1,0},{{1,2},{1,3},{2,3}}] outputs {{1,2},{1,3},{2,3}} 
where the number of new vertices at each step is 2 then 1 then 0. *)
*)
Module[{curCt=Length[curlst],curVerts=DeleteDuplicates[Flatten[curlst]],
	clqCt=Length[ptn],lastpos=If[curlst=={},0,Position[facetsLst,curlst[[-1]]][[1,1]]],cands},


Catch[

	If[curCt==clqCt,Throw[curlst],


		cands=Select[facetsLst[[lastpos+1;;-1]],
				Length[Complement[#,curVerts]]==ptn[[curCt+1]]&];
				(*all possible new facets that increase the number of vertices by the 
				value of the partition at the next position to curCt*)


		If[cands!={},
			addOneClique[Append[curlst,#],ptn,facetsLst]&/@cands(*recursive call*)
		,Nothing]

	]
]

];

result=Table[addOneClique[{},j,facets],{j,allowedPartitionPermutations}];


result=Apply[Join,result,{0,facetOrder-1}];

result+1(* the +1 turns the vertex labels from {0,1,...n-1} to {1,2,...,n}*)
];

(* Overload Pattern *)
GenerateAllVertexLabeledPureSimplicialComplexes[{purity_Integer,facetOrder_Integer}]:=
With[{nmin=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]]},
	Join@@ParallelTable[
		GenerateAllVertexLabeledPureSimplicialComplexes[{purity,facetOrder,n}],
	{n,nmin,purity*facetOrder},DistributedContexts->{$Context,"ECGrav`Private`"}]
]


(* Catch-all Pattern *)
GenerateAllVertexLabeledPureSimplicialComplexes[args___]:=(Message[GenerateAllVertexLabeledPureSimplicialComplexes::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Generate all facet-labeled pure simplicial complexes*)


(* ::Item::Closed:: *)
(*GenerateAllFacetLabeledPureSimplicialComplexes*)


(* :Code Section: *)

(* Primary Pattern *)
GenerateAllFacetLabeledPureSimplicialComplexes[{purity_Integer,facetOrder_Integer}]:=
(*
(****************************************)
(*   (* Last updated 10/11/2024. *) *)
(*   (* Version : 1 .*) *) 
(*   (* Note:  .*) *) 
(****************************************)
(*A code to generate all  pure facet labeled(possibly disconnected) simplicial complex 
of a given purity and facet order through iterative addition of vertices labeled by 
their facets. Input is two integers, p (purity) and q (facet order).
E.g. ECGrav`GenerateAllFacetLabeledPureSimplicialComplexes[2,4] generates all 2-pure 
facet-labeled simplicial complex with 4 facets.
It successively adds a vertex from a decreasing list of available vertex space.  *)

*)
Module[{complexLst,facetSubsets=Subsets[Range[facetOrder]][[2;;-1]],AddOneVertex,BuildComplex,GetFacetsFromVList,result},


AddOneVertex[curLstOfVertices_List]:=
	Module[{availableVertices,unAvailableVertices={},vertexDegree,newVertex={}},


	vertexDegree=<|Table[i->Length[Select[curLstOfVertices,SubsetQ[#,i]&]],{i,facetSubsets}]|>;


	unAvailableVertices=Table[If[(Length[i]==1&&vertexDegree[[Key[i]]]==purity)||(Length[i]>1&&vertexDegree[[Key[i]]]==purity-1),i,Nothing],{i,facetSubsets}];


	availableVertices=Fold[Function[{x,y},Select[x,!SubsetQ[#,y]&]],facetSubsets,unAvailableVertices];


	availableVertices=Select[availableVertices,Order[#,curLstOfVertices[[-1]]]>=0&];


	If[availableVertices=={},Nothing,Join[curLstOfVertices,{#}]&/@availableVertices]


];


BuildComplex[]:=
	Module[{finishedSXLst={},curSXLst={}, nextSxLst={#}&/@facetSubsets[[facetOrder;;-1]],iterCount=0},


	While[nextSxLst!={}&&iterCount<=purity*facetOrder,
		iterCount++;
		curSXLst=nextSxLst;

		nextSxLst=Join@@(
			ParallelMap[AddOneVertex[#]&,curSXLst,DistributedContexts->{$Context,"ECGrav`Private`"}]);
		finishedSXLst=Join[finishedSXLst,Select[nextSxLst,Length[Join@@#]==purity*facetOrder&]];
		nextSxLst=Complement[nextSxLst,finishedSXLst];

		(*Print[""];
		Print[""];
		Print["     In BuildComplex while loop, after adding one vertex, curSX ",curSXLst, " nextSx ",nextSxLst," finishedSXLst ",finishedSXLst];*)

	];

	finishedSXLst

];


complexLst=BuildComplex[];


GetFacetsFromVList[cmpxVlst_List]:=With[{cmlxAsn=<|Table[i->cmpxVlst[[i]],{i,1,Length[cmpxVlst]}]|>},
	Table[
		With[{k=Select[cmpxVlst,MemberQ[#,q]&]},
			Keys[Select[cmlxAsn,MemberQ[k,#]&]]
		]
	,{q,1,facetOrder}]
];

result=Reverse[GetFacetsFromVList[#]]&/@complexLst;


result

];

(* Catch-all Pattern *)
GenerateAllFacetLabeledPureSimplicialComplexes[args___]:=(Message[GenerateAllFacetLabeledPureSimplicialComplexes::argerr, args];
$Failed);


(* ::Section:: *)
(*Generate A Random Pure Simplicial Complex*)


(* ::Subsection::Closed:: *)
(*Random [purity*facetOrder]-labeled Pure Simplicial Complex labeled by *)


(* ::Item::Closed:: *)
(*RandomPQLabeledPureSimplicialComplex*)


(* :Code Section: *)

(* Primary Pattern *)
RandomPQLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
(*
(****************************************)
(*   (* Last updated 08/01/2024. *) *)
(*   (* Virsion: 2 *) *)
(*   (* Note: The first facet is free to be anything compared to version 1. this 
			version gives the correct distribution for pqLabeled complexes. *) *) 
(****************************************)
(*A code to generate a random pure (possibly disconnected) simplicial complex of a 
	given purity and clique order. Labeled using the range of integers 1 through 
	purity*facetOrder. Input is a list of two numbers, {a, b} where a is the purity 
	number and b is the number of cliques. E.g. RandomPureSimplicialComplexSimp[{2,4}] 
	generates a 2-pure random simplicial complex with 4 edges. It successively choses 
	a random facet from a growing list of available set space.  *)
*)
Module[{facets={},labelsToChoseFrom=Range[purity*facetOrder],newfacet={},iterNum},


Do[
newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

(*Print["newfacet ",newfacet];
Print["memberQ ",MemberQ[facets,newfacet]];*)

iterNum=0;

While[MemberQ[facets,newfacet]&&iterNum<100,
iterNum++;
newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

];


facets=Join[facets,{newfacet}],
{i,1,facetOrder}];

Sort[facets]
];

(* Overload Pattern *)

RandomPQLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer},numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 10/03/2025. *) *)
(*   (* Virsion: 2 *) *)
(*   (* Note: This version gives the correct distribution for pqLabeled complexes. *) *) 
(****************************************)
(*A code to generate numSamples random pure (possibly disconnected) simplicial complexes of a 
	given purity and clique order. Labeled using the range of integers 1 through 
	purity*facetOrder. Input is a list of two numbers, {a, b} where a is the purity 
	number and b is the number of cliques. E.g. RandomPureSimplicialComplexSimp[{2,4}] 
	generates a 2-pure random simplicial complex with 4 edges. It successively choses 
	a random facet from a growing list of available set space.  *)
*)
Module[{facets={},labelsToChoseFrom=Range[purity*facetOrder],newfacet={},iterNum, result},


If[numSamples<10^5, 
	result = 
	Table[
		Do[
			newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

			(*Print["newfacet ",newfacet];
			Print["memberQ ",MemberQ[facets,newfacet]];*)

			iterNum=0;

			While[MemberQ[facets,newfacet]&&iterNum<100,
			iterNum++;
			newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

			];


			facets=Join[facets,{newfacet}],
			{i,1,facetOrder}
			];

		Sort[facets]
	
		,{numSamples}
	],

	result = 
	ParallelTable[
		Do[
			newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

			(*Print["newfacet ",newfacet];
			Print["memberQ ",MemberQ[facets,newfacet]];*)

			iterNum=0;

			While[MemberQ[facets,newfacet]&&iterNum<100,
				iterNum++;
				newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

			];


			facets=Join[facets,{newfacet}],
			{i,1,facetOrder}
			];

		Sort[facets]
	
	,{numSamples},DistributedContexts->{$Context,"ECGrav`Private`"}]
	];

result
];

(* Catch-all Pattern *)
RandomPQLabeledPureSimplicialComplex[args___]:=(Message[RandomPQLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random Vertex labeled Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomVertexLabeledPureSimplicialComplex*)


(* Private helpers for RandomVertexLabeledPureSimplicialComplex. -----------------------------
   These are deliberately top-level rather than Module-local closures inside the generator, and
   take purity/facetOrder/nV as arguments instead of capturing them. A fresh Module per call
   would mint new temporary symbols each time, and ParallelTable re-serialises those to every
   subkernel on every call: measured at 0.06 s for {2,3,4} up to 0.26 s for {4,20,40}, growing
   with the weight table and still 45% of the total at 20000 samples. Stable symbols distribute
   once per session instead, which drops the fixed cost to a flat 0.010 s and is what lets the
   parallel threshold below sit at 1000 rather than 10^4. *)

(* Private helper *)
RandomPureComplexCDFRow[purity_Integer,q_Integer,v_Integer]:=
(*Cumulative composition weights for the qth facet, given that the first q facets cover v
vertices. The weight of "the qth facet brings k new vertices" is
	NumVertexLabeledPureComplexes[p,q-1,v-k] (Binomial[v,v-k] Binomial[v-k,p-k] - [k==0](q-1)),
i.e. the sub-complex on the other q-1 facets, times the choice of which v-k vertices they
cover, times the choice of the qth facet itself; the [k==0] term drops the facets that would
repeat one already present. Accumulated so a composition can be drawn with one RandomInteger
and a scan, and kept in exact integers because the weights are counts running to hundreds of
digits.*)
	Accumulate@If[q>2,
		Table[NumVertexLabeledPureComplexes[purity,q-1,v-k]*(Binomial[v,v-k]*Binomial[v-k,purity-k]-If[k==0,q-1,0]),{k,0,purity}],
		Table[NumVertexLabeledPureComplexes[purity,q-1,v-k]*(Binomial[v,v-k]*Binomial[v-k,purity-k]),{k,1,purity}]];

(* Private helper: facet order above which composition weights are recomputed instead of kept.
   A row depends on (purity,q,v) alone -- never on the sample, and never on the vertex count of
   the call -- so it is worth keeping: rebuilding these per sample was about 85% of the old
   runtime, and they are now shared across calls as well, which is worth a further 1.1x to 1.3x
   on repeated calls. Capping bounds what that retains, since a cached q also bounds v (a row
   for facet index q is only ever asked for v <= p q). Measured cumulative size at this cap:
   45 KB after sampling {3,16,30}, 200 KB by {4,20,40}, 2.6 MB by {4,100,200}. Same shape as
   $NumPCMaxCachedQ, and cleared by the same NumPCClearCache[]. *)

$RandomComplexMaxCachedQ=150;

(* Private helper *)
RandomPureComplexCDF[purity_Integer,q_Integer,v_Integer]:=
	If[q<=$RandomComplexMaxCachedQ,
		RandomPureComplexCDF[purity,q,v]=RandomPureComplexCDFRow[purity,q,v],
		RandomPureComplexCDFRow[purity,q,v]];

(* Private helper *)
RandomPureComplexComposition[purity_Integer,facetOrder_Integer,nV_Integer]:=
(*How many new vertices each facet brings, drawn back to front. The qth entry is sampled from
the cumulative weights above by inverse transform: one RandomInteger in [1,total] and a scan
for the first cumulative weight reaching it. q>2 offers k=0..purity, q==2 offers k=1..purity,
and the first facet is all new by definition.*)
Module[{newkTable=ConstantArray[0,facetOrder],curVCount=nV,cdf,r},

	Do[
		cdf=RandomPureComplexCDF[purity,q,curVCount];
		r=RandomInteger[{1,Last[cdf]}];
		newkTable[[q]]=LengthWhile[cdf,#<r&]+If[q>2,0,1];
		curVCount-=newkTable[[q]];
	,{q,facetOrder,2,-1}];

	newkTable[[1]]=purity;

	newkTable
];

(* Private helper *)
RandomPureComplexFacets[purity_Integer,facetOrder_Integer,nV_Integer]:=
(*One complex. The composition needs no feasibility re-draw: the old generator looped while any
prefix of j facets covering v_j vertices had Binomial[v_j,purity]<j -- too few distinct
p-subsets to fill j facets -- but a k is only ever drawn with positive weight, and a positive
weight carries the factor NumVertexLabeledPureComplexes[purity,q-1,v-k], which vanishes unless a complex
with q-1 facets on those v-k vertices exists. That is the very condition being retested, so
with the caller's guard pinning the top of the chain the loop could never fire, and it cost an
O(facetOrder^2) rescan per sample.

Vertices come from one random permutation of [nV] instead of repeated sampling of the uncovered
set: the unused tail of a uniform permutation is a uniform random subset in uniform random
order, so taking its next k entries is the same draw as RandomSample of the uncovered vertices,
without the Complement bookkeeping.*)
Module[{composition=RandomPureComplexComposition[purity,facetOrder,nV],
		perm=RandomSample[Range[nV]],facets=ConstantArray[{},facetOrder],coveredCount=0,k,newfacet},

	Do[
		k=composition[[j]];
		newfacet=Sort[Join[RandomSample[Take[perm,coveredCount],purity-k],
			Take[perm,{coveredCount+1,coveredCount+k}]]];

		(*A facet bringing no new vertex may repeat an earlier one; one bringing a new vertex
		cannot. The k=0 weight above excludes the q-1 facets already placed, so it is only ever
		drawn when an unused p-subset remains and this loop terminates.*)
		If[k==0,
			While[MemberQ[facets,newfacet],
				newfacet=Sort[RandomSample[Take[perm,coveredCount],purity]]
			];
		];

		coveredCount+=k;
		facets[[j]]=newfacet;

	,{j,1,facetOrder}];

	facets
];

(* Private helper: sample count above which the draw is farmed out to ParallelTable. ---------
   Set from the measured crossover on 11 subkernels with kernels already up: with the helpers
   above hoisted, parallel overtakes serial at about 500 samples for {3,16,30}, 1000 for
   {5,10,30} and 2000 for {3,6,10}. (Before hoisting the same measurement gave 6000 to 15000,
   which is why the threshold used to be 10^4.) Tune with
	ECGrav`Private`$RandomComplexParallelThreshold = 5000;
   or set it to Infinity to keep the generator serial.

   The $KernelCount test is the other half of the decision. ParallelTable launches kernels on
   its own if none are running, which costs about 3.7 s -- more than forty times the whole
   serial draw at the threshold, and a surprising bill for a call that looks small. With kernels
   already up, distribution costs about 10 ms. So parallelise only when the session has already
   paid for its kernels, and otherwise stay serial rather than quietly launching 11 of them. The
   notebook workflow this package is used from calls LaunchKernels[] up front, as does
   Tests/TestPrelude.wl, so the parallel path is live where it is wanted. *)

$RandomComplexParallelThreshold=1000;

(* Private helper *)
RandomPureComplexSample[purity_Integer,facetOrder_Integer,nV_Integer,numSamples_Integer]:=
	If[numSamples>$RandomComplexParallelThreshold&&$KernelCount>0,
		ParallelTable[RandomPureComplexFacets[purity,facetOrder,nV],{numSamples},
			DistributedContexts->{$Context,"ECGrav`","ECGrav`Private`"}],
		Table[RandomPureComplexFacets[purity,facetOrder,nV],{numSamples}]];


(* Primary Pattern *)

RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 03/30/2026. *) *)
(*   (* Note: A totally new approach from 
			previous versions based onrecursive picking of the number of new vertices 
			in the last facet. Very efficient. Gives the correct distribution for 
			vertex labeled pure complexes. *) 
(****************************************)
(*A code to generate a random sample of (possibly disconnected)vertex labeled pure 
simplicial complexes of a given purity, facet order, and number of vertices 
through iterative addition of facets. Outputs a list numSamples of them. Input is a 
list with numbers, {p,q,n} and an integer, numSamples; 
where p is the purity number, q is the facet order, n is the number of vertices
and numSamples is the number of samples required. E.g. 
RandomVertexLabeledPureSimplicialComplex[{2,4,6},10] generates 10, 2-pure random 
simplicial complexes with 4 edges and 6 vertices.     *)
*)
*)
Catch[

If[numSamples<0,
	Message[RandomVertexLabeledPureSimplicialComplex::argerr,{purity,facetOrder,nV},numSamples];
	Throw[$Failed, "ECGravReturn$63"]];

(*Refuse parameters with nothing to sample. Every composition weight carries a factor
NumVertexLabeledPureComplexes[purity,q-1,...], so on an empty sample space they are all zero; the draw then
fails and its unevaluated result feeds the next stage, which used to return "complexes" with
RandomSample[...] expressions sitting inside them rather than failing.*)
If[purity<1||facetOrder<1||NumVertexLabeledPureComplexes[purity,facetOrder,nV]==0,
	Message[RandomVertexLabeledPureSimplicialComplex::empty,purity,facetOrder,nV];
	Throw[$Failed, "ECGravReturn$63"]];

RandomPureComplexSample[purity,facetOrder,nV,numSamples]

, "ECGravReturn$63"];


(* Overload Pattern *)

RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 03/30/2026. *) *)
(*   (* Note: A totally new approach from 
			previous versions based onrecursive picking of the number of new vertices 
			in the last facet. Very efficient. Gives the correct distribution for 
			vertex labeled pure complexes. *) 
(****************************************)
(*A code to generate a random sample of (possibly disconnected)vertex labeled pure 
simplicial complexes of a given purity and clique order through iterative addition of 
facets. Outputs a list numSamples of them. Input is a list of three numbers, p,q, 
numSamples, where p is the purity number, q is the facet order, and numSamples is 
the number of samples required. numSamples is set to 1 by default. E.g. 
RandomUnlLabeledPureSimplicialComplex[{2,4},10] generates 10, 2-pure random simplicial 
complex with 4 edges.     *)
(*Algorithm:, 
1. Pick the number of vertices n between n=nmin and n=p*q weighted by sv(p,q,n)., 
2. Once numVertices n is picked, it simply randomly chooses q p-combinations (the facets) 
	successively one at a time by randomly sampling [n]. If at the kth stage the random facet 
	picked has already been picked previously, it tries again., 
3. After generating q facets, if their union doesn't cover [n], it goes back to step 1 and repeats.*)
*)
*)
Module[{nmin,svTable},

Catch[

If[numSamples<0,
	Message[RandomVertexLabeledPureSimplicialComplex::argerr,{purity,facetOrder},numSamples];
	Throw[$Failed, "ECGravReturn$64"]];

(*See the three-argument pattern above: with no vertex count to pin down there is still nothing
to sample unless purity and facet order are at least 1.*)
If[purity<1||facetOrder<1,
	Message[RandomVertexLabeledPureSimplicialComplex::empty,purity,facetOrder,"any"];
	Throw[$Failed, "ECGravReturn$64"]];

(*Set nmin*)
nmin=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]];


(*table of sv(p,q,n)*)
svTable=Table[NumVertexLabeledPureComplexes[purity,facetOrder,n]*1.0,{n,nmin,purity*facetOrder}];

(*Normalize svTable by the geometric mean for easier computation*)
svTable=svTable/GeometricMean[svTable];

(*Each sample picks its own vertex count from sv(p,q,n) and then hands off to the same shared
builder as the three-argument pattern. The two tables are substituted in by With rather than
referenced as Module locals, so the parallel branch ships them as literal data instead of
having to distribute a pair of temporary symbols.*)
With[{weights=svTable,vertexCounts=Range[nmin,purity*facetOrder]},
	If[numSamples>$RandomComplexParallelThreshold&&$KernelCount>0,
		ParallelTable[
			RandomPureComplexFacets[purity,facetOrder,RandomChoice[weights->vertexCounts]],
			{numSamples},DistributedContexts->{$Context,"ECGrav`","ECGrav`Private`"}],
		Table[
			RandomPureComplexFacets[purity,facetOrder,RandomChoice[weights->vertexCounts]],
			{numSamples}]]]

, "ECGravReturn$64"]

];


(* Overload Pattern *)

(*First only if there is something to take: the one-sample forms above return $Failed on an
empty sample space, and First[$Failed] would bury that under a Part error.*)
RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
With[{sample=RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder}, 1]},
	If[ListQ[sample],First[sample],sample]];


(* Overload Pattern *)

RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}]:=
With[{sample=RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder, nV}, 1]},
	If[ListQ[sample],First[sample],sample]];


(* Catch-all Pattern *)
RandomVertexLabeledPureSimplicialComplex[args___]:=(Message[RandomVertexLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random Uniform Unlabeled Pure Simplicial Complex by Burnside pair sampling*)


(* Exact uniform sampling of fully unlabeled complexes. ---------------------------------------
   Neither the vertices nor the facets are labelled, so the objects are the isomorphism classes:
   the S_n orbits on the covering M-element sets of distinct p-subsets, which is what
   NumUnlabeledPureComplexes counts.

   Pair sampling, as for the facet-labeled generator. Draw uniformly from the fixed pairs
   {(sigma,S) : sigma in S_n, S a sigma-invariant covering M-set} and discard sigma. Every orbit O
   contributes Sum_{S in O} |Stab(S)| = |O| (n!/|O|) = n! pairs -- the same number whatever the
   orbit's size -- so the class comes out exactly uniform, with no automorphism weighting and no
   rejection. A sigma-invariant set is a union of <sigma>-orbits of p-subsets, giving five steps:

	1. a cycle type lambda, weighted by (permutations of that type) x (covering sets they fix);
	2. the canonical sigma_lambda, whose cycles are consecutive blocks;
	3. a size PROFILE a, where a_d is how many orbits of size d the set uses;
	4. the orbits themselves, one at a time, each weighted by its number of completions;
	5. a uniform relabelling of [n].

   Steps 2 and 5 stand in for "realise a uniform sigma of cycle type lambda": conjugating the
   canonical sigma_lambda by a uniform permutation is the same thing, and it makes every orbit
   computation a function of lambda alone, so it is memoised across samples rather than redone per
   sample.

   Step 3 is what makes step 4 tractable. Reading the count the other way round -- as orbits of
   separating incidence tableaux under S_M rather than S_n orbits of facet sets -- the acting
   permutation tau of the facet labels has a cycle type that IS the size profile, since a tau-cycle
   of length d pins one sigma-orbit of size d. So the profile is not hidden state to be inferred
   during the draw; it is sampled up front from its own marginal. Without it the completion count
   depends on the profile of what has been used so far, and step 4's state would have to carry it.

   Covering is imposed during the draw rather than by rejecting afterwards. Rejecting is simpler
   and exact, but its expected trial count is exactly A(p,M,n)/U(p,M,n) with A the non-covering
   count -- tolerable mid-range and hopeless at the top, where the facets must nearly partition
   [n]: 66 trials at {3,4,12}, 4279 at {3,6,18}, over 10^7 at {4,7,28}.

   Every count here is one inclusion-exclusion over DELETED cycles, because the sets avoiding a
   cycle set T are exactly the sets of the permutation with those cycles removed -- which the
   counter's own NumULPCFixSets already handles at the reduced cycle type. The same argument gives
   the completion counts: an orbit already chosen has its support inside the covered cycles, so it
   avoids every deleted set and survives into every term, and only the number used OF EACH SIZE
   enters. The total alone is not enough, which is where this differs from the facet-labeled
   sampler's lemma: there every chosen object cost one unit of the budget, so the profile was the
   total; here orbit sizes exceed one and two histories can spend the same budget in different
   denominations. *)

(* Private helper: sample count past which the parallel branch pays. ------------------------
   Measured, not guessed, as for the other generators -- and it lands far higher than theirs.
   Almost all the work here is in memo tables keyed on the cycle type, and every subkernel has
   to build its own, so the parallel branch pays a fixed cost the serial one pays once. On 11
   subkernels with kernels already up, serial/parallel runs 0.24-0.53 at 100 samples, 0.31-0.81
   at 200, 0.78-1.48 at 400, and 1.47-2.71 at 800; by 3200 it is 2.7-4.2. So the crossover sits
   between 400 and 800 and 500 is the conservative side of it -- against 100 for the
   facet-labeled generator, whose per-sample state is far cheaper to rebuild. Gated on
   $KernelCount as everywhere else, so a small call never quietly launches kernels. *)
$RandomUnlabeledParallelThreshold=500;

(* Private helper *)
RandULPCGoParallel[numSamples_Integer]:=
	numSamples>$RandomUnlabeledParallelThreshold&&$KernelCount>0;

(* Private helper *)
RandULPCOrbits[lambda_List,p_Integer]:=RandULPCOrbits[lambda,p]=
(*The <sigma>-orbits on the p-subsets of [n] for the canonical sigma of cycle type lambda, each as
{size, support, members}, the support being the set of cycle indices the subset meets. This is the
only place the C(n,p) subsets are touched, and being keyed on lambda alone it is shared by every
sample that draws this cycle type.*)
Module[{pi={},cyc=<||>,off=0,seen=<||>,out={},orb,cur},
	Do[pi=Join[pi,RotateLeft[Range[off+1,off+k]]];off+=k,{k,lambda}];
	off=0;
	Do[Scan[(cyc[#]=j)&,Range[off+1,off+lambda[[j]]]];off+=lambda[[j]],{j,Length[lambda]}];
	Do[
		If[!KeyExistsQ[seen,s],
			orb={};cur=s;
			While[!MemberQ[orb,cur],AppendTo[orb,cur];cur=Sort[pi[[cur]]]];
			Scan[(seen[#]=True)&,orb];
			AppendTo[out,{Length[orb],Union[cyc/@s],orb}]]
	,{s,Subsets[Range[Total[lambda]],{p}]}];
	out
];

(* Private helper *)
RandULPCNd[lambda_List,p_Integer]:=RandULPCNd[lambda,p]=
(*d -> the number of <sigma>-orbits of size d on the p-subsets, for cycle type lambda. Taken from
the counter's f by Moebius inversion rather than by enumerating orbits, so the reduced cycle types
the inclusion-exclusion asks for cost nothing extra. Callers pass lambda sorted, so that the
reduced types they build and the ones here agree as memo keys.*)
If[lambda==={},
	<||>,
	Module[{parts=Tally[lambda],divs,f},
		divs=Divisors[LCM@@lambda];
		f=Association[Table[e->NumULPCFixedSubsets[parts,e,p],{e,divs}]];
		Select[
			Association[Table[d->Total[Table[MoebiusMu[d/e]*f[e],{e,Divisors[d]}]]/d,{d,divs}]],
			#>0&]
	]
];

(* Private helper *)
RandULPCSubMulti[lambda_List]:=RandULPCSubMulti[lambda]=
(*Every sub-multiset of lambda with its inclusion-exclusion sign and the number of ways to choose
which cycles of each length survive. Deleting cycles is how covering is imposed throughout.*)
Module[{t=Tally[lambda],mult},
	mult=t[[All,2]];
	Table[
		With[{sub=DeleteCases[Transpose[{t[[All,1]],j}],{_,0}]},
			{Sort[Flatten[Table[ConstantArray[q[[1]],q[[2]]],{q,sub}]]],
			 (-1)^(Total[mult]-Total[j])*(Times@@MapThread[Binomial,{mult,j}])}]
	,{j,Tuples[Range[0,#]&/@mult]}]
];

(* Private helper *)
RandULPCFixCov[lambda_List,p_Integer,MM_Integer]:=RandULPCFixCov[lambda,p,MM]=
(*The number of sigma-invariant COVERING M-sets for sigma of cycle type lambda. The counter never
forms this -- it reaches covering by differencing over the vertex count, which says nothing per
cycle type -- so this is the piece the sampler adds.*)
Total[Function[sm,sm[[2]]*NumULPCFixSets[Tally[sm[[1]]],p,MM]]/@RandULPCSubMulti[lambda]];

(* Private helper *)
RandULPCTypeWeights[p_Integer,MM_Integer,n_Integer]:=RandULPCTypeWeights[p,MM,n]=
(*Step-1 weights: {cycle types, (permutations of that type) x (covering sets they fix)}. These sum
to n! NumUnlabeledPureComplexes[p,MM,n] exactly, which the suite checks: a wrong cycle-type weight
still emits plausible-looking complexes, so the identity is what to verify, not the output.*)
Module[{parts=IntegerPartitions[n]},
	{parts,
	 Table[(n!/(Times@@(#[[1]]^#[[2]]*#[[2]]!&/@Tally[lambda])))*RandULPCFixCov[lambda,p,MM],
		{lambda,parts}]}
];

(* Private helper *)
RandULPCProfiles[lambda_List,p_Integer,MM_Integer]:=RandULPCProfiles[lambda,p,MM]=
(*The admissible size profiles: how many orbits of each size, the sizes weighted by d summing to
MM, and no more of a size taken than exist.*)
Module[{nd=Counts[#[[1]]&/@RandULPCOrbits[lambda,p]],sizes,ranges},
	sizes=Sort[Keys[nd]];
	If[sizes==={},
		{},
		ranges=Table[Range[0,Min[nd[d],Quotient[MM,d]]],{d,sizes}];
		Select[Association[Thread[sizes->#]]&/@Tuples[ranges],
			Total[Times@@@Transpose[{sizes,Values[#]}]]===MM&]
	]
];

(* Private helper *)
RandULPCCCov[lambda_List,p_Integer,MM_Integer,a_Association]:=RandULPCCCov[lambda,p,MM,a]=
(*Step-3 weights: covering sets with exactly this size profile, by the same deleted-cycle
inclusion-exclusion. Summed over the profiles these give RandULPCFixCov back, which the suite
checks.*)
Total[Function[sm,
	sm[[2]]*(Times@@KeyValueMap[
		Function[{d,m},Binomial[Lookup[RandULPCNd[sm[[1]],p],d,0],m]],a])
	]/@RandULPCSubMulti[lambda]];

(* Private helper *)
RandULPCComps[p_Integer,sizes_List,aVec_List,alphaVec_List,lens_List,mv_List,uv_List]:=
	RandULPCComps[p,sizes,aVec,alphaVec,lens,mv,uv]=
(*The number of ways to finish a partial draw, given how many orbits of each size are already used
(alphaVec) and which cycle lengths are still uncovered (uv). Already-chosen orbits sit inside the
covered cycles, so they avoid every deleted set and survive into every term -- only their count per
size enters, which is what keeps the state this small. Indexed by the uncovered LENGTH MULTISET
rather than the uncovered cycle set, since the orbit counts of a reduced cycle type depend only on
its shape: states of equal shape share an entry, and the sum runs over sub-multisets instead of all
2^(uncovered cycles) subsets.*)
Module[{tot=0,surv,ndm},
	Do[
		surv=mv-tv;
		ndm=RandULPCNd[Sort[Flatten[MapThread[ConstantArray[#2,#1]&,{surv,lens}]]],p];
		tot+=(-1)^Total[tv]*(Times@@MapThread[Binomial,{uv,tv}])*
			(Times@@MapThread[Function[{d,al,aa},Binomial[Lookup[ndm,d,0]-al,aa-al]],
				{sizes,alphaVec,aVec}])
	,{tv,Tuples[Range[0,#]&/@uv]}];
	tot
];

(* Private helper *)
RandULPCDraw[lambda_List,p_Integer,MM_Integer,a_Association]:=
(*Step 4. The orbits are drawn one at a time, each weighted by the number of ways the draw can
still be completed, which is exact by the chain rule. Candidates are bucketed by (size, support
meeting the uncovered cycles): orbits sharing that pair lead to the same state, so one completion
count serves a whole bucket. The draw produces an ORDERED sequence and forgets the order -- every
valid set has the same number of orderings, so the set is uniform and the ordering factor cancels
out of the weights.*)
Module[{orbs=RandULPCOrbits[lambda,p],nOrb,lens,mv,uv,uncov,sizes,aVec,alphaVec,sIdx,
		used,chosen={},tot,buckets,keys,w,kpick,opick},
	nOrb=Length[orbs];
	lens=Sort[DeleteDuplicates[lambda]];
	mv=Table[Count[lambda,k],{k,lens}];
	uv=mv;
	uncov=Range[Length[lambda]];
	sizes=Sort[Keys[a]];
	aVec=a[#]&/@sizes;
	sIdx=Association[Thread[sizes->Range[Length[sizes]]]];
	alphaVec=ConstantArray[0,Length[sizes]];
	used=ConstantArray[False,nOrb];
	tot=Total[aVec];
	Do[
		buckets=<||>;
		Do[
			If[!used[[i]]&&alphaVec[[sIdx[orbs[[i,1]]]]]<a[orbs[[i,1]]],
				With[{key={orbs[[i,1]],Intersection[orbs[[i,2]],uncov]}},
					buckets[key]=Append[Lookup[buckets,Key[key],{}],i]]]
		,{i,nOrb}];
		keys=Keys[buckets];
		w=Table[
			Length[buckets[key]]*RandULPCComps[p,sizes,aVec,
				MapAt[#+1&,alphaVec,sIdx[key[[1]]]],lens,mv,
				uv-Table[Count[lambda[[key[[2]]]],L],{L,lens}]]
		,{key,keys}];
		kpick=RandomChoice[w->keys];
		opick=RandomChoice[buckets[kpick]];
		used[[opick]]=True;
		AppendTo[chosen,opick];
		alphaVec=MapAt[#+1&,alphaVec,sIdx[kpick[[1]]]];
		uv=uv-Table[Count[lambda[[kpick[[2]]]],L],{L,lens}];
		uncov=Complement[uncov,kpick[[2]]]
	,{tot}];
	Catenate[orbs[[#,3]]&/@chosen]
];

(* Private helper *)
RandULPCOne[p_Integer,MM_Integer,n_Integer]:=
(*One uniform isomorphism class, by the five steps above, as a facet list.*)
Module[{tw,lambda,profs,a,facets,relabel},
	tw=RandULPCTypeWeights[p,MM,n];
	lambda=RandomChoice[tw[[2]]->tw[[1]]];
	profs=RandULPCProfiles[lambda,p,MM];
	a=RandomChoice[(RandULPCCCov[lambda,p,MM,#]&/@profs)->profs];
	facets=RandULPCDraw[lambda,p,MM,a];
	relabel=RandomSample[Range[n]];
	Sort[Sort/@Map[relabel[[#]]&,facets]]
];


(* ::Item::Closed:: *)
(*RandomUniformUnlabeledPureSimplicialComplex*)


(* Primary Pattern *)

RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}, numSamples_Integer]:=
(*Gives numSamples fully unlabeled pure simplicial complexes of purity p, facet order M and vertex
count nV, each a facet list. Uniform over the isomorphism classes counted by
NumUnlabeledPureComplexes, and produced without rejection.*)
Catch[

	If[numSamples<0,
		Message[RandomUniformUnlabeledPureSimplicialComplex::argerr,{purity,facetOrder,nV},numSamples];
		Throw[$Failed, "ECGravReturn$65"]];

	If[purity<1||facetOrder<1||NumUnlabeledPureComplexes[purity,facetOrder,nV]==0,
		Message[RandomUniformUnlabeledPureSimplicialComplex::empty,purity,facetOrder,nV];
		Throw[$Failed, "ECGravReturn$65"]];

	If[RandULPCGoParallel[numSamples],
		ParallelTable[RandULPCOne[purity,facetOrder,nV],{numSamples},
			DistributedContexts->{$Context,"ECGrav`","ECGrav`Private`"}],
		Table[RandULPCOne[purity,facetOrder,nV],{numSamples}]]

, "ECGravReturn$65"];


(* Overload Pattern *)

RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}, numSamples_Integer]:=
(*The same with the vertex count free. The sample space is the unlabeled complexes, so it is
NumUnlabeledPureComplexes that supplies the vertex-count distribution; the old body drew n from the
VERTEX-labeled counts and let rejection repair the difference, which is not the distribution that
was wanted.*)
Module[{counts,ns},
Catch[

	If[numSamples<0,
		Message[RandomUniformUnlabeledPureSimplicialComplex::argerr,{purity,facetOrder},numSamples];
		Throw[$Failed, "ECGravReturn$66"]];

	If[purity<1||facetOrder<1,
		Message[RandomUniformUnlabeledPureSimplicialComplex::empty,purity,facetOrder,"any"];
		Throw[$Failed, "ECGravReturn$66"]];

	ns=Range[purity,purity*facetOrder];
	counts=NumUnlabeledPureComplexes[purity,facetOrder,#]&/@ns;

	If[Total[counts]==0,
		Message[RandomUniformUnlabeledPureSimplicialComplex::empty,purity,facetOrder,"any"];
		Throw[$Failed, "ECGravReturn$66"]];

	If[RandULPCGoParallel[numSamples],
		ParallelTable[RandULPCOne[purity,facetOrder,RandomChoice[counts->ns]],{numSamples},
			DistributedContexts->{$Context,"ECGrav`","ECGrav`Private`"}],
		Table[RandULPCOne[purity,facetOrder,RandomChoice[counts->ns]],{numSamples}]]

, "ECGravReturn$66"]
];


(* Overload Pattern *)

(*First only if there is something to take: the sampling forms return $Failed on an empty
sample space, and First[$Failed] would bury that under a Part error.*)
RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer,nV_Integer}]:=
With[{sample=RandomUniformUnlabeledPureSimplicialComplex[{purity,facetOrder,nV}, 1]},
	If[ListQ[sample],First[sample],sample]];


(* Overload Pattern *)

(*First only if there is something to take: the sampling forms return $Failed on an empty
sample space, and First[$Failed] would bury that under a Part error.*)
RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
With[{sample=RandomUniformUnlabeledPureSimplicialComplex[{purity,facetOrder}, 1]},
	If[ListQ[sample],First[sample],sample]];


(* Catch-all Pattern *)
RandomUniformUnlabeledPureSimplicialComplex[args___]:=(Message[RandomUniformUnlabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random uniform Facet-labeled Pure Simplicial Complex by Burnside pair sampling*)


(* Exact uniform sampling of facet-labeled complexes. -----------------------------------------
   RandomVertexLabeledPureSimplicialComplex reads its weights straight off the counting
   recursion, because with labeled vertices "how many complexes have this many new vertices in
   the next facet" is a plain product of binomials. Nothing like that exists here: with the
   vertices unlabeled the count is a Burnside average over cycle types, not a sequential
   decomposition, so there is no recursion to walk facet by facet.

   What replaces it is pair sampling. Let X be the ordered M-tuples of distinct p-subsets of [n]
   covering [n], on which S_n acts by relabeling vertices; facet-labeled complexes are its
   orbits. Sample uniformly from the FIXED PAIRS {(sigma,x) : sigma.x = x}. Every orbit O
   contributes Sum over x in O of |Stab(x)| = |O| (n!/|O|) = n! pairs -- the same number for
   every orbit, regardless of its size -- so discarding sigma leaves the orbit exactly uniform.
   No rejection, and no automorphism-order weighting to go wrong.

   sigma fixes a tuple exactly when every facet is sigma-invariant, i.e. a union of cycles whose
   lengths sum to p. Calling the cycles "blocks", one sample is:

	1  draw a cycle type with weight (number of permutations of that type) x (number of
	   tuples it fixes). Only partitions of n into parts <= p occur: a cycle longer than p
	   fits in no facet, so it could never be covered.
	2  cut a uniform random permutation of [n] into blocks of those sizes
	3  draw the M distinct weight-p block-subsets one at a time, each weighted by how many
	   ways the rest can still be finished
	4  facet i is the union of the blocks in subset i

   Step 3 is exact for a reason worth stating: the number of completions depends on the draws so
   far only through (how many subsets are used, which blocks are still uncovered), never on WHICH
   subsets were used. Every used subset lies inside the already-covered blocks, so it is excluded
   from any inclusion-exclusion term automatically and only its count enters. That is the same
   observation that lets the vertex-labeled generator subtract (j-1) for a repeated facet.

   The state is indexed by the multiplicity vector of uncovered block sizes rather than by the
   set of uncovered blocks. Blocks of equal size are interchangeable at this point, so that turns
   a 2^(number of blocks) inclusion-exclusion into O(n^p) -- the same polynomial scaling, and the
   same independence of M, that NumFacetLabeledPureComplexes has. *)

(* Private helper *)
RandFLPCWeightCounts[nu_List]:=RandFLPCWeightCounts[nu]=
(*Coefficients of Product[(1+x^k)^nu[[k]],{k,1,p}] out to x^p. Entry w+1 counts the sub-multisets
of a block collection holding nu[[k]] blocks of size k whose sizes sum to w; the last entry is
therefore the number of admissible facets. Built by truncated convolution, so the degree-n
polynomial is never formed.*)
Module[{p=Length[nu],acc},
	acc=PadRight[{1},p+1];
	Do[
		acc=Table[
			Sum[If[Divisible[d-i,k],acc[[i+1]]*Binomial[nu[[k]],(d-i)/k],0],{i,0,d}]
		,{d,0,p}]
	,{k,1,p}];
	acc
];

(* Private helper *)
RandFLPCCompletions[lambda_List,MM_Integer,j_Integer,nu_List]:=
	RandFLPCCompletions[lambda,MM,j,nu]=
(*Ways to finish a partly drawn sample: choose MM-j further weight-p subsets, pairwise distinct
and none of the j already spent, covering every block still uncovered. Inclusion-exclusion over
which uncovered blocks are left uncovered -- for a given omega the admissible facets are the
weight-p subsets avoiding omega, and j of those are already spent, since every used subset sits
inside the covered blocks and so avoids omega for free.*)
Module[{p=Length[lambda]},
	Total[Table[
		(-1)^Total[om]*(Times@@Binomial[nu,om])*
			FactorialPower[Last[RandFLPCWeightCounts[lambda-om]]-j,MM-j]
	,{om,Tuples[Range[0,#]&/@nu]}]]
];

(* Private helper *)
RandFLPCTypeWeights[p_Integer,MM_Integer,n_Integer]:=RandFLPCTypeWeights[p,MM,n]=
(*The contributing cycle types, as both partitions and multiplicity vectors, with their Burnside
weights. The weights sum to n! NumFacetLabeledPureComplexes[p,MM,n], which the test suite checks.*)
Module[{types,lambdas},
	types=IntegerPartitions[n,All,Range[p]];
	lambdas=Table[Count[t,#]&/@Range[p],{t,types}];
	{types,lambdas,
		Table[
			(n!/Times@@Table[k^lam[[k]]*lam[[k]]!,{k,1,p}])*RandFLPCCompletions[lam,MM,0,lam]
		,{lam,lambdas}]}
];

(* Private helper: sample count above which the facet-labeled draw goes to ParallelTable. -----
   Kept separate from $RandomComplexParallelThreshold because the per-sample costs differ by
   more than an order of magnitude: a vertex-labeled sample is tens of microseconds, a
   facet-labeled one a few hundred, so parallel pays off far sooner here. Measured crossover on
   11 subkernels with kernels already up: 50 samples for {3,4,8}, 100 for {2,3,4} and {2,5,8},
   200 for {3,3,5}. At 100 the ratios run 0.93 to 1.77 and by 200 they are all above 1, so 100
   sits at the crossover. Gated on $KernelCount as everywhere else, so a small call never
   quietly launches kernels. *)

$RandomFacetLabeledParallelThreshold=100;

(* Private helper *)
RandFLPCGoParallel[numSamples_Integer]:=
	numSamples>$RandomFacetLabeledParallelThreshold&&$KernelCount>0;

(* Private helper *)
RandFLPCOne[p_Integer,MM_Integer,n_Integer]:=
(*One uniform facet-labeled complex, by the four steps above.*)
Module[{types,lambdas,ws,pick,sizes,lambda,blocks,cands,used={},uncov,nuOf,avail,weights,choice,facets},

	{types,lambdas,ws}=RandFLPCTypeWeights[p,MM,n];
	pick=RandomChoice[ws->Range[Length[types]]];
	sizes=types[[pick]];
	lambda=lambdas[[pick]];

	(*a uniform partition of [n] into blocks of these sizes*)
	blocks=TakeList[RandomSample[Range[n]],sizes];

	(*every admissible facet, as a set of block indices; blocks have size >= 1, so at most p*)
	cands=Select[Subsets[Range[Length[sizes]],{1,p}],Total[sizes[[#]]]==p&];

	nuOf[idxs_]:=Count[sizes[[idxs]],#]&/@Range[p];
	uncov=Range[Length[sizes]];

	facets=Table[
		avail=Complement[cands,used];
		weights=Table[RandFLPCCompletions[lambda,MM,j+1,nuOf[Complement[uncov,S]]],{S,avail}];
		choice=RandomChoice[weights->avail];
		AppendTo[used,choice];
		uncov=Complement[uncov,choice];
		Sort[Catenate[blocks[[choice]]]]
	,{j,0,MM-1}];

	facets
];


(* ::Item::Closed:: *)
(*RandomUniformFacetLabeledPureSimplicialComplex*)


(* Primary Pattern *)

RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 04/29/2026. *) *)
(*   (* Note:  *) 
(****************************************)
(*A code to generate a random sample of (possibly disconnected) unlabeled pure 
simplicial complexes of a given purity, facet order, and number of vertices by 
rejection based on number of ways to label from vertex-labeled complexes. Outputs a 
list numSamples of them. Input is a list of three numbers, p,q,n, 
numSamples, where p is the purity number, q is the facet order, n is the number of vertices
and numSamples is the number of samples required. E.g. 
RandomUniformFacetLabeledPureSimplicialComplex[{2,4,6},10] generates 10, 2-pure facet-labeled
random simplicial complex with 4 edges and 6 vertices with a uniform distribution.     *)
*)
*)
Catch[

If[numSamples<0,
	Message[RandomUniformFacetLabeledPureSimplicialComplex::argerr,{purity,facetOrder,nV},numSamples];
	Throw[$Failed, "ECGravReturn$67"]];

(*The sample space is the facet-labeled complexes, so it is NumFacetLabeledPureComplexes that
has to be non-empty, not the vertex-labeled count. The two vanish together here -- both need
p <= n <= p M -- but the facet-labeled one also needs Binomial[n,p] >= M distinct facets, so it
is the right test.*)
If[purity<1||facetOrder<1||NumFacetLabeledPureComplexes[purity,facetOrder,nV]==0,
	Message[RandomUniformFacetLabeledPureSimplicialComplex::empty,purity,facetOrder,nV];
	Throw[$Failed, "ECGravReturn$67"]];

(*Exact Burnside pair sampling; see RandFLPCOne above. The old body proposed from the
vertex-labeled generator and accepted with probability Min[1,(Min[M,n]!/n!) x facet stabiliser
order], which is correct but collapses as n grows -- 15 s per hundred samples at {3,4,6}. There
is no rejection here at all.*)
If[RandFLPCGoParallel[numSamples],
	ParallelTable[RandFLPCOne[purity,facetOrder,nV],{numSamples},
		DistributedContexts->{$Context,"ECGrav`","ECGrav`Private`"}],
	Table[RandFLPCOne[purity,facetOrder,nV],{numSamples}]]

, "ECGravReturn$67"];


(* Overload Pattern *)

RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 05/01/2026. *) *)
(*   (* Note: . *) 
(****************************************)
(*A code to generate a random sample of (possibly disconnected)vertex labeled pure 
simplicial complexes of a given purity and clique order through iterative addition of 
facets. Outputs a list numSamples of them. Input is a list of three numbers, p,q, 
numSamples, where p is the purity number, q is the facet order, and numSamples is 
the number of samples required. numSamples is set to 1 by default. E.g. 
RandomUniformFacetLabeledPureSimplicialComplex[{2,4},10] generates 10, 2-pure random 
facet-labeled simplicial complex with 4 edges.     *)
*)
*)
Catch[

If[numSamples<0,
	Message[RandomUniformFacetLabeledPureSimplicialComplex::argerr,{purity,facetOrder},numSamples];
	Throw[$Failed, "ECGravReturn$68"]];

If[purity<1||facetOrder<1,
	Message[RandomUniformFacetLabeledPureSimplicialComplex::empty,purity,facetOrder,"any"];
	Throw[$Failed, "ECGravReturn$68"]];

(*With the vertex count free, draw it from the facet-labeled counts themselves and then sample
at that n. The old body drew n from the VERTEX-labeled counts and let the rejection step repair
the difference; weighting by NumFacetLabeledPureComplexes is the distribution that was wanted
all along, and needs no repair.*)
With[{verts=Range[purity,purity*facetOrder],
		weights=Table[NumFacetLabeledPureComplexes[purity,facetOrder,n],{n,purity,purity*facetOrder}]},
	If[RandFLPCGoParallel[numSamples],
		ParallelTable[RandFLPCOne[purity,facetOrder,RandomChoice[weights->verts]],{numSamples},
			DistributedContexts->{$Context,"ECGrav`","ECGrav`Private`"}],
		Table[RandFLPCOne[purity,facetOrder,RandomChoice[weights->verts]],{numSamples}]]]

, "ECGravReturn$68"];


(* Overload Pattern *)

(*First only if there is something to take: the sampling forms return $Failed on an empty
sample space, and First[$Failed] would bury that under a Part error.*)
RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer,nV_Integer}]:=
With[{sample=RandomUniformFacetLabeledPureSimplicialComplex[{purity,facetOrder,nV}, 1]},
	If[ListQ[sample],First[sample],sample]];


(* Overload Pattern *)

(*First only if there is something to take: the sampling forms return $Failed on an empty
sample space, and First[$Failed] would bury that under a Part error.*)
RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
With[{sample=RandomUniformFacetLabeledPureSimplicialComplex[{purity,facetOrder}, 1]},
	If[ListQ[sample],First[sample],sample]];


(* Catch-all Pattern *)
RandomUniformFacetLabeledPureSimplicialComplex[args___]:=(Message[RandomUniformFacetLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random facet labeled Pure Simplicial Complex (not uniform)*)


(* ::Item::Closed:: *)
(*RandomFacetLabeledPureSimplicialComplex*)


(* :Code Section *)

(* Primary Pattern *)
RandomFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
(*
(****************************************)
(*   (* Last updated 08/15/2024. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*A code to generate a random pure facet-labeled(possibly disconnected) simplicial 
complex of a given purity and clique order through iterative addition of vertices 
labeled by the facets. Input is a list of two numbers, {p, q} where p is the purity 
number and q is the number of facets. E.g. 
ECGrav`RandomFacetLabeledPureSimplicialComplex[{2,4}] generates a 2-pure facet-labeled 
random simplicial complex with 4 edges. It successively choses a random vertex from a 
decreasing list of available vertex space. The resulting complexes are not generated 
from a uniform distribution over all facet-labeled complexes. *)

*)
Module[{complex={},result,curAvailableVertices=Subsets[Range[facetOrder]][[2;;-1]],
	curUnAvailableVertices={},vertexDegree,newVertex={}},

(*randomSubset[lst_List]:=Map[If[RandomReal[]<=1/2,#,##&[]]&,lst];*)

vertexDegree=<|Table[i->0,{i,curAvailableVertices}]|>;

While[curAvailableVertices!={},

(*Print[""];
Print["In while loop"];
Print["Before updating, complex ",complex," vertexDegree, ",vertexDegree, " curAvailableVertices ",curAvailableVertices, " curUnAvailableVertices ",curUnAvailableVertices];*)

newVertex=RandomChoice[curAvailableVertices];

Do[vertexDegree[[Key[i]]]++,{i,Subsets[newVertex][[2;;-1]]}];
complex=Join[complex,{newVertex}];

curUnAvailableVertices=Union[curUnAvailableVertices,
	Table[
		If[(Length[i]==1&&vertexDegree[[Key[i]]]==purity)||
			(Length[i]>1&&vertexDegree[[Key[i]]]==purity-1),i,Nothing],{i,curAvailableVertices}]];

curAvailableVertices=Fold[Function[{x,y},Select[x,!SubsetQ[#,y]&]],curAvailableVertices,curUnAvailableVertices];
(*Print[""];
Print["After updating, newVertex ",newVertex," complex ",complex];
Print[" vertex degrees ",vertexDegree];

Print[" curUnAvailableVertices ",curUnAvailableVertices," curAvailableVertices ",curAvailableVertices];*)

];


complex=Sort[complex];

(*Join@@(Map[LexicographicSort[#]&,
GatherBy[ReverseSortBy[Table[Keys[Select[cmlxAsn,MemberQ[#,v]&]],{v,vlist}],Length],Length],{1}])
*)
result=With[{cmlxAsn=<|Table[i->complex[[i]],{i,1,Length[complex]}]|>},
Table[With[{k=Select[complex,MemberQ[#,q]&]},
Keys[Select[cmlxAsn,MemberQ[k,#]&]]
],{q,1,facetOrder}]
];


(*result=Table[With[{k=Select[complex,MemberQ[#,v]&]},


Flatten[SubsetPosition[complex,k]]],{v,1,facetOrder}];*)


Sort[Sort/@result]

];

(* Overload Pattern *)
RandomFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer},numSamples_Integer]:=
If[numSamples<10^4,
	Table[RandomFacetLabeledPureSimplicialComplex[{purity,facetOrder}],{numSamples}],
	ParallelTable[RandomFacetLabeledPureSimplicialComplex[{purity,facetOrder}],{numSamples},DistributedContexts->{$Context,"ECGrav`Private`"}]
];

(* Catch-all Pattern *)
RandomFacetLabeledPureSimplicialComplex[args___]:=(Message[RandomFacetLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection:: *)
(*MCMC Random Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexIterMCMCSweep*)


(* :Code Section *)

(* Primary Pattern *)
RandomPureSimplicialComplexMCMCSweep[seedComplex_List,labelingChoise_Integer,NN_Integer]:=
(*
(****************************************)
(*   (* Last updated 05/14/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*A code to perform NN MC sweeps on a seed pure simplicial 
complex through MC chain. Inputs are a seed complex, labelingChoise (which is 
0 for vertex-labeled, 1 for facet-labeled, and 2 for unlabeled), and the number of 
sweeps.

Outputs an association with the probability, number of vertices, and last complex *)

*)
Module[{purity=Length[seedComplex[[1]]],facetOrder=Length[seedComplex],data,
	allVertexLabels,step},

allVertexLabels=Range[purity*facetOrder];

data=<|"complex"->Sort[Sort/@seedComplex],
		"vertexCount"->Length[DeleteDuplicates[Flatten[seedComplex]]],
		"weight"->0.0,"energy"->0.0|>;

data[[Key["weight"]]]=
	Which[labelingChoise==0,
			1.0/Binomial[purity*facetOrder,data[[Key["vertexCount"]]]],
		labelingChoise==1,
			(1.0/Binomial[purity*facetOrder,data[[Key["vertexCount"]]]])
				*PureComplexFacetStabilizerGroupOrder[seedComplex]*(facetOrder!/(data[[Key["vertexCount"]]]!)),
		labelingChoise==2,
			(1.0/Binomial[purity*facetOrder,data[[Key["vertexCount"]]]])
				*PureComplexAutomorphismGroupOrder[seedComplex]/((data[[Key["vertexCount"]]])!)];

data[[Key["energy"]]]=-Log[data[[Key["weight"]]]];

step[]:=(*Performs one spin flip step*)
Module[{
	curComplexWeight=data[[Key["weight"]]],
	curVCount=data[[Key["vertexCount"]]],
	ranFacet=RandomChoice[data[[Key["complex"]]]],
	intermediateComplex,
	newFacet,
	newComplex,
	newVCount,
	newComplexWeight,
	selectionProb,
	accept},

	intermediateComplex=Complement[data[[Key["complex"]]],{ranFacet}];


	newFacet=Sort[RandomSample[allVertexLabels,purity]];

	While[MemberQ[data[[Key["complex"]]],newFacet],
		newFacet=Sort[RandomSample[allVertexLabels,purity]]
	];

	newComplex=Join[intermediateComplex,{newFacet}];
	newVCount=Length[DeleteDuplicates[Flatten[newComplex]]];

	Which[labelingChoise==0,
			newComplexWeight=1.0/Binomial[purity*facetOrder,newVCount],
		labelingChoise==1,
			newComplexWeight=(1.0/Binomial[purity*facetOrder,newVCount])
				*PureComplexFacetStabilizerGroupOrder[newComplex]*(facetOrder!/newVCount!),
		labelingChoise==2,
			newComplexWeight=(1.0/Binomial[purity*facetOrder,newVCount])
				*PureComplexAutomorphismGroupOrder[newComplex]/(newVCount!)];

	selectionProb=Min[1.0,1.0*newComplexWeight/curComplexWeight];

	accept = 0;

	If[RandomReal[]<selectionProb,accept =1];

	If[accept==1,
		data[[Key["complex"]]]=Sort[newComplex];
		data[[Key["vertexCount"]]]=newVCount;
		data[[Key["weight"]]]=newComplexWeight;
		data[[Key["energy"]]]= -Log[1.0*newComplexWeight];
	];

];


Do[step[],{NN*purity* facetOrder}];

data
]


(* Overload Pattern *)

RandomPureSimplicialComplexMCMCSweep[seedComplex_List,HoldNumberOfVerticesFixed_?BooleanQ,labelingChoise_Integer,NN_Integer]:=
(*
(****************************************)
(*   (* Last updated 06/12/2026. *) *)
(*   (* Version: 1 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*A code to perform NN MC sweeps on a seed pure simplicial 
complex through MC chain. Inputs are a seed complex, HoldNumberOfVerticesFixed is a 
Boolean with True keeping the number of vertices fixed, labelingChoise (which is 
0 for vertex-labeled, 1 for facet-labeled, and 2 for unlabeled), and the number of 
sweeps.

Outputs an association with the probability, number of vertices, and last complex *)

*)
Module[{result,purity=Length[seedComplex[[1]]],
		facetOrder=Length[seedComplex],
		nV=Length[DeleteDuplicates[Flatten[seedComplex]]],
		minNv,data,allVertexLabels,step},

Catch[If[HoldNumberOfVerticesFixed==False,Throw[RandomPureSimplicialComplexMCMCSweep[seedComplex,labelingChoise,NN], "ECGravReturn$54"]];

minNv=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]];

allVertexLabels=Union@@seedComplex;

data=<|"complex"->Sort[Sort/@seedComplex],"edgeCount"->Length[Union@@(Subsets[#,{2}]&/@(Sort/@seedComplex))],"weight"->0.0,"energy"->0.0|>;

data[[Key["weight"]]]=
	Which[
		labelingChoise==0,
			1.0,
		labelingChoise==1,
			PureComplexFacetStabilizerGroupOrder[seedComplex]*(1.0*(facetOrder!)/(nV!)),
		labelingChoise==2,
			1.0*PureComplexAutomorphismGroupOrder[seedComplex]]/(nV!);

data[[Key["energy"]]]=-Log[data[[Key["weight"]]]];

step[]:=(*Performs one spin flip step*)
	Block[{
		curComplexWeight=data[[Key["weight"]]],
		ranFacetPairs=RandomSample[data[[Key["complex"]]],2],
		ranFacetPairVertices,
		vertexPairSwapRule,
		intermediateComplex,
		newFacetPairs,
		newComplex,
		newComplexWeight,
		selectionProb,
		accept},

	intermediateComplex=Complement[data[[Key["complex"]]],ranFacetPairs];

	ranFacetPairVertices=Union@@ranFacetPairs;
	vertexPairSwapRule=Thread[ranFacetPairVertices->RandomSample[ranFacetPairVertices]];

	newFacetPairs=Sort/@(ranFacetPairs/.vertexPairSwapRule);

	While[Intersection[newFacetPairs,intermediateComplex]!={},

		vertexPairSwapRule=Thread[ranFacetPairVertices->RandomSample[ranFacetPairVertices]];

		newFacetPairs=Sort/@(ranFacetPairs/.vertexPairSwapRule);

	];

	newComplex=Sort[Join[intermediateComplex,newFacetPairs]];


	Which[
		labelingChoise==0,
			newComplexWeight=1.0,
		labelingChoise==1,
			newComplexWeight=PureComplexFacetStabilizerGroupOrder[newComplex]*(1.0*(facetOrder!)/(nV!)),
		labelingChoise==2,
			newComplexWeight=1.0*PureComplexAutomorphismGroupOrder[newComplex]/(nV!)];

	selectionProb=Min[1.0,1.0*newComplexWeight/curComplexWeight];

	accept = 0;

	If[RandomReal[]<selectionProb,accept =1];

	If[accept==1,

		data[[Key["complex"]]]=Sort[newComplex];
		data[[Key["edgeCount"]]]=Length[Union@@(Subsets[#,{2}]&/@(newComplex))];
		data[[Key["weight"]]]=newComplexWeight;
		data[[Key["energy"]]]= -Log[1.0*newComplexWeight];
		];

	];


	Do[step[];
	,{i,NN*purity* facetOrder}];

	result=data;

	result, "ECGravReturn$54"]
];


(* Catch-all Pattern *)
RandomPureSimplicialComplexMCMCSweep[args___]:=(Message[RandomPureSimplicialComplexMCMCSweep::argerr, args];
$Failed);


(* ::Item:: *)
(*RandomPureSimplicialComplexMCMCEquilibriate*)


(* :Code Section *)

(* Primary Pattern *)
RandomPureSimplicialComplexMCMCEquilibriate[seedComplex_List, labelingChoise_Integer]:=
(*
(****************************************)
(*   (* Last updated 05/05/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*runs enough MC sweeps to equilibriate the input seedComplex. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. Equilibriation 
	is established when four different complexes (the input, the complex with the minumum 
number	 of vertices, the complex with the maximum number of vertices, and another random 
complex	) all start to have close enough probabilities such that the fluctions of each is
	of the order of the fluctuation across them.*)

*)

Module[{result,purity=Length[seedComplex[[1]]],facetOrder=Length[seedComplex],data,
	sweepOutput,Entable,outWinLength,inWinLength,AllEnMat,sqMeanEMat,
	sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime=20000,maxNumSweeps=25000},

data=<|"state"-><|"complex"->Sort[Sort/@seedComplex],
				"vertexCount"->Length[DeleteDuplicates[Flatten[seedComplex]]],
				"weight"->1.0,"energy" ->0.0|>,
		"maxNumV"-><|"complex"->Partition[Range[purity*facetOrder],purity],
				"vertexCount"->purity*facetOrder,"weight"->1.0,"energy" ->0.0|>,
		"minNumV"-><|"complex"->With[{minN=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]]},RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder,minN}]],
				"vertexCount"->Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]],"weight"->1.0,"energy" ->0.0|>,
		"random"-><|"complex"->RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder}],
				"vertexCount"->0,"weight"->1.0,"energy" ->0.0|>|>;

data[["random","vertexCount"]]=Length[DeleteDuplicates[Flatten[data[["random","complex"]]]]];


Which[
	labelingChoise==1,(*facet-labeled*)
		data[["state","weight"]]=PureComplexFacetStabilizerGroupOrder[seedComplex]*(facetOrder!/(data[["state","vertexCount"]]!));
		data[["maxNumV","weight"]]=((purity!)^facetOrder)*(facetOrder!/(purity*facetOrder)!);
		data[["minNumV","weight"]]=PureComplexFacetStabilizerGroupOrder[data[["minNumV","complex"]]]*(facetOrder!/data[["minNumV","vertexCount"]]!);
		data[["random","weight"]]=PureComplexFacetStabilizerGroupOrder[data[["random","complex"]]]*(facetOrder!/data[["random","vertexCount"]]!),
	labelingChoise==2,(*unlabeled*)	
		data[["state","weight"]]=PureComplexAutomorphismGroupOrder[seedComplex]/((data[["state","vertexCount"]])!);
		data[["maxNumV","weight"]]=(facetOrder!*(purity!)^facetOrder)/((purity*facetOrder)!);
		data[["minNumV","weight"]]=PureComplexAutomorphismGroupOrder[data[["minNumV","complex"]]]/(data[["minNumV","vertexCount"]]!);
		data[["random","weight"]]=PureComplexAutomorphismGroupOrder[data[["random","complex"]]]/(data[["random","vertexCount"]]!);
];

data[["state","energy"]]=-Log[1.0*data[["state","weight"]]];
data[["maxNumV","energy"]]=-Log[1.0*data[["maxNumV","weight"]]];
data[["minNumV","energy"]]=-Log[1.0*data[["minNumV","weight"]]];
data[["random","energy"]]=-Log[1.0*data[["random","weight"]]];


numsweeps=0;
outWinLength =400 ;
	(* length of a table to store running energy values to test equilibriation*)
inWinLength=10;
	(* A segment to be averaged over at the beginning and the end of the table*)

Entable=Table[0.0,{m,4},{n,outWinLength}];(*A table to store energy values of the seed, random and empty tracks for tunning calculation of tests of equilibriation*)

(*Main calculation*)
AllEnMat=
	(*AllEnMat is for storing all energy values for diagnostic*)
Reap[
	While[numsweeps<maxNumSweeps,
		numsweeps++;

		If[Mod[numsweeps,100]==0,
			PrintTemporary["Running equilibriation at sweep number ",numsweeps, " complex ",data[[Key["state"],Key["complex"]]]]
		];

		Sow[{data[[Key["state"],Key["energy"]]],data[[Key["maxNumV"],Key["energy"]]],data[[Key["minNumV"],Key["energy"]]],data[[Key["random"],Key["energy"]]]}];
		Do[

			sweepOutput=RandomPureSimplicialComplexMCMCSweep[data[[Key[replicaName],Key["complex"]]],
				labelingChoise,1];

			data[[Key[replicaName]]]=sweepOutput;

		,{replicaName,{"state","maxNumV","minNumV","random"}}];

		Entable[[All,1]]={data[[Key["state"],Key["energy"]]],
			data[[Key["maxNumV"],Key["energy"]]],data[[Key["minNumV"],Key["energy"]]],data[[Key["random"],Key["energy"]]]};

		If[numsweeps>outWinLength,

			sqMeanEMat=Table[(Mean[Entable[[i]][[1;;inWinLength]]]-Mean[Entable[[j]][[(outWinLength-inWinLength);;outWinLength]]])^2,{i,4},{j,4}];
				(*sqMeanEMat is a four by cour matrix of the squared means of difference in energy within and across the 
					tracks at the beginning and end of outWinLength. At equilibrium these 12 mumbers should be 
					randomly distributed with mean of 0. Their fluctuation from 0 should be within the variance of the 
					newer(late time) variance in energy for the equilibriation to exit. 
					sqMeanPairwiseDiff is the mean of these 12 numbers*)

			sqMeanPairwiseDiff=Mean[Flatten[sqMeanEMat]];

			meanLateVar=Mean[Table[Variance[Entable[[i]][[1;;inWinLength]]],{i,4}]];(*Mean variance of the newer energies*)

			If[Abs[sqMeanPairwiseDiff]<meanLateVar,
				eqlTime=numsweeps-outWinLength+inWinLength;
				Break[]
			];(*Exit the while loop and go to the Do loop*)
			If[Abs[Tr[sqMeanEMat]]<0.00001,	(*Print["exiting because stuck in a metastable state, sqMeanEMat is ",MatrixForm[sqMeanEMat]];*)eqlTime=numsweeps-outWinLength+inWinLength;
				Break[]
			](*Exit the while loop and go to the Do loop*)
		];

		Entable=RotateRight[#,1]&/@Entable;(*Shifts every entry to the right by one cyclically. The new data will be written on the first slot so newer data is at the beginning.*)
	];
][[2,1]];
(********************
*  For diagnostics *
********************)
(*Print["AllEnMat ",Transpose[AllEnMat]];*)

Print[ListLinePlot[Transpose[AllEnMat[[1;;Min[Length[AllEnMat],2*eqlTime]]]],
		PlotRange->All,
	PlotLabel->"t from 1 to 2 times eqlT ",PlotLegends->{"seedComplex","maxNumV","minNumV","random"}]
];


result=<| "eqlT"->eqlTime,"state"->data[[Key["state"]]]|>;

result

]


(* Overload Pattern *)

RandomPureSimplicialComplexMCMCEquilibriate[seedComplex_List, HoldNumberOfVerticesFixed_?BooleanQ, labelingChoise_Integer]:=
(*
(****************************************)
(*   (* Last updated 06/12/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*runs enough MC sweeps to equilibriate the input seedComplex. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. Equilibriation 
	is established when four different complexes (the input, the complex with the minumum 
number	 of vertices, the complex with the maximum number of vertices, and another random 
complex	) all start to have close enough probabilities such that the fluctions of each is
	of the order of the fluctuation across them.*)

*)

Module[{result,
		purity=Length[seedComplex[[1]]],
		facetOrder=Length[seedComplex],
		nV=Length[DeleteDuplicates[Flatten[seedComplex]]],
		makeLeastEvenlyConnectedComplex,
		data,sweepOutput,
		Entable,outWinLength,inWinLength,AllEnMat,
		sqMeanEMat,
		sqMeanPairwiseDiff,
		meanLateVar,
		numsweeps,
		eqlTime=20000,
		maxNumSweeps=25000},

Catch[If[HoldNumberOfVerticesFixed==False,
	Throw[RandomPureSimplicialComplexMCMCEquilibriate[seedComplex, labelingChoise], "ECGravReturn$55"]
];

(************)

makeLeastEvenlyConnectedComplex[]:=(*construct the complex with the least even connection made of a highly connected island and least connected facets*)
	Module[{facetOrderpairs,
			allowedFacetOrderPair,
			highlyConnectedIsland,
			unusedVertices,
			remainingFacets},
		facetOrderpairs=
			Table[
				With[{nm=Catch[Do[If[Binomial[n,purity]>=q1,Throw[n]],
					{n,purity,purity*q1}]]},{q1,Floor[(nV-nm)/purity]}],
			{q1,facetOrder,1,-1}];
		
		allowedFacetOrderPair=First@(Select[facetOrderpairs,Total[#]==facetOrder&]);
		
		highlyConnectedIsland=
			RandomVertexLabeledPureSimplicialComplex[{purity,allowedFacetOrderPair[[1]],nV-purity*allowedFacetOrderPair[[2]]}];
		
		unusedVertices=Complement[Range[nV],Union@@highlyConnectedIsland];
		
		remainingFacets=Partition[unusedVertices,{purity}];
		Join[highlyConnectedIsland,remainingFacets]
	];

data=<|
	"state"-><|"complex"->Sort[Sort/@seedComplex],"edgeCount"->Length[Union@@(Subsets[#,{2}]&/@(Sort/@seedComplex))],"weight"->1.0,"energy" ->0.0|>,
	"leastEvenlyConnected"-><|"complex"->makeLeastEvenlyConnectedComplex[],"edgeCount"->0,"weight"->1.0,"energy" ->0.0|>,
	"random1"-><|"complex"->RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder,nV}],"edgeCount"->0,"weight"->1.0,"energy" ->0.0|>,
	"random2"-><|"complex"->RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder,nV}],"edgeCount"->0,"weight"->1.0,"energy" ->0.0|>|>;


data[["leastEvenlyConnected","edgeCount"]]=Length[Union@@(Subsets[#,{2}]&/@(Sort/@data[["leastEvenlyConnected","complex"]]))];
data[["random1","edgeCount"]]=Length[Union@@(Subsets[#,{2}]&/@(Sort/@data[["random1","complex"]]))];
data[["random2","edgeCount"]]=Length[Union@@(Subsets[#,{2}]&/@(Sort/@data[["random2","complex"]]))];

Which[

	labelingChoise==1,

	data[["state","weight"]]=PureComplexFacetStabilizerGroupOrder[seedComplex]*(1.0*(facetOrder!)/(nV!));
	data[["leastEvenlyConnected","weight"]]=PureComplexFacetStabilizerGroupOrder[data[["leastEvenlyConnected","complex"]]]*(1.0*(facetOrder!)/(nV!));
	data[["random1","weight"]]=PureComplexFacetStabilizerGroupOrder[data[["random1","complex"]]]*(1.0*(facetOrder!)/(nV!));
	data[["random2","weight"]]=PureComplexFacetStabilizerGroupOrder[data[["random2","complex"]]]*(1.0*(facetOrder!)/(nV!)),

	labelingChoise==2,

	data[["state","weight"]]=1.0*PureComplexAutomorphismGroupOrder[seedComplex]/(nV!);
	data[["leastEvenlyConnected","weight"]]=1.0*PureComplexAutomorphismGroupOrder[data[["leastEvenlyConnected","complex"]]]/(nV!);
	data[["random1","weight"]]=1.0*PureComplexAutomorphismGroupOrder[data[["random1","complex"]]]/(nV!);
	data[["random2","weight"]]=1.0*PureComplexAutomorphismGroupOrder[data[["random2","complex"]]]/(nV!)

];

data[["state","energy"]]=-Log[1.0*data[["state","weight"]]];
data[["leastEvenlyConnected","energy"]]=-Log[1.0*data[["leastEvenlyConnected","weight"]]];
data[["random1","energy"]]=-Log[1.0*data[["random1","weight"]]];
data[["random2","energy"]]=-Log[1.0*data[["random2","weight"]]];


numsweeps=0;
outWinLength =400 ;
(* length of a table to store running energy values to test equilibriation*)inWinLength=20;(* A segment to be averaged over at the beginning and the end of the table*)

Entable=Table[0.0,{m,4},{n,outWinLength}];(*A table to store energy values of the seed, random and empty tracks for tunning calculation of tests of equilibriation*)

(*Main calculation*)
AllEnMat=
(*AllEnMat is for storing all energy values for diagnostic*)
Reap[
	While[numsweeps<maxNumSweeps,
		numsweeps++;

		If[Mod[numsweeps,200]==0,
			PrintTemporary["Running equilibriation at sweep number ",numsweeps, " complex ",data[[Key["state"],Key["complex"]]]]
		];

		Sow[{data[[Key["state"],Key["energy"]]],data[[Key["leastEvenlyConnected"],Key["energy"]]],data[[Key["random1"],Key["energy"]]],data[[Key["random2"],Key["energy"]]]}];
		Do[
			sweepOutput=RandomPureSimplicialComplexMCMCSweep[data[[Key[replicaName],Key["complex"]]],HoldNumberOfVerticesFixed,
				labelingChoise,1];

			data[[Key[replicaName]]]=sweepOutput;

		,{replicaName,{"state","leastEvenlyConnected","random1","random2"}}];

		Entable[[All,1]]={data[[Key["state"],Key["energy"]]],
			data[[Key["leastEvenlyConnected"],Key["energy"]]],data[[Key["random1"],Key["energy"]]],data[[Key["random2"],Key["energy"]]]};

		If[numsweeps>outWinLength,

			sqMeanEMat=Table[(Mean[Entable[[i]][[1;;inWinLength]]]-Mean[Entable[[j]][[(outWinLength-inWinLength);;outWinLength]]])^2,{i,4},{j,4}];

			(*sqMeanEMat is a four by four matrix of the squared means of difference in energy within and across the tracks at the beginning and end of outWinLength. At equilibrium these 9 mumbers should be randomly distributed with mean of 0. Their fluctuation from 0 should be within the variance of the newer(late time) variance in energy for the equilibriation to exit. sqMeanPairwiseDiff is the mean of these 9 numbers*)

			sqMeanPairwiseDiff=Mean[Flatten[sqMeanEMat]];

			meanLateVar=Mean[Table[Variance[Entable[[i]][[1;;inWinLength]]],{i,4}]];(*Mean variance of the newer energies*)

			If[Abs[sqMeanPairwiseDiff]<meanLateVar,
				
				eqlTime=numsweeps-outWinLength+inWinLength;
				Break[]
			];(*Exit the while loop and go to the Do loop*)
			If[Abs[Tr[sqMeanEMat]]<0.000001,	eqlTime=numsweeps-outWinLength+inWinLength;
				Break[]
			](*Exit the while loop and go to the Do loop*)
		];

		Entable=RotateRight[#,1]&/@Entable;(*Shifts every entry to the right by one cyclically. The new data will be written on the first slot so newer data is at the beginning.*)
	];
][[2,1]];

(********************
*  For diagnostics *
********************)

Print[ListLinePlot[Transpose[AllEnMat[[1;;Min[Length[AllEnMat],2*eqlTime]]]],
	PlotRange->All,
PlotLabel->"t from 1 to 2 times eqlT ",PlotLegends->{"seedComplex","leastEvenlyConnected","random1","random2"}]];

(*Print["data",data];*)

result=<| "eqlT"->eqlTime,"state"->data[[Key["state"]]]|>;

result, "ECGravReturn$55"]

];


(* Catch-all Pattern *)
RandomPureSimplicialComplexMCMCEquilibriate[args___]:=(Message[RandomPureSimplicialComplexMCMCEquilibriate::argerr, args];
$Failed);


(* ::Item:: *)
(*RandomPureSimplicialComplexMCMCCorrelationTime*)


(* :Code Section *)

(* Primary Pattern *)
RandomPureSimplicialComplexMCMCCorrelationTime[seedComplex_List,eqlT_Integer, operators_/;MatchQ[operators,{__Function}],
	labelingChoise_Integer]:=
(*
(****************************************)
(*   (* Last updated 05/05/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*takes as an input equilibriated complex, the equlibriation time, a list of functions, and 
	a choice of labeleing and determines the correltion time. Correlation time is the maximum 
	of the correlation times for energy = -Log[prob], the operator which counts the number
	of edges, and the input list of operators. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. *)

*)

Module[{purity=Length[seedComplex[[1]]],facetOrder=Length[seedComplex],
		data,sweepOutput,observablesTable,
		fluctuatingObservableIndices,norm,corrTable,tmaxVals,corrTValues,numsweeps},

data=<|"eqlT"->eqlT,"corrT"->2,"corrTValues"->Table[2,{Length[operators]+2}],
		"state"-><|"complex"->Sort[Sort/@seedComplex],"vertexCount"->Length[DeleteDuplicates[Flatten[seedComplex]]],
		"weight"->1.0,"energy" ->0.0|>|>;

numsweeps=10*eqlT;

PrintTemporary["computing correlation time with numsweeps ",numsweeps];

observablesTable=
Transpose[
	Table[
		If[Mod[i,Ceiling[numsweeps/5.0]]==0,PrintTemporary[" sweepno ",i]];
		sweepOutput=RandomPureSimplicialComplexMCMCSweep[data[[Key["state"],Key["complex"]]],labelingChoise,1];
		data[[Key["state"]]]=sweepOutput;
		Flatten[{data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["vertexCount"]]],
					Through[operators[data[[Key["state"],Key["complex"]]]]]}]
		,{i,1,numsweeps}]
	];

(*If the measured observable values have all equal value
then correlation time can not be computed, so fluctuating selects the measurements whose correlation time can be computed. 
If none of the measurements fluctuate and are all stuck at a single value, the program 
will return the default values corrT -> 2 and corrTValues->{2,2,...,2} *)

fluctuatingObservableIndices=
Table[
	If[Length[DeleteDuplicates[observablesTable[[i]]]]>1,i,
		Which[i==1,

				Message[RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the energy", observablesTable[[1,1]]],

			i==2,
				Message[RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the magnetization", observablesTable[[2,1]]],

			i>2,
				Message[RandomPureSimplicialComplexMCMCCorrelationTime::stuck, operators[[i-2]], observablesTable[[i-2,1]]]];
	Nothing],
{i,1,Length[observablesTable]}];

If[fluctuatingObservableIndices=={},
		Message[RandomPureSimplicialComplexMCMCCorrelationTime::alldefault,
			DeleteDuplicates[#]&/@observablesTable,data[[Key["corrTValues"]]]],
	corrTable=Table[
		norm=CorrelationTime[0,observablesTable[[i]]];
		If[norm==0.0,norm=1.0];
			(*the time 0 correlation can be 0 sometimes*)
		Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,PrintTemporary["      computing corrT at t = ",t]];
			CorrelationTime[t,observablesTable[[i]]],
			{t,0,numsweeps-10}]/norm
		,{i,fluctuatingObservableIndices}];
	
	tmaxVals=Table[(FirstPosition[cT,_?(#<=0&),numsweeps-10])[[1]],{cT,corrTable}]; 
		(* A place to stop the integration for calculation of correlation time is when the autocorrelation value first becomes negative*)
	
	
	corrTValues=Table[Max[Ceiling[Sum[corrTable[[i,t]],{t,1,tmaxVals[[i]]}]],2],
		{i,1,Length[corrTable]}];
		
	data[[Key["corrT"]]]=Max[corrTValues];

	Do[data[[Key["corrTValues"],j]]=
		corrTValues[[First@Flatten[Position[fluctuatingObservableIndices,j]]]],
	{j,fluctuatingObservableIndices}];
];

data

];


(* Overload Pattern *)
RandomPureSimplicialComplexMCMCCorrelationTime[seedComplex_List,HoldNumberOfVerticesFixed_?BooleanQ,eqlT_Integer, operators_/;MatchQ[operators,{__Function}],
	labelingChoise_Integer]:=
(*
(****************************************)
(*   (* Last updated 06/12/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*takes as an input equilibriated complex, the equlibriation time, a list of functions, and 
	a choice of labeleing and determines the correltion time. Correlation time is the maximum 
	of the correlation times for energy = -Log[prob], the operator which counts the number
	of edges, and the input list of operators. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. *)

*)

Module[{purity=Length[seedComplex[[1]]],facetOrder=Length[seedComplex],nV=Length[DeleteDuplicates[Flatten[seedComplex]]],data,sweepOutput,observablesTable,
	fluctuatingObservableIndices,norm,corrTable,tmaxVals,corrTValues,numsweeps},

Catch[If[HoldNumberOfVerticesFixed==False,Throw[RandomPureSimplicialComplexMCMCCorrelationTime[seedComplex,eqlT, operators,labelingChoise], "ECGravReturn$56"]];

data=<|"eqlT"->eqlT,"corrT"->2,"corrTValues"->Table[2,{Length[operators]+2}],
	"state"-><|"complex"->Sort[Sort/@seedComplex],"edgeCount"->Length[Union@@(Subsets[#,{2}]&/@(Sort/@seedComplex))],"weight"->1.0,"energy" ->0.0|>|>;

numsweeps=10*eqlT;

PrintTemporary["computing correlation time with numsweeps ",numsweeps];


observablesTable=
Transpose[
	Table[
		If[Mod[i,Ceiling[numsweeps/5.0]]==0,PrintTemporary[" sweepno ",i]];
		sweepOutput=RandomPureSimplicialComplexMCMCSweep[data[[Key["state"],Key["complex"]]],HoldNumberOfVerticesFixed,labelingChoise,1];
		data[[Key["state"]]]=sweepOutput;
		Flatten[{data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["edgeCount"]]],
					Through[operators[data[[Key["state"],Key["complex"]]]]]}]
	,{i,1,numsweeps}]
];

(*If the measured observable values have all equal value
then correlation time can not be computed, so fluctuating selects the measurements whose correlation time can be computed. If none of the measurements fluctuate and are all stuck at a single value, the program will return the default values corrT -> 2 and corrTValues->{} *)

fluctuatingObservableIndices=
Table[
	If[Length[DeleteDuplicates[observablesTable[[i]]]]>1,i,
		Which[i==1,
			Message[RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the energy", observablesTable[[1,1]]],
			i==2,
			Message[RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the magnetization", observablesTable[[2,1]]],
			i>2,
			Message[RandomPureSimplicialComplexMCMCCorrelationTime::stuck, operators[[i-2]], observablesTable[[i-2,1]]]];
		Nothing],
{i,1,Length[observablesTable]}];

If[fluctuatingObservableIndices=={},
		Message[RandomPureSimplicialComplexMCMCCorrelationTime::alldefault,
			DeleteDuplicates[#]&/@observablesTable,data[[Key["corrTValues"]]]],
	corrTable=
		Table[
			norm=CorrelationTime[0,observablesTable[[i]]];If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)
			Table[
				If[Mod[t,Ceiling[numsweeps/5.0]]==0,PrintTemporary["      computing corrT at t = ",t]];
				CorrelationTime[t,observablesTable[[i]]],
			{t,0,numsweeps-10}]/norm
		,{i,fluctuatingObservableIndices}];
	
	tmaxVals=Table[(FirstPosition[cT,_?(#<=0&),numsweeps-10])[[1]],{cT,corrTable}]; 
	(* A place to stop the integration for calculation of correlation time is when the autocorrelation value first becomes negative*)
	
	
	corrTValues=Table[Max[Ceiling[Sum[corrTable[[i,t]],{t,1,tmaxVals[[i]]}]],2],{i,1,Length[corrTable]}];
		
	data[[Key["corrT"]]]=Max[corrTValues];

	Do[data[[Key["corrTValues"],j]]=corrTValues[[First@Flatten[Position[fluctuatingObservableIndices,j]]]],{j,fluctuatingObservableIndices}];
];


data, "ECGravReturn$56"]

];


(* Catch-all Pattern *)
RandomPureSimplicialComplexMCMCCorrelationTime[args___]:=(Message[RandomPureSimplicialComplexMCMCCorrelationTime::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexMCMC*)


(* :Code Section *)

(* Primary Pattern *)
RandomPureSimplicialComplexMCMC[seedComplex_List, operators_/;MatchQ[operators,{__Function}],
	labelingChoise_Integer, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 05/05/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*takes as an input seed complex, a list of functions, a choice of labeleing and 
	number of samples and outputs a list of numSamples independent random pure simplicial 
	complexes from a uniform disribution. it first runs equilibriation run and determines correlation times
	labelingChoise is an integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. *)

*)

Module[{data,Tempoutput,numsweeps,stopnum,printCase,measurements},

PrintTemporary[" Starting RandomPureSimplicialComplexMCMC "];
(*
(********************)
(*   Equilibriate   *)
(********************)
*)

data=RandomPureSimplicialComplexMCMCEquilibriate[seedComplex,labelingChoise];

(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=RandomPureSimplicialComplexMCMCCorrelationTime[data[["state","complex"]],data[["eqlT"]], operators,
	labelingChoise];

data[["corrT"]]=Tempoutput[[Key["corrT"]]];
data[["corrTValues"]]=Tempoutput[[Key["corrTValues"]]];
data[[Key["state"]]]=Tempoutput[[Key["state"]]];

(*(******************************************
* Take numSamples samples **
******************************************)*)

numsweeps=1;
printCase=Floor[(numSamples*1.0)/5.0];

measurements=Reap[
	While[numsweeps<=numSamples,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		data[[Key["state"]]]=RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],labelingChoise,data[["corrT"]]];

		If[Mod[numsweeps,1]==0,
			Sow[
				Join[
					Flatten[{numsweeps,data[[Key["state"],Key["energy"]]],
						data[[Key["state"],Key["vertexCount"]]],
						Through[operators [data[[Key["state"],Key["complex"]]]]]}],
					{data[[Key["state"],Key["complex"]]]}
				]
			];
		];

		numsweeps++;

	]
][[2,1]];

{measurements,data}

]


(* Overload Pattern *)
RandomPureSimplicialComplexMCMC[seedComplex_List, HoldNumberOfVerticesFixed_?BooleanQ, operators_/;MatchQ[operators,{__Function}],
	labelingChoise_Integer, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 06/12/2026. *) *)
(*   (* Version: 1 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*takes as an input seed complex, a list of functions, a choice of labeleing and 
	number of samples and outputs a list of numSamples independent random pure simplicial 
	complexes from a uniform disribution. it first runs equilibriation run and determines correlation times
	labelingChoise is an integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. *)

*)

Module[{data,Tempoutput,numsweeps,stopnum,printCase,measurements},

Catch[If[HoldNumberOfVerticesFixed==False,Throw[RandomPureSimplicialComplexMCMC[seedComplex, operators,labelingChoise, numSamples], "ECGravReturn$57"]];

PrintTemporary[" Starting RandomPureSimplicialComplexMCMC "];

(*
(********************)
(*   Equilibriate   *)
(********************)
*)


data=RandomPureSimplicialComplexMCMCEquilibriate[seedComplex,HoldNumberOfVerticesFixed,labelingChoise];


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

	
Tempoutput=RandomPureSimplicialComplexMCMCCorrelationTime[data[["state","complex"]],HoldNumberOfVerticesFixed,data[["eqlT"]], operators,
	labelingChoise];

data[["corrT"]]=Tempoutput[[Key["corrT"]]];
data[["corrTValues"]]=Tempoutput[[Key["corrTValues"]]];
data[[Key["state"]]]=Tempoutput[[Key["state"]]];


(*(******************************************
* Take numSamples samples **
******************************************)*)

numsweeps=1;
printCase=Floor[(numSamples*1.0)/5.0];

measurements=
	Reap[
		While[numsweeps<=numSamples,

			If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

			data[[Key["state"]]]=RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],
										HoldNumberOfVerticesFixed,labelingChoise,data[["corrT"]]
								];

			If[Mod[numsweeps,1]==0,
				Sow[
					Join[
						Flatten[{numsweeps,data[[Key["state"],Key["energy"]]],
								data[[Key["state"],Key["edgeCount"]]],
								Through[operators [data[[Key["state"],Key["complex"]]]]]}],
						{data[[Key["state"],Key["complex"]]]}]
				];
			];

			numsweeps++;

		]
	][[2,1]];

{measurements,data}, "ECGravReturn$57"]

];


(* Overload Pattern *)
RandomPureSimplicialComplexMCMC[priorRunInput_Association, operators_/;MatchQ[operators,{__Function}],
	labelingChoise_Integer, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 05/05/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*An overload that takes the input of a previous RandomPureSimplicialComplexMCMC and 
	a new RandomPureSimplicialComplexMCMC. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. It assumes
	the previous run output has computed all the necessary equilibriation and correlation
	time estimations so only runs the final sample generation. *)

*)
Module[{data,Tempoutput,numsweeps,stopnum,printCase,measurements},

data=priorRunInput;

numsweeps=1;
printCase=Floor[(numSamples*1.0)/5.0];

measurements=Reap[
	While[numsweeps<=numSamples,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		data[[Key["state"]]]=
			RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],
				labelingChoise,data[["corrT"]]];

		If[Mod[numsweeps,1]==0,
			Sow[Join[Flatten[{numsweeps,data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["vertexCount"]]],Through[operators [data[[Key["state"],Key["complex"]]]]]}],{data[[Key["state"],Key["complex"]]]}]];
		];

	numsweeps++;

	]
][[2,1]];


{measurements,data}

];


(* Overload Pattern *)
RandomPureSimplicialComplexMCMC[priorRunInput_Association,HoldNumberOfVerticesFixed_?BooleanQ, operators_/;MatchQ[operators,{__Function}],
	labelingChoise_Integer, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 06/12/2026. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: .*) *) 
(****************************************)
(*An overload that takes the input of a previous RandomPureSimplicialComplexMCMC and 
	a new RandomPureSimplicialComplexMCMC. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. It assumes
	the previous run output has computed all the necessary equilibriation and correlation
	time estimations so only runs the final sample generation. *)

*)
Module[{data,Tempoutput,numsweeps,stopnum,printCase,measurements},


Catch[If[HoldNumberOfVerticesFixed==False,Throw[RandomPureSimplicialComplexMCMC[priorRunInput, operators,labelingChoise, numSamples], "ECGravReturn$58"]];

data=priorRunInput;

numsweeps=1;
printCase=Floor[(numSamples*1.0)/5.0];

measurements=Reap[
	While[numsweeps<=numSamples,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		data[[Key["state"]]]=RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],
										HoldNumberOfVerticesFixed,labelingChoise,data[["corrT"]]
							];


		If[Mod[numsweeps,1]==0,
			Sow[Join[Flatten[{numsweeps,data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["edgeCount"]]],Through[operators [data[[Key["state"],Key["complex"]]]]]}],{data[[Key["state"],Key["complex"]]]}]];
		];

		numsweeps++;

	]
][[2,1]];

{measurements,data}, "ECGravReturn$58"]

];


(* ::Chapter:: *)
(*Generate A Random Pseudo Manifold*)


(* ::Subsection::Closed:: *)
(*Random Pseudo Manifold Through Successive Facet Addition*)


(* ::Item::Closed:: *)
(*RandomPseudoManifold*)


(* :Code Section *)

(* Primary Pattern *)
RandomUnlabeledPseudoManifold[{purity_Integer,facetOrder_Integer,connected_Integer:0}]:=
(*
(****************************************)
(*   (* Last updated 12/28/2024. *) *)
(*   (* Note:  .*) *) 
(****************************************)
(*A code to generate a random pseudomanifold of a given purity and facet order 
through iterative addition of vertices labeled by the facets. Input is the purity, 
the facet order, and integer "connected" that is 1 for connected pseudomanifolds and 
0 for disconnected.  E.g. RandomPseudoManifold[2,4,1] generates a 2-pure, connected 
random pseudomanifold with 4 edges. It successively choses a random vertex from a 
decreasing list of available vertex space.  *)

*)
Module[{facets={Range[purity]},pstingSites=Subsets[Range[purity],{purity-1}]},

Catch[If[purity<2,
If[facetOrder==1,Throw[facets, "ECGravReturn$59"]];
Throw[Null, "ECGravReturn$59"]
];

If[facetOrder==0,
Throw[Null, "ECGravReturn$59"]];

If[facetOrder==1,Throw[facets, "ECGravReturn$59"]
];


Do[

{facets,pstingSites} =AddRandomUnlabeledFacetToPseudoManifold[facets,pstingSites,connected],{n,1,facetOrder-1}

];

{facets,pstingSites}, "ECGravReturn$59"]

];

(* Catch-all Pattern *)
RandomUnlabeledPseudoManifold[args___]:=(Message[RandomUnlabeledPseudoManifold::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Add  a random unlabeled facet to a pseudo-manifold*)


(* ::Item::Closed:: *)
(*AddRandomUnlabeledFacetToPseudoManifold*)


(* :Code Section *)

(* Primary Pattern *)
AddRandomUnlabeledFacetToPseudoManifold[facetsLst_List,apastingSites_List,connected_Integer:0]:=
(*
(****************************************)
(*   (* Last updated 10/11/2025. *) *)
(*   (* Note:  .*) *) 
(****************************************)
(*A program that adds one more facet to existing labeled pseudo manifold at random. 
Input is the current pseudo manifold as a list of facets, a second list of  codimension 
1 faces available for binding, and an integer, connected, taking the values 0 or 1; 0 
is for disconnected pseudomanifolds is the default. 1 is for connected pseudomanifolds. *)

*)
Module[{purity=Length[facetsLst[[1]]],glen,g,relabelingRule,autGroup,facets,newfacets,
	newPastingSites,newclqVset,IsPartOfAFacet,CheckMatchingFunctionness,
	CheckVertexDistance,bindingSiteSpace,bindingSite,AllowedPartnerSitesQ,
	FindAllowedPastingSites,PasteSites,randomOrderOfPastingSites,allAllowedPSites,
	success,relabelingRules},

Catch[If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Throw[Null, "ECGravReturn$60"]];If[DeleteDuplicates[Length/@facetsLst][[1]]!=purity,
Throw[Null, "ECGravReturn$60"]];


glen=Length[DeleteDuplicates[Flatten[facetsLst]]];

(*facetAutGroup=PureComplexFacetAutomorphismGroup[facets];*)

{autGroup,relabelingRule}={#1,Reverse/@#2}&@@PureComplexAutomorphismGroup[facetsLst];


g=CanonicalGraph[GraphFromCliques[facetsLst]];

(*relabelingRule=Normal[FindGraphIsomorphism[GraphFromCliques[facetsLst],g][[1]]];*)


facets=Sort[Sort/@(facetsLst/.relabelingRule)];

newPastingSites=Sort[Sort/@(apastingSites/.relabelingRule)];


newclqVset=Range[glen+1,glen+purity];(*vertex set for the new clique to be added*)


If[apastingSites=={},Throw[{Join[facets,{newclqVset}],Subsets[newclqVset,{purity-1}]}, "ECGravReturn$60"]];
(*If apastingSites=={}, the manifold has no boundaries, so the program will return the 
manifold and the new disconnected facet*)

(*bindingSiteSpace=Subsets[Subsets[newclqVset,{purity-1}],purity];*)

bindingSiteSpace=First/@Table[Subsets[Subsets[newclqVset,{purity-1}],{i}],{i,0,purity}];


bindingSiteSpace=bindingSiteSpace/.Thread[Range[purity]->newclqVset];


If[connected==1,bindingSiteSpace=Complement[bindingSiteSpace,{{}}]];


IsPartOfAFacet[a_List]:=(*Given a list of codimension 1 faces, this checks whether any subset of a pair is contained within a single facet in the pseudo manifold. Returns true if there is a facet containing any two of the faces in the list and false otherwise*)
With[{pairs=Apply[Union,Subsets[a,{2}],1]},
(*Or@@(Union@@Table[SubsetQ[i,#]&/@pairs,{i,facets}])*)
AnyTrue[Table[Table[{i,j},{i,facets}],{j,pairs}],SubsetQ[#[[1]],#[[2]]]&,2]
];

(*Given a matching, decide whether or not it is a function, i.e., every member of the range has a unique preimage. It does so by finding the image of each element of the domain, calculating the size of such sets, and checking to see if any have cardinality \[GreaterEqual]2. Returns true if the matching is a proper function. It will be used to determine whether the pasting from partner sites in the graph to the binding sites in the new clique is a function.*)
CheckMatchingFunctionness[ps_List,bs_List]:=Block[{asn},
asn=<||>;
Table[asn[[Key[ps[[i]]]]]=Union[Lookup[asn,Key[ps[[i]]],{}],{bs[[i]]}],{i,1,Length[ps]}];
SelectFirst[Length/@asn,#>=2&,-1]==-1
];

(* Check that for a given set of vertices in the graph or psite that are mapped to the same vertex in bsite (the new clique) the distance in the original graph between any pair of them is greater than 2 so that they don't form multiloops or self loops. Returns true if no such case is found so that the pasting doesn't create any accidental self loops or multiloops. *)
CheckVertexDistance[ps_List,bs_List]:=Block[{asn,pairs,distances},
asn=<||>;
(*Table[asn[[Key[bs[[i]]]]]=Append[Lookup[asn,Key[bs[[i]]],{}],ps[[i]]],{i,1,Length[bs]}];*)
Table[asn[[Key[bs[[i]]]]]=Union[Lookup[asn,Key[bs[[i]]],{}],{ps[[i]]}],{i,1,Length[bs]}];

(*Print["            In CheckVertexDistance, ps ",ps, " bs ",bs];
Print["            In CheckVertexDistance, asn ",asn];*)

pairs=Subsets[#,{2}]&/@asn;


distances=Flatten[Table[GraphDistance[g,#[[1]],#[[2]]]&/@i,{i,pairs}]];


SelectFirst[distances,#<=2&,-1]==-1 (*If there are no pairs with distance of 2 or less, SelectFirst ourputs -1 and the check returns true*)
];


AllowedPartnerSitesQ[candPsite_List]:=(*Given candPsite, a candidate partnersite list of codim 1 faces in the current complex( which will end up being part of the new facet that will be added) it checks whether or not a "slit" will be formed after pasting the new facet. A slit is basically two codim 1 faces that end up being identified at all their vertices. Returns True if the candPsite is ok, and False otherwise. *)
Block[{vset=Union@@candPsite,facetComplements},
Catch[If[Length[candPsite]<2,Throw[True, "ECGravReturn$61"]];
facetComplements=Complement[Subsets[vset,{purity-1}],candPsite];


NoneTrue[Tuples[{facets,facetComplements}],SubsetQ[#[[1]],#[[2]]]&], "ECGravReturn$61"]
];

FindAllowedPastingSites[ ]:=(*Finds the list of all lists of codimension 1 faces in the current allowed by the pasting rules sites in the complex.*)Module[{allPastingSites,validPastingSites,validPastingSiteOrbits,allUniquePastingSites},


allPastingSites=Subsets[newPastingSites,{1,purity}];


allPastingSites=Select[allPastingSites,!IsPartOfAFacet[#]&];


validPastingSites=Select[allPastingSites,AllowedPartnerSitesQ[#]&];


validPastingSiteOrbits=GroupOrbits[autGroup,validPastingSites];


(*Sorting is needed because Mathematica treats {1,2} as different from {2,1}*)
validPastingSiteOrbits=Map[Sort,validPastingSiteOrbits,3];


validPastingSiteOrbits=DeleteDuplicates[(DeleteDuplicates/@validPastingSiteOrbits)];

(*removing duplicates after sorting.*)


allUniquePastingSites = validPastingSiteOrbits[[All,1]];
(*Picks a representative from each orbit of the automorphism group of the graph*)


RandomSample[allUniquePastingSites]

];


PasteSites[bsite_List,psite_List]:=(*Given a list of binding site p-1)-faces (bsite) in the new facet and partnering pasting sites psite in the complex, this method pasts the two. It does so by creating a rule which relabels all vertices in the graph with their corresponding label in bsite. This pasting is a surjective (many to one) function from the vertex list of the partnering sites in psite to the vertex set of the binding sites in the clique. The method has side effect, it modifies newfacets. It returns 1 if succesful and 0 otherwise.*)
Block[{psiteApexSet,bsiteApexSet,psiteSet,bsiteSet,rules,tgraph},

Catch[newfacets=facets;


(*care must be taken to make sure vertices that are in the intersection of two or more cliques in psite are also mapped to vertices that are intersections in bsite. psiteApexSet and bsiteApexSet collect such vertices. Then these vertices are mapped to eachother.*)

psiteApexSet=Union@@(Apply[Intersection,Subsets[psite,{2}],1]);


bsiteApexSet=Union@@(Apply[Intersection,Subsets[bsite,{2}],1]);


(*psiteSet=Flatten[psite];
bsiteSet=Flatten[bsite];

Print["      before sorting, psiteSet ",psiteSet, " bsiteSet ",bsiteSet];*)

psiteSet=Flatten[Table[SortBy[i,!MemberQ[psiteApexSet,#]&],{i,psite}]];
bsiteSet=Flatten[Table[SortBy[i,!MemberQ[bsiteApexSet,#]&],{i,bsite}]];


rules=DeleteDuplicates[Thread[psiteSet->bsiteSet]];

(*(** Check 1. Check that the matching from psite to bsite is a function, i.e., each vertex in psite is mapped to atmost one vertex in bsite **)*)

If[CheckMatchingFunctionness[psiteSet,bsiteSet]==False,(*Print["Can not paste. No function from psite to bsite. Returning 0. "];*)Throw[0, "ECGravReturn$62"]
];

(*(** Check 2.  For each member of bsiteApexSet, check that its preimage, i.e, all vertices in psite mapped to the same vertex in bsiteApexSet have distance of greater than 2. **)*)

If[CheckVertexDistance[psiteSet,bsiteSet]==False,(*Print["Two or more vertices of separations less than 3 will be identified by this pasting. Returning 0 "];*)
Throw[0, "ECGravReturn$62"]
];


newfacets=Sort[Sort/@(newfacets/.rules)];


AppendTo[newfacets,newclqVset];


newPastingSites=Complement[newPastingSites,psite];


newPastingSites=Sort[Sort/@(newPastingSites/.rules)];


newPastingSites=Join[newPastingSites,Complement[Subsets[newclqVset,{purity-1}],bsite]];


newPastingSites=Sort[Sort/@newPastingSites];


(*tgraph=GraphFromCliques[newfacets];*)


1, "ECGravReturn$62"]

];


(*****************************
*        main loop          *
****************************
*)

(*randomOrderOfPastingSites=RandomSample[Select[bindingSiteSpace,Length[#]\[LessEqual]Length[facets]&]];*)


(*Print[""];
Print[""];
Print[""];
Print[""];
Print["Main loop "];
(*Print[ "bg ",bg, " Length[facets] ",Length[facets]];*)Print["facets ",facets," graph, ",GraphPlot[GraphFromCliques[facets],VertexLabels->Automatic,ImageSize->{65,65}]];
Print["relabeled pastingSites ",newPastingSites];
Print[""];
Print[""];*)


(******************************************************
*     Select pasting sites from                      *
*     the existing graph                             *
*                                                    *
******************************************************)

allAllowedPSites=FindAllowedPastingSites[];


(********************************
*   Commence pasting           *
*******************************)

Do[


success=0;

bindingSite=SelectFirst[bindingSiteSpace,Length[#]==Length[randomPsite]&];


success = PasteSites[SortBy[bindingSite,Length],SortBy[randomPsite,Length]];


If[success==1,


Break[],Continue[]
];


,{randomPsite,allAllowedPSites}] ;

(*Print[" before relabel, newfacets ",newfacets, " relabel rules ",
	Thread[DeleteDuplicates[Flatten[newfacets]]\[Rule]Range[Length[DeleteDuplicates[
			Flatten[newfacets]]]]] ];
	Print[ " after relabel, newfacets ",newfacets/.Thread[DeleteDuplicates[Flatten[
		newfacets]]\[Rule]Range[Length[DeleteDuplicates[Flatten[newfacets]]]]]];*)

(*Print[" before relabeling, manifold is ",newfacets, " as graph ",
	GraphPlot[GraphFromCliques[newfacets],VertexLabels->Automatic,ImageSize->{65,65}]];*)

relabelingRules=Thread[DeleteDuplicates[Flatten[newfacets]]
	->Range[Length[DeleteDuplicates[Flatten[newfacets]]]]];

(*Print[" After relabeling, manifold is ",newfacets/.relabelingRules, 
	" as graph ",GraphPlot[GraphFromCliques[newfacets/.relabelingRules],
	VertexLabels->Automatic,ImageSize->{65,65}]];*)

(*newfacets/.Thread[DeleteDuplicates[Flatten[newfacets]]\[Rule]Range[Length[DeleteDuplicates[Flatten[newfacets]]]]]*)

(*Print["AddRandomFacetToPseudoManifold returning new facets ",newfacets," as graphs ",
	GraphPlot[GraphFromCliques[newfacets],VertexLabels->Automatic,ImageSize->{70,70}], 
	" newPastingSites ",newPastingSites];*)

(*If[Length[Union@(Sort/@newfacets)]<=Length[facetsLst],Print["WARNING!! WARNING!! 
	current complex ",facetsLst," new complex ",newfacets ]];*)

{Sort[Sort/@(newfacets/.relabelingRules)], Sort[Sort/@(newPastingSites/.relabelingRules)]}, "ECGravReturn$60"]

];

(* Catch-all Pattern *)
AddRandomUnlabeledFacetToPseudoManifold[args___]:=(Message[AddRandomUnlabeledFacetToPseudoManifold::argerr, args];
$Failed);


(* ::Title:: *)
(*PureComplexes Protect and End*)


(* End private context *)

(* Protect exported symbols *)
