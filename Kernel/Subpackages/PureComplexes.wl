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


BeginPackage["ECGrav`PureComplexes`"];


(* ::Title:: *)
(*PureComplexes Public Functions*)


(* ::Subtitle:: *)
(*Main PureComplexes Public Symbols*)


(* ::Title:: *)
(*PureComplexes Private*)


Begin["`Private`"] (* Begin private context *)


(* ::Chapter:: *)
(*Helper Functions*)


(* ::Section::Closed:: *)
(*General*)


(* ::Subsection::Closed:: *)
(*Basic Simplicial Complex Constructions*)


(* ::Item::Closed:: *)
(*FacetIncidenceMatrix*)


(* Primary Pattern *)
ECGrav`FacetIncidenceMatrix[facetsLst_List]:=
With[{facetOrder=Length[facetsLst],vlist=DeleteDuplicates[Flatten[facetsLst]]},
	Table[Table[If[MemberQ[i,j],1,0],{j,vlist}],{i,facetsLst}]
];

(* Overload Pattern *)
ECGrav`FacetIncidenceMatrix[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},
	ECGrav`FacetIncidenceMatrix[clqs]
];

(* Catch-all Pattern *)
ECGrav`FacetIncidenceMatrix[args___]:=(Message[ECGrav`FacetIncidenceMatrix::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetAdjacencyMatrix*)


(* Primary Pattern *)
ECGrav`FacetAdjacencyMatrix[facetsLst_List]:=Module[{clqOrder=Length[facetsLst],facetAdjMat=DiagonalMatrix[Length/@facetsLst]},
Do[
Do[
facetAdjMat[[i,j]]=facetAdjMat[[j,i]]=Length[Intersection[facetsLst[[i]],facetsLst[[j]]]]
,{j,i+1,clqOrder}]
,{i,1,clqOrder-1}];
facetAdjMat
];

(* Overload Pattern *)
ECGrav`FacetAdjacencyMatrix[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},
ECGrav`FacetAdjacencyMatrix[clqs]
];

(* Catch-all Pattern *)
ECGrav`FacetAdjacencyMatrix[args___]:=(Message[ECGrav`FacetAdjacencyMatrix::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*GraphFromCliques*)


(* Primary Pattern *)
ECGrav`GraphFromCliques[clqs_List]:=With[{vertexlist=DeleteDuplicates[Flatten[clqs]],edgelist=UndirectedEdge@@@(DeleteDuplicates[Flatten[Subsets[Sort[#],{2}]&/@clqs,1]])},
Graph[vertexlist,edgelist]
];

(* Catch-all Pattern *)
ECGrav`GraphFromCliques[args___]:=(Message[ECGrav`GraphFromCliques::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*GraphFromFacetIncidence*)


(* Primary Pattern *)
ECGrav`GraphFromFacetIncidence[incidenceMat_List/;MatrixQ[incidenceMat]]:=
With[{clqsLst=Table[Pick[Range[Length[incidenceMat[[1]]]],i,1],{i,incidenceMat}]},
ECGrav`GraphFromCliques[clqsLst]
];

(* Catch-all Pattern *)
ECGrav`GraphFromFacetIncidence[args___]:=(Message[ECGrav`GraphFromFacetIncidence::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliquesFromFacetIncidence*)


(* Primary Pattern *)
ECGrav`CliquesFromFacetIncidence[incidenceMat_List/;MatrixQ[incidenceMat]]:=
With[{clqsLst=Table[Pick[Range[Length[incidenceMat[[1]]]],i,1],{i,incidenceMat}]},
clqsLst
];

(* Catch-all Pattern *)
ECGrav`CliquesFromFacetIncidence[args___]:=(Message[ECGrav`CliquesFromFacetIncidence::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ComplexFromFacetLabeledVertexList*)


(* Primary Pattern *)
ECGrav`ComplexFromFacetLabeledVertexList[facetLabeledVertices_List]:=
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
ECGrav`ComplexFromFacetLabeledVertexList[args___]:=(Message[ECGrav`ComplexFromFacetLabeledVertexList::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetLabeledVertexListFromComplex*)


(* Primary Pattern *)
ECGrav`FacetLabeledVertexListFromComplex[facetsLst_List]:=
With[{vlist=DeleteDuplicates[Flatten[facetsLst]],
cmlxAsn=<|Table[i->facetsLst[[i]],{i,1,Length[facetsLst]}]|>},
(*ReverseSort[Table[Keys[Select[cmlxAsn,MemberQ[#,v]&]],{v,vlist}]]*)
Join@@(Map[LexicographicSort[#]&,
GatherBy[ReverseSortBy[Table[Keys[Select[cmlxAsn,MemberQ[#,v]&]],{v,vlist}],Length],Length],{1}])
];

(* Catch-all Pattern *)
ECGrav`FacetLabeledVertexListFromComplex[args___]:=(Message[ECGrav`FacetLabeledVertexListFromComplex::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexQ*)


(* Primary Pattern *)
ECGrav`PureComplexQ[facets_List]:=If[Length[DeleteDuplicates[Length/@facets]]==1,Return[True],Return[False]];

(* Catch-all Pattern *)
ECGrav`PureComplexQ[args___]:=(Message[ECGrav`PureComplexQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureGraphQ*)


(* Primary Pattern *)
ECGrav`PureGraphQ[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},
If[Length[DeleteDuplicates[Map[Length,clqs]]]==1,Return[True],Return[False]]
];

(* Catch-all Pattern *)
ECGrav`PureGraphQ[args___]:=(Message[ECGrav`PureGraphQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueComplexQ*)


(* Primary Pattern *)
ECGrav`CliqueComplexQ[facetsLst_List]:=
(*Given a list of maximal cliques, it checks whether or not it passes the clique condition test*)
Module[{threeSubsets=Subsets[facetsLst,{3}],fIsInAClique},
If[Length[facetsLst]<3,Return[True]];
fIsInAClique[candClq_List]:=AnyTrue[facetsLst,SubsetQ[#,candClq]&];
AllTrue[threeSubsets,fIsInAClique[Union[Intersection[#[[1]],#[[2]]],Intersection[#[[1]],#[[3]]],Intersection[#[[2]],#[[3]]]]]&]
];

(* Catch-all Pattern *)
ECGrav`CliqueComplexQ[args___]:=(Message[ECGrav`CliqueComplexQ::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Spheres and Balls, Links and Stars*)


(* ::Item::Closed:: *)
(*Sph*)


(* Primary Pattern *)
ECGrav`Sph[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer]:=
With[{size = Length[Amat],rowcolsToKeep=Flatten[Position[Amat[[i]],1]]},
	If[i>size,Message[ECGrav`Sph::vtxnotfound, i];
	Return[{}]];
	If[rowcolsToKeep=={},Return[{}]];
	If[size==0,Return[{}]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep]
];

(* Overload Pattern *)
ECGrav`Sph[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer,r_Integer]:=
Module[{rowcolsToKeep,dm},
	If[Length[Amat]==0,Return[{}]];
	dm=GraphDistanceMatrix[AdjacencyGraph[Amat]];
	rowcolsToKeep=Flatten[Position[dm[[i]],r]];
	If[rowcolsToKeep=={},Return[{}]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep]
];

(* Overload Pattern *)
ECGrav`Sph[g_Graph,i_Integer]:=
Module[{},
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Message[ECGrav`Sph::vtxnotfound, i];{}
	];

	Subgraph[g,AdjacencyList[g,i]]
];

(* Overload Pattern *)
ECGrav`Sph[g_Graph,i_Integer,r_Integer]:=
With[{sphVertices=Select[VertexList[g],GraphDistance[g,i,#]==r&]},
	
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Message[ECGrav`Sph::vtxnotfound, i];{}
	];

	Subgraph[g,sphVertices]
];

(* Catch-all Pattern *)
ECGrav`Sph[args___]:=(Message[ECGrav`Sph::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Bll*)


(* Primary Pattern *)

ECGrav`Bll[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer]:=
With[{size = Length[Amat],rowcolsToKeep=Flatten[Position[Amat[[i]],1]]},
	If[i>size,Message[ECGrav`Bll::vtxnotfound, i];
	Return[{}]];
	If[size==0,Return[{}]];
	If[rowcolsToKeep=={},Return[{}]];
	Join[
		Transpose[
			Join[
				Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep],
				Table[1,{Length[rowcolsToKeep]}]
			]
		],
		PadRight[Table[1,{Length[rowcolsToKeep]}],Length[rowcolsToKeep]+1]
	]
];

(* Overload Pattern *)
ECGrav`Bll[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer,r_Integer]:=
Module[{rowcolsToKeep,dm},
	If[Length[Amat]==0,Return[{}]];
	dm=GraphDistanceMatrix[AdjacencyGraph[Amat]];
	rowcolsToKeep=Flatten[Position[dm[[1]],_?(#<=r&)]];
	If[rowcolsToKeep=={},Return[{}]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep]
];

(* Overload Pattern *)
ECGrav`Bll[g_Graph,i_Integer]:=
(*Unit ball in a graph g at vertex i*)
With[{unitBallVertices=Union[{i},AdjacencyList[g,i]]},
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Message[ECGrav`Bll::vtxnotfound, i];{}
	];

Subgraph[g,unitBallVertices]

];

(* Overload Pattern *)
ECGrav`Bll[g_Graph,i_Integer,r_Integer]:=
(*Ball of radius r in a graph g at vertex i*)
With[{ballVertices=Select[VertexList[g],GraphDistance[g,i,#]<=r&]},
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Message[ECGrav`Bll::vtxnotfound, i];{}
	];

	Subgraph[g,Union[{i},ballVertices]]

];

(* Catch-all Pattern *)
ECGrav`Bll[args___]:=(Message[ECGrav`Bll::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Lnk*)


ECGrav`Lnk[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1}),q_List/;Depth[q]==2]:=
(*Returns the link of the face q in the complex given by the list of facets facetsLst*)
With[{link=Select[facetsLst,SubsetQ[#,q]&],rules=Table[i->Nothing,{i,q}]},
	link/.rules
];

(* Catch-all Pattern *)
ECGrav`Lnk[args___]:=(Message[ECGrav`Lnk::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Str*)


ECGrav`Str[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1}),q_List/;Depth[q]==2]:=
(*Returns the link of the face q in the complex given by the list of facets facetsLst*)
Select[facetsLst,SubsetQ[#,q]&];

(* Catch-all Pattern *)
ECGrav`Str[args___]:=(Message[ECGrav`Str::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Basic Graph Observables*)


(* ::Item::Closed:: *)
(*Deg*)


(* Primary Pattern *)
ECGrav`Deg[Amat_List,i_Integer]:=(*Given an adjacency matrix of a graph Amat, and a vertex i, it gives the degree of the vertex*)
Total[Amat[[i]]];

(* Catch-all Pattern *)
ECGrav`Deg[args___]:=(Message[ECGrav`Deg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*AvgDeg*)


(* Primary Pattern *)
ECGrav`AvgDeg[Amat_List/;!GraphQ[Amat]]:=
(*The average degree of a graph given as an adjacency matrix Amat*)
(1/Length[Amat])*Total[Amat,2];

(* Overload Pattern *)
ECGrav`AvgDeg[g_Graph]:=
(*The average degree of a graph given as a graph object g*)
(1/VertexCount[g])*Total[VertexDegree[g]];

(* Catch-all Pattern *)
ECGrav`AvgDeg[args___]:=(Message[ECGrav`AvgDeg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetDeg*)


(* Primary Pattern *)
ECGrav`FacetDeg[facetsLst_List,v_Integer]:=
(*Given a simplicial complex as a list facetsLst of facets and a vertex v, 
it computes the facet degree (number of facets containing the vertex) *)
Length[Select[facetsLst,MemberQ[#,v]&]];

(* Overload Pattern *)
ECGrav`FacetDeg[facetsLst_List]:=
(*Given a simplicial complex, it lists the facet degrees of all vertices  *)
With[{vertices=Sort[DeleteDuplicates[Flatten[facetsLst]]]},
ECGrav`FacetDeg[facetsLst,#]&/@vertices
];

(* Catch-all Pattern *)
ECGrav`FacetDeg[args___]:=(Message[ECGrav`FacetDeg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HyperDeg*)


(* Primary Pattern *)
ECGrav`HyperDeg[g_Graph,clq_List]:=
(*Computes the degree of the clique clq given as a list of vertices. 
Checks whether or not the input is a clique in the graph*)
With[{spheres=Table[AdjacencyList[g,v],{v,clq}]},
If[CompleteGraphQ[Subgraph[g,clq]]==False,Return[]];
Length[Intersection@@spheres]
];

(* Overload Pattern *)
ECGrav`HyperDeg[Amat_List/;(SymmetricMatrixQ[Amat]&&Sort[DeleteDuplicates[Flatten[Amat]]]=={0,1}),clq_List]:=
(*Computes the degree of the clique clq given as a list of vertices. 
Checks whether or not the input is a clique in the graph*)
With[{g=AdjacencyGraph[Amat]},
If[CompleteGraphQ[Subgraph[g,clq]]==False,Return[]];
ECGrav`HyperDeg[g,clq]
];

(* Overload Pattern *)
ECGrav`HyperDeg[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&
	Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1}),clq_List]:=
(*Computes the hyperdegree of the face clq given as a list of vertices. 
Checks whether or not the input is a face in the graph*)
With[{lnk=ECGrav`Lnk[facetsLst,clq]},
	If[NoneTrue[facetsLst,SubsetQ[#,clq]&],Return[]];
	Total[Length/@lnk]
];

(* Catch-all Pattern *)
ECGrav`HyperDeg[args___]:=(Message[ECGrav`HyperDeg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FVector*)


(* Primary Pattern *)
ECGrav`FVector[g_Graph]:=
(*Computes the f-vector of the graph g which is the vector {f0,f1,...} of 
the number of vertices, number of edges, etc. *)
With[{n0=VertexCount[g],n1=EdgeCount[g],maxClqs=FindClique[g,\[Infinity],All]},
	If[n1==0,Return[{n0}]];
	Join[{n0,n1},
		Table[Length[Union@@(Subsets[#,{q}]&/@maxClqs)],{q,3,Max[Length/@maxClqs]}]
	]
];

(* Overload Pattern *)
ECGrav`FVector[Amat_List/;(SymmetricMatrixQ[Amat]&&Sort[DeleteDuplicates[Flatten[Amat]]]=={0,1})]:=
(*Computes the f-vector of the graph given by its adjacency matrix Amat which is the 
vector {f0,f1,...} of the number of vertices, number of edges, etc. *)
With[{n0=Length[Amat],n1=Total[Amat,2]/2,maxClqs=FindClique[AdjacencyGraph[Amat],\[Infinity],All]},
If[n1==0,Return[{n0}]];
Join[{n0,n1},Table[Length[Union@@(Subsets[#,{q}]&/@maxClqs)],{q,3,Max[Length/@maxClqs]}]]
];

(* Overload Pattern *)
ECGrav`FVector[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&
	Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1})]:=
(*Computes the hyperdegree of the face clq given as a list of vertices. 
Checks whether or not the input is a face in the graph*)
With[{n0=Length[DeleteDuplicates[Flatten[facetsLst]]],maxFacetSize=Max[Length/@facetsLst]},
If[n0==0,Return[{}]];
Join[{n0},Table[Length[Union@@(Subsets[#,{q}]&/@facetsLst)],{q,2,maxFacetSize}]]
];

(* Catch-all Pattern *)
ECGrav`FVector[args___]:=(Message[ECGrav`FVector::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*GVolume*)


(* Primary Pattern *)
ECGrav`GVolume[Amat_List,Dmat_List,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, where the distance 
between nodes is given by Dmat, centered at node i and having radius r*)
With[{n=Length[Amat],inball=Table[If[Dmat[[i,j]]<=r,1,0],{j,n}]},
Total[inball]-1
];

(* Overload Pattern *)
ECGrav`GVolume[g_Graph,Dmat_List,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, where the distance 
between nodes is given by Dmat, centered at node i and having radius r*)
With[{n=VertexCount[g],inball=Table[If[Dmat[[i,j]]<=r,1,0],{j,n}]},
Total[inball]-1
];

(* Overload Pattern *)
ECGrav`GVolume[Amat_List,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, centered at node i 
and having radius r*)
With[{n=Length[Amat],dmat=GraphDistanceMatrix[AdjacencyGraph[Amat]]},
Total[Table[If[dmat[[i,j]]<=r,1,0],{j,n}]]-1
];

(* Overload Pattern *)
ECGrav`GVolume[g_Graph,r_Integer,i_Integer]:=
(*Volume of the ball in the graph with adjacency matrix Amat, centered at node i 
and having radius r*)
With[{n=VertexCount[g],dmat=GraphDistanceMatrix[g]},
Total[Table[If[dmat[[i,j]]<=r,1,0],{j,n}]]-1
];

(* Catch-all Pattern *)
ECGrav`GVolume[args___]:=(Message[ECGrav`GVolume::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*KFaceDistance*)


(* Primary Pattern *)
ECGrav`KFaceDistance[facets_List/;Depth[facets]==3,clq1_List/;Depth[clq1]==2,clq2_List/;Depth[clq2]==2]:=

(*Given a list of maximal cliques (facets) of a clique complex or a graph and two equal
 cardinality faces, clq1 and clq2 of cardinality d in the graph,  
 this method computes the minimum clique distance between them where 
 the path steps on only (d+1)-dimensional faces joined along d-dimensional subfaces.  
 It  uses the built in Mathematica GraphDistance function on a graph constructed 
 such that the vertices are all d-faces and edges between them are 1 if two faces 
 are contained in a bigger clique and 0 otherwise. *)

Module[{d=Length[clq1],dests,pos1,pos2,clqGraphAmat},
If[Complement[clq1,clq2]=={},Return[0]];
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

GraphDistance[AdjacencyGraph[clqGraphAmat],pos1,pos2]

];

(* Catch-all Pattern *)
ECGrav`KFaceDistance[args___]:=(Message[ECGrav`KFaceDistance::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*KFaceDistanceMatrix*)


(* Primary Pattern *)
ECGrav`KFaceDistanceMatrix[facets_List, facedim_Integer]:=
(*Given a list of maximal cliques (facets) of a clique complex or a graph and 
dimensionality clqdim, this method computes the clique distance matrix between 
every pair of faces of cardinality clqdim, where each path steps on only 
(facedim+1)-dimensional faces.  It  uses the built in Mathematica 
GraphDistanceMatrix function on a graph constructed such that the vertices 
are all clqdim-faces and edges between them are 1 if two faces are contained in a 
bigger clique and 0 otherwise. *)

Module[{dests,clqGraphAmat},

If[facedim==1,Return[GraphDistanceMatrix[ECGrav`GraphFromCliques[facets]]]];

dests=Union@@(Subsets[#,{facedim}]&/@facets);

clqGraphAmat = 
	Table[
		Table[
			If[Complement[i,j]=={},0,
				If[(Length[Union[i,j]]==facedim+1)&&AnyTrue[facets,SubsetQ[#,Union[i,j]]&],1,0]
			],{i,dests}
		],
{j,dests}];

GraphDistanceMatrix[AdjacencyGraph[clqGraphAmat]]

];

(* Catch-all Pattern *)
ECGrav`KFaceDistanceMatrix[args___]:=(Message[ECGrav`KFaceDistanceMatrix::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*KpathConnectedComponents*)


(* Primary Pattern *)
ECGrav`KpathConnectedComponents[facets_List,k_Integer]:=

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

If[k==1,Return[ConnectedComponents[ECGrav`GraphFromCliques[facets]]]];

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
	]

];

(* Overload Pattern *)
ECGrav`KpathConnectedComponents[g_Graph,k_Integer]:=
With[{clqs=FindClique[g,\[Infinity],All]},
	ECGrav`KpathConnectedComponents[clqs,k]
];

(* Catch-all Pattern *)
ECGrav`KpathConnectedComponents[args___]:=(Message[ECGrav`KpathConnectedComponents::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ConnectedComplexComponents*)


(* Primary Pattern *)
ECGrav`ConnectedComplexComponents[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it gives a list of simplicial
 complexes which are connected. 
 E.g.ECGrav`ConnectedComplexComponents[{{1,2},{3,4}}] = {{{1,2}},{{3,4}}} *)
With[{connectedvertices=ConnectedComponents[ECGrav`GraphFromCliques[facelst]]},
Table[Select[facelst,SubsetQ[i,#]&],{i,connectedvertices}]
];

(* Catch-all Pattern *)
ECGrav`ConnectedComplexComponents[args___]:=(Message[ECGrav`ConnectedComplexComponents::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FractionInLargestComponent*)


(* Primary Pattern *)
ECGrav`FractionInLargestComponent[g_Graph]:=
(*
(*Notes: *)
(*Given a graph, it computes the ratio of the number of vertices in the largest 
connected component to th etotal numbe rof vertices in the graph.  *)*)
With[{numV=VertexCount[g],largestComponentLength=Max[Length/@(ConnectedComponents[g])]},
largestComponentLength/numV
];

(* Overload Pattern *)
ECGrav`FractionInLargestComponent[amat_List]:=ECGrav`FractionInLargestComponent[AdjacencyGraph[amat]];

(* Catch-all Pattern *)
ECGrav`FractionInLargestComponent[args___]:=(Message[ECGrav`FractionInLargestComponent::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FractionInLargestKPathComponent*)


(* Primary Pattern *)
ECGrav`FractionInLargestKPathComponent[g_Graph,k_Integer]:=
(*
(*Given a graph, it computes the ratio of the number of vertices in the largest 
connected component to the total number of vertices in the graph.  *)*)
With[{kpcomps=ECGrav`KpathConnectedComponents[g,k]},
Max[Length/@kpcomps]/VertexCount[g]
];

(* Overload Pattern *)
ECGrav`FractionInLargestKPathComponent[g_Graph]:=With[{k=Length@@FindClique[g]},
If[k==1,Return[1/VertexCount[g]]];(*graph is fully isolated*)
ECGrav`FractionInLargestKPathComponent[g,k-1]
];

(* Overload Pattern *)
ECGrav`FractionInLargestKPathComponent[amat_List,k_Integer]:=ECGrav`FractionInLargestKPathComponent[AdjacencyGraph[amat],k];

(* Overload Pattern *)
ECGrav`FractionInLargestKPathComponent[amat_List]:=ECGrav`FractionInLargestKPathComponent[AdjacencyGraph[amat]];

(* Catch-all Pattern *)
ECGrav`FractionInLargestKPathComponent[args___]:=(Message[ECGrav`FractionInLargestKPathComponent::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueOrder*)


(* Primary Pattern *)
ECGrav`CliqueOrder[g_Graph]:=With[{clqs=FindClique[g,\[Infinity],All]},Length[clqs]];

(* Catch-all Pattern *)
ECGrav`CliqueOrder[args___]:=(Message[ECGrav`CliqueOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*FacetOrder*)


(* Primary Pattern *)
ECGrav`FacetOrder[facetsLst_List]:=Length[facetsLst];

(* Catch-all Pattern *)
ECGrav`FacetOrder[args___]:=(Message[ECGrav`FacetOrder::argerr, args];
$Failed);


(* ::Item:: *)
(*EulerChi*)


(* Primary Pattern *)
ECGrav`EulerChi[facetsLst_List/;(Depth[facetsLst]==3&&Max[facetsLst]>1)]:=
(*(*****************************)
(* Last Updated: 04/02/2026  *)
(*****************************)*)

(*Given the list of facets of a simplicial complex (not necessarily 
clique complex), it finds the Euler Characteristic by counting all simplices.*)
Module[{connectedComponents=ECGrav`ConnectedComplexComponents[facetsLst],sortedfacetsLst,maxDim,altSignVector,fVector},

	If[Length[facetsLst]==1,Return[1]];
	If[Length[connectedComponents]==1,
		sortedfacetsLst=Sort[#]&/@facetsLst;
		maxDim=Max[Length/@facetsLst];
		altSignVector=Table[(-1)^(i+1),{i,1,maxDim}];
		fVector=Table[Length[DeleteDuplicates[Join@@(Subsets[#,{i}]&/@sortedfacetsLst)]],{i,1,maxDim}];
		Return[altSignVector . fVector],
		Return[Total[ECGrav`EulerChi[#]&/@connectedComponents]]
	]
];

(* Overload Pattern *)
ECGrav`EulerChi[Amat_List/;(Depth[Amat]==3&&Max[Amat]<=1)]:=
(*Given adjacency matrix of a graph, it computes its Euler Characteristic by counting all simplices.*)
With[{maxCliques=FindClique[AdjacencyGraph[Amat],\[Infinity],All]},
If[Amat=={{1}},Return[1]];
ECGrav`EulerChi[maxCliques]];

(* Overload Pattern *)
ECGrav`EulerChi[g_Graph]:=
(*Given a graph it computes its Euler Characteristic by counting all simplices.*)
With[{maxCliques=FindClique[g,\[Infinity],All]},
ECGrav`EulerChi[maxCliques]
];

(* Catch-all Pattern *)
ECGrav`EulerChi[args___]:=(Message[ECGrav`EulerChi::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*RmsPurity*)


(* Primary Pattern *)
ECGrav`RmsPurity[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&
	Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1})]:=
(*Given a complex, it computes its mean purity = 1/|F| sum_{fi,fj}(dim(f_i)-dim(fj))^2, 
where fi, fj run over all facets. If there's just one facet it 
outputs 0. *)
If[Length[facetsLst]==1,Return[0],StandardDeviation[Length/@facetsLst]]

(* Overload Pattern *)
ECGrav`RmsPurity[g_Graph]:=
(*Given a graph, it computes its mean purity = 1/|F| sum_{fi,fj}(dim(f_i)-dim(fj))^2, 
where fi, fj run over all maximal cliques. If there's just one maximal clique it 
outputs 0. *)

With[{clqsLst=FindClique[g,\[Infinity],All]},
	If[Length[clqsLst]==1,Return[0]];StandardDeviation[Length/@clqsLst]
];

(* Overload Pattern *)
ECGrav`RmsPurity[Amat_List/;(SymmetricMatrixQ[Amat]&&Sort[DeleteDuplicates[Flatten[Amat]]]=={0,1})]:=
ECGrav`RmsPurity[AdjacencyGraph[Amat]];

(* Catch-all Pattern *)
ECGrav`RmsPurity[args___]:=(Message[ECGrav`RmsPurity::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Branchedness*)


(* Primary Pattern *)
ECGrav`Branchedness[g_Graph,n_Integer]:=
(*
(*Given a graph and a dimension n, it computes the rms value of the degree of 
the n facets, minus 3/2. This quantity will be zero for pseudo manifolds 
if n is the purity of the graph minus one since all codimension 1 
faces have degree 1 or 2. ratio of the number of vertices in the largest 
connected component to th etotal numbe rof vertices in the graph.  *)*)
Module[{clqsLst=FindClique[g,\[Infinity],All],maxDim,nfaces},
maxDim=Max[Length/@clqsLst];
nfaces=Union@@(Subsets[#,{n}]&/@clqsLst);
If[maxDim<=1,Return[2]];(*a bunch of isolated vertices have degree 0 each so the output would be 2*)
Mean[((ECGrav`HyperDeg[g,#]&/@nfaces-3/2)^2-1/4)]
];

(* Overload Pattern *)
ECGrav`Branchedness[g_Graph]:=
With[{dMad=Length@@FindClique[g]},
	ECGrav`Branchedness[g,dMad-1]
];

(* Overload Pattern: adjacency-matrix input *)
ECGrav`Branchedness[amat_List/;(SymmetricMatrixQ[amat]&&SubsetQ[{0,1},DeleteDuplicates[Flatten[amat]]]),n_Integer]:=
	ECGrav`Branchedness[AdjacencyGraph[amat],n];

(* Overload Pattern: adjacency-matrix input *)
ECGrav`Branchedness[amat_List/;(SymmetricMatrixQ[amat]&&SubsetQ[{0,1},DeleteDuplicates[Flatten[amat]]])]:=
	ECGrav`Branchedness[AdjacencyGraph[amat]];

(* Overload Pattern: facet-list input, interpreted as a clique complex *)
ECGrav`Branchedness[facetsLst_List,n_Integer]:=ECGrav`Branchedness[ECGrav`GraphFromCliques[facetsLst],n];

(* Overload Pattern: facet-list input, interpreted as a clique complex *)
ECGrav`Branchedness[facetsLst_List]:=ECGrav`Branchedness[ECGrav`GraphFromCliques[facetsLst]];

(* Catch-all Pattern *)
ECGrav`Branchedness[args___]:=(Message[ECGrav`Branchedness::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Graph Dimensions*)


(* ::Item::Closed:: *)
(*AvgKDim*)


(* Primary Pattern *)
ECGrav`AvgKDim[g_Graph]:=
(*Computes the inductive (Knill) dimension of a graph as defined by Oliver Knill.*)
Block[{n,ne,isolatedVertices,ni,components,joinComponents,joinCount,cliques,temp,clqAssoc,numClqSizes,w,gm,degassoc,iassoc,rassoc,resultassoc,ComputeDeg,ComputeI,ComputeR,curClqs,curClqFaces,curI,curR,curval,curdeg,dimg},
n=VertexCount[g];
ne=EdgeCount[g];
If[n==0,Return[-1]];
If[n==1||ne==0,Return[0]];
If[ne==n (n-1)/2,Return[ n-1]];
If[ne==n (n-1)/2-1, Return[ n-2]];

isolatedVertices=Select[VertexList[g],VertexDegree[g,#]==0&];
ni=Length[isolatedVertices];
If[ni>0,
Return[ECGrav`AvgKDim[Subgraph[g,Complement[VertexList[g],isolatedVertices]]]*(n-ni)/n]
];
If[TreeGraphQ[g],Return[1]];

components=ConnectedGraphComponents[g];
If[Length[components]>1,Return[Sum[VertexCount[components[[i]]]*ECGrav`AvgKDim[components[[i]]],{i,Length[components]}]/n]];

(************************************************************
Check if the graph is a join of subgraphs.
*******************************************************)
joinComponents=ConnectedGraphComponents[GraphComplement[g]];
joinCount=Length[joinComponents];

If[joinCount>1,Return[joinCount-1+Sum[ECGrav`AvgKDim[Subgraph[g,VertexList[i]]],{i,joinComponents}]
]
];


cliques=FindClique[g,\[Infinity],All];
temp=Reap[Do[Sow[i,Length[i]],{i,cliques}]][[2]];
clqAssoc=<|Table[Length[i[[1]]]->i,{i,temp}]|>;
numClqSizes=Length[clqAssoc];
w=Max[Keys[clqAssoc]];


If[numClqSizes==1,Return[w-1]];

degassoc=<||>;
iassoc=<||>;
rassoc=<||>;
resultassoc=<||>;

ComputeDeg[clq1_List]:=Block[{spheres,deg},
deg=Lookup[degassoc,Key[clq1],-1];
If[deg>=0,Return[deg]];
spheres=Table[AdjacencyList[g,v],{v,clq1}];
(*Print["spheres",spheres]*);
deg=Length[Intersection@@spheres];
degassoc[clq1]=deg;
deg
];

ComputeI[clq2_List]:=Block[{spheres,sph,isolatedNodes,isolatedNodesCt},
isolatedNodesCt=Lookup[iassoc,Key[clq2],-1];
If[isolatedNodesCt>=0,Return[isolatedNodesCt]];
spheres=Table[AdjacencyList[g,v],{v,clq2}];
sph=Subgraph[g,Intersection@@spheres];
isolatedNodes=Select[VertexList[sph],VertexDegree[sph,#]==0&];
isolatedNodesCt=Length[isolatedNodes];
iassoc[clq2]=isolatedNodesCt;
isolatedNodesCt
];

ComputeR[clq3_List]:=Block[{k,rval,deg,faces,faceRvals},
rval=Lookup[rassoc,Key[clq3],-1];
If[rval>=0,Return[rval]];

k=Length[clq3];
If[k==1,rval=1/VertexDegree[g,clq3[[1]]];rassoc[clq3]=rval;Return[rval]
];
deg=Lookup[degassoc,Key[clq3],ComputeDeg[clq3]];
faces=Subsets[clq3,{k-1}];
faceRvals=Table[Lookup[rassoc,Key[i],ComputeR[i]],{i,faces}];
rval=Total[faceRvals]/deg;
rassoc[clq3]=rval;
rval
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
dimg
];

(* Overload Pattern *)
ECGrav`AvgKDim[Amat_List]:=With[{g=AdjacencyGraph[Amat]},ECGrav`AvgKDim[g]];

(* Catch-all Pattern *)
ECGrav`AvgKDim[args___]:=(Message[ECGrav`AvgKDim::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*SpectralDim*)


(* Primary Pattern *)
ECGrav`SpectralDim[Amat_List/;ArrayDepth[Amat]>1,s_Integer/;s>1]:=
(* Spectral dimension for all nodes at step s *)
Module[{n=Length[Amat],AmatPows= MatrixPower[Amat,s],totTable,Ps,avgP},

totTable=Table[Total[AmatPows[[i]]],{i,n}];

Ps=Table[If[totTable[[i]]==0,0,(AmatPows[[i,i]]/totTable[[i]])],{i,n}];
avgP=Total[Ps]/n;

-2.0*Log[avgP*1.0]/Log[s*1.0]
];

(* Overload Pattern *)
ECGrav`SpectralDim[g_Graph,s_Integer/;s>1]:=
(* Spectral dimension for all nodes at step s *)
With[{Amat=Normal[AdjacencyMatrix[g]]},
	ECGrav`SpectralDim[Amat,s]
];

(* Catch-all Pattern *)
ECGrav`SpectralDim[args___]:=(Message[ECGrav`SpectralDim::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HausdorffDim*)


(* Primary Pattern *)
ECGrav`HausdorffDim[Amat_List/;ArrayDepth[Amat]>1,Dmat_List, s_Integer/;s>1]:=
(* Hausdorff dimension for the first node at step s when the distance between nodes is 
given by the matrix Dmat*)
Module[{n=Length[Amat],vol,result},
vol=Sum[ECGrav`GVolume[Amat,Dmat,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,Dmat,s,1];(*Vertex 1 picked at random*)*)
result=Log[vol*1.0]/Log[s*1.0];
result
];

(* Overload Pattern *)
ECGrav`HausdorffDim[g_Graph,Dmat_List, s_Integer/;s>1]:=
(* Hausdorff dimension for the first node at step s when the distance between nodes is 
given by the matrix Dmat*)
Module[{n=VertexCount[g],vol,result},
vol=Sum[ECGrav`GVolume[g,Dmat,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,Dmat,s,1];(*Vertex 1 picked at random*)*)
result=Log[vol*1.0]/Log[s*1.0];
result
];

(* Overload Pattern *)
ECGrav`HausdorffDim[Amat_List/;ArrayDepth[Amat]>1,s_Integer/;s>1]:=
(* Hausdorff dimension for all nodes at step s *)
Module[{n=Length[Amat],vol},
n=Length[Amat];
vol=Sum[ECGrav`GVolume[Amat,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,s,1]; (*Vertex 1 picked at random*) *)
Log[vol*1.0]/Log[s*1.0]
];

(* Overload Pattern *)
ECGrav`HausdorffDim[g_Graph,s_Integer/;s>1]:=
(* Hausdorff dimension for all nodes at step s *)
Module[{n=VertexCount[g],vol},
vol=Sum[ECGrav`GVolume[g,s,i],{i,n}]/n; (*The average volume over all nodes*)
(*vol=GVolume[Amat,s,1]; (*Vertex 1 picked at random*) *)
Log[vol*1.0]/Log[s*1.0]
];

(* Catch-all Pattern *)
ECGrav`HausdorffDim[args___]:=(Message[ECGrav`HausdorffDim::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Pure and Geometric Complexes*)


(* ::Item::Closed:: *)
(*DSphereQ*)


(* Primary Pattern *)
ECGrav`DSphereQ[g_Graph]:=
(*Outputs true if the graph is a d-sphere*)
Module[{n=VertexCount[g],vlist=VertexList[g],clqs=FindClique[g,\[Infinity],All],numClqSizes,dim,spheres},
If[EdgeCount[g]==0,If[n==2,Return[True],Return[False]]];
numClqSizes=DeleteDuplicates[Length/@clqs];
If[Length[numClqSizes]>1,Return[False]];
dim=numClqSizes[[1]];
If[dim==2,Return[PathGraphQ[g]&&!TreeGraphQ[g]]];
spheres=Table[ECGrav`Sph[g,i],{i,vlist}];
AllTrue[spheres,ECGrav`DSphereQ[#]&]
];

(* Catch-all Pattern *)
ECGrav`DSphereQ[args___]:=(Message[ECGrav`DSphereQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*DGraphQ*)


(* Primary Pattern *)
ECGrav`DGraphQ[g_Graph]:=
(*Returns true of the graph is a geometric graph and false if not.*)
Module[{n=VertexCount[g],vlist=VertexList[g],clqs=FindClique[g,\[Infinity],All],numClqSizes,dim,spheres},
If[EdgeCount[g]==0,If[n<=2,Return[True],Return[False]]];
If[CompleteGraphQ[g],Return[True]];(*Complete graphs are d-graphs*)
numClqSizes=DeleteDuplicates[Length/@clqs];
If[Length[numClqSizes]>1,Return[False]];
dim=numClqSizes[[1]];
If[dim==2,Return[PathGraphQ[g]]];
spheres=Table[ECGrav`Sph[g,i],{i,vlist}];
AllTrue[spheres,ECGrav`DGraphQ[#]&]
];

(* Catch-all Pattern *)
ECGrav`DGraphQ[args___]:=(Message[ECGrav`DGraphQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*DGraphBoundary*)


(* Primary Pattern *)
ECGrav`DGraphBoundary[g_Graph]:=
(*Outputs the induced subgraph which is the boundary of the graph g. The boudnary of a geometric graph is the induced subgraph over boundary nodes, i.e. those whose unit spheres are contractible, or paths. *)
Module[{vlist=VertexList[g],interiorPoints,bdryPoints},
(*If[!DGraphQ[g],Return["Error! The graph is not geometric."]];*)
interiorPoints=Select[vlist,ECGrav`DSphereQ[ECGrav`Sph[g,#]]&];
bdryPoints=Complement[vlist,interiorPoints];
Return[Subgraph[g,bdryPoints]]
];

(* Overload Pattern *)
ECGrav`DGraphBoundary[Amat_List]:=
(*Outputs the induced subgraph which is the boundary of the graph with adjacency matrix 
Amat. The boudnary of a geometric graph is the induced subgraph over boundary nodes, i.e. those whose unit spheres are contractible, or paths. *)
With[{bdryGraph=ECGrav`DGraphBoundary[AdjacencyGraph[Amat]]},
	If[VertexCount[bdryGraph]==0,Return[{}]];
	Normal[AdjacencyMatrix[bdryGraph]]
];

(* Catch-all Pattern *)
ECGrav`DGraphBoundary[args___]:=(Message[ECGrav`DGraphBoundary::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CombinatorialManifoldQ*)


(* Primary Pattern *)
ECGrav`CombinatorialManifoldQ[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it checks whether or not it is a combinatorial manifold. 
Returns true if a combinatorial manifold and false if not.*)
Module[{purity,codim1faces,links},
If[facelst=={},Return[True]];(*The empty graph is the -1 sphere*)
If[Length[facelst]==1,Return[True]];
(*an n-simplex is a combinatorial manifold, an n-ball *)
If[Length[DeleteDuplicates[Length/@facelst]]>1,Return[False]];(*the complex has to be pure*)

purity=Length[facelst[[1]]];
If[purity==1,
If[Length[facelst]<=2,Return[True],Return[False]]
];
If[purity==2,
If[PathGraphQ[ECGrav`GraphFromCliques[facelst]],Return[True],Return[False]]
];
If[Length[ECGrav`ConnectedComplexComponents[facelst]]>1,Return[False]];(*the complex has to be connected*)
codim1faces=Union@@(Subsets[#,{purity-1}]&/@facelst);
If[AnyTrue[codim1faces,ECGrav`HyperDeg[facelst,#]>2&],Return[False]];
links=ECGrav`Lnk[facelst,{#}]&/@DeleteDuplicates[Flatten[facelst]];
AllTrue[links,ECGrav`CombinatorialManifoldQ[#]&]
];

(* Catch-all Pattern *)
ECGrav`CombinatorialManifoldQ[args___]:=(Message[ECGrav`CombinatorialManifoldQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CombinatorialSphereQ*)


(* Primary Pattern *)
ECGrav`CombinatorialSphereQ[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it checks whether or not it 
is a combinatorial sphere. Returns true if a combinatorial manifold and false if not.*)
Module[{purity,codim1faces,links},
If[facelst=={},Return[True]];(*The empty graph is the -1 sphere*)
If[Length[facelst]==1,Return[False]];
(*an n-simplex is an n-ball, not a sphere *)
If[Length[DeleteDuplicates[Length/@facelst]]>1,Return[False]];(*the complex has to be pure*)

purity=Length[facelst[[1]]];

If[purity==1,
	If[Length[facelst]==2,Return[True],Return[False]]
];

If[purity==2,Return[With[{g=ECGrav`GraphFromCliques[facelst]},
	If[PathGraphQ[g]&&!AcyclicGraphQ[g],True,False]]]
];

If[Length[ECGrav`ConnectedComplexComponents[facelst]]>1,Return[False]];(*the complex has to be connected*)

codim1faces=Union@@(Subsets[#,{purity-1}]&/@facelst);
If[AnyTrue[codim1faces,ECGrav`HyperDeg[facelst,#]!=2&],Return[False]];

links=ECGrav`Lnk[facelst,{#}]&/@DeleteDuplicates[Flatten[facelst]];
AllTrue[links,ECGrav`CombinatorialSphereQ[#]&]
];

(* Catch-all Pattern *)
ECGrav`CombinatorialSphereQ[args___]:=(Message[ECGrav`CombinatorialSphereQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*OrientableCombinatorialManifoldQ*)


(* Primary Pattern *)
ECGrav`OrientableCombinatorialManifoldQ[inputfacetsLst_List]:=
(*Checks whether or not a pseudo-combinatorial manifold given as a set of facets is 
orientable. Does not check whether or not the input is a combinatorial pseudomanifold *)
Module[{p=Length[inputfacetsLst[[1]]],q=Length[inputfacetsLst],facetsLst=Sort/@inputfacetsLst,orientations,RelativeOrientation,edgeList,graph,spanningtree,m1,m2,unexploredEdges},
If[q<=2,Return[True]];
If[p<=2,Return[True]];

RelativeOrientation[simplex1_List,simplex2_List]:=
(*gives the relative orientation between neighboring facets*)
Block[{sharedFace=Intersection[simplex1,simplex2],remaining1,pos1,remaining2,pos2},(*Find the vertex not in the shared face for each simplex*)
If[Length[sharedFace]!=Length[simplex1]-1,Return[0]];(*relative orientation of non-adjacent is set to 0*)
remaining1=First@Complement[simplex1,sharedFace];
remaining2=First@Complement[simplex2,sharedFace];
pos1=Position[simplex1,remaining1][[1,1]];
pos2=Position[simplex2,remaining2][[1,1]];

((-1)^pos1)*Signature[Delete[simplex1,pos1]]*(-1)^pos2*Signature[Delete[simplex2,pos2]]
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

AllTrue[unexploredEdges,RelativeOrientation[orientations[[Key[#[[1]]]]],orientations[[Key[#[[2]]]]]]==-1&]
];


(* Overload Pattern *)
ECGrav`OrientableCombinatorialManifoldQ[g_Graph]:=
(*Checks whether or not the clique complex of the graph g is an orientable. 
Does not check whether or not the input is a combinatorial pseudomanifold *)
Module[{facetsLst=Sort/@FindClique[g,\[Infinity],All],p,q,orientations,RelativeOrientation,edgeList,graph,spanningtree,m1,m2,unexploredEdges},
p=Length[facetsLst[[1]]];
q=Length[facetsLst];
If[q<=2,Return[True]];
If[p<=2,Return[True]];

RelativeOrientation[simplex1_List,simplex2_List]:=
(*gives the relative orientation between neighboring facets*)
Block[{sharedFace=Intersection[simplex1,simplex2],remaining1,pos1,remaining2,pos2},(*Find the vertex not in the shared face for each simplex*)
If[Length[sharedFace]!=Length[simplex1]-1,Return[0]];(*relative orientation of non-adjacent is set to 0*)
remaining1=First@Complement[simplex1,sharedFace];
remaining2=First@Complement[simplex2,sharedFace];
pos1=Position[simplex1,remaining1][[1,1]];
pos2=Position[simplex2,remaining2][[1,1]];

((-1)^pos1)*Signature[Delete[simplex1,pos1]]*(-1)^pos2*Signature[Delete[simplex2,pos2]]
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

AllTrue[unexploredEdges,RelativeOrientation[orientations[[Key[#[[1]]]]],orientations[[Key[#[[2]]]]]]==-1&]
];


(* Catch-all Pattern *)
ECGrav`OrientableCombinatorialManifoldQ[args___]:=(Message[ECGrav`OrientableCombinatorialManifoldQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CombinatorialBoundary*)


(* Primary Pattern *)
ECGrav`CombinatorialBoundary[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it returns the boundary. 
It first checks whether or not the manifold is a combinatorial manifold and returns 
$Failed if the complex is not a combinatorial manifold with a warning. *)
Module[{purity=Length[facelst[[1]]],codim1faces},

If[ECGrav`CombinatorialManifoldQ[facelst],
codim1faces=Union@@(Subsets[#,{purity-1}]&/@facelst);
Select[codim1faces,ECGrav`HyperDeg[facelst,#]==1&],
Message[ECGrav`CombinatorialBoundary::notmanifold];$Failed]

];

(* Catch-all Pattern *)
ECGrav`CombinatorialBoundary[args___]:=(Message[ECGrav`CombinatorialBoundary::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CountHoles*)


(* Primary Pattern *)
ECGrav`CountHoles[g_Graph,k_Integer]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a graph. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@FindClique[g,Infinity,All]},
ResourceFunction["BettiNumbers"][simplices][[k]]
];

(* Overload Pattern *)
ECGrav`CountHoles[g_Graph]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a graph. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@FindClique[g,Infinity,All]},
ResourceFunction["BettiNumbers"][simplices]
];

(* Overload Pattern *)
ECGrav`CountHoles[maxClqLst_List,k_Integer]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a simplicial complex. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@maxClqLst},
ResourceFunction["BettiNumbers"][simplices][[k]]
];

(* Overload Pattern *)
ECGrav`CountHoles[maxClqLst_List]:=
(*Gives the number of k-holes, i.e., the kth Betti number of a simplicial complex. 
Uses Mathematica's built in ResourceFunction["BettiNumbers"]*)
With[{simplices=Simplex/@maxClqLst},
ResourceFunction["BettiNumbers"][simplices]
];

(* Catch-all Pattern *)
ECGrav`CountHoles[args___]:=(Message[ECGrav`CountHoles::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Visualize2DComplex*)


(* Primary Pattern *)
ECGrav`Visualize2DComplex[facelst_List]:=
(*Given an input simplicial complex as a list of facets, it returns a 3D embedding. *)
Module[{facets,cmap,coord,scoords},

cmap=FindGraphIsomorphism[ECGrav`GraphFromCliques[facelst],CanonicalGraph[ECGrav`GraphFromCliques[facelst]]];

facets=facelst/.cmap[[1]];

coord=GraphEmbedding[CanonicalGraph[ECGrav`GraphFromCliques[facets]],"SpringEmbedding",3];
(*"SpringElectricalEmbedding","SpringEmbedding","HighDimensionalEmbedding","CircularEmbedding","SpiralEmbedding","RandomEmbedding","RadialEmbedding","SpectralEmbedding","StarEmbedding"}*)
(*Graphics3D[{Opacity[1.7],Yellow,GraphicsComplex[coord,Polygon[facets]]}];*)
scoords=Map[Part[coord,#]&,facets];
Graphics3D[{Opacity[0.9],Table[Simplex[i],{i,scoords}]}]

];

(* Catch-all Pattern *)
ECGrav`Visualize2DComplex[args___]:=(Message[ECGrav`Visualize2DComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Counting Pure Complexes*)


(* ::Item::Closed:: *)
(*RankComb*)


(* Primary Pattern *)

ECGrav`RankComb[set_List,numLabels_Integer]:=
(***************************************
Given a sorted set which is a subset of the set {0,1,,...,numLabels-1}, it assigns a 
unique integer to the set between 0 and numLabels choose length(set).
E.g., Rank[{0,1},3] = 0, and  Rank[{1,2},3]=3.
***************************************)
Module[{setSize=Length[set],result},
If[AnyTrue[set,#>=numLabels&],Return[]];
result=Sum[Sum[Binomial[numLabels-k-j,setSize-k],{j,1,set[[k]]-k+1}],{k,1,setSize}];
result
];

(* Catch-all Pattern *)
ECGrav`RankComb[args___]:=(Message[ECGrav`RankComb::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*UnrankComb*)


(* Primary Pattern *)
ECGrav`UnrankComb[l_Integer,numLabels_Integer,setSize_Integer]:=
(***************************************
Given an integer l between 0 and (numLabels choose setSize)-1, it assigns a unique sorted 
set which is a sorted setSize-subset of {0,1,,...,numLabels-1}.
E.g., UnrankComb[0,3,2] = {0,1}, and  UnrankComb[2,3,2] = {1,2}.
***************************************)
Module[{result={},lCur,iCur,iNext,lowVal,highVal,stopnum},
If[Binomial[numLabels,setSize]<=l,
Return[];
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

result
];

(* Catch-all Pattern *)
ECGrav`UnrankComb[args___]:=(Message[ECGrav`UnrankComb::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*NumPureComplexes*)


(* Primary Pattern *)
ECGrav`NumPureComplexes[pp_Integer,qq_Integer,nn_Integer]:=
(*Gives the number of vertex labeled pure simplicial complexes of purity p, facet order q, 
and number of vertices n. Computes recursively and uses memoization*)
Block[{ECGrav`NumPureComplexes},

ECGrav`NumPureComplexes[p_Integer,q_Integer,n_Integer]:=ECGrav`NumPureComplexes[p,q,n]=
Which[q<0,0,q==1&&n==p,1,q==1&&n!=p,0,True,
(((Binomial[n,p]-(q-1))/q)*ECGrav`NumPureComplexes[p,q-1,n]
	+(Binomial[n,p]/q)*Sum[Binomial[p,k]*ECGrav`NumPureComplexes[p,q-1,n-k],{k,1,p}])];
	
ECGrav`NumPureComplexes[pp,qq,nn]

];

(* Overload Pattern *)
ECGrav`NumPureComplexes[p_Integer,q_Integer]:=
(*The number of pure simplicial complexes of purity p, clique order q, and number of vertices n*)
With[{n0=Catch[Do[If[Binomial[nn,p]>=q,Throw[nn]],{nn,p,p*q}]]},
	Total[Table[ECGrav`NumPureComplexes[p,q,n],{n,n0,p*q}]]
];

(* Catch-all Pattern *)
ECGrav`NumPureComplexes[args___]:=(Message[ECGrav`NumPureComplexes::argerr, args];
$Failed);


(* ::Section::Closed:: *)
(*Choosing Isomorphism Classes*)


(* ::Subsection::Closed:: *)
(*Graph Isomorphism Classes*)


(* ::Item::Closed:: *)
(*GraphIsContained*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`GraphIsContained[glist_List/;GraphQ[glist[[1]]],g_Graph/;GraphQ[g]]:=
(*
(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)
(*Notes: *)
(*Given a list of graphs glist and a new graph g, it returns True if there is a graph 
isomorphic to g in glist, False otherwise. It stops checking at the first occurence 
of isimorphic graph.*)
If[VertexCount[g]==0,Return[True],
Return[AnyTrue[glist,IsomorphicGraphQ[g,#]&]]
];

(* Overload Pattern *)
ECGrav`GraphIsContained[Amats_List/;SquareMatrixQ[Amats[[1]]],Am_List/;SquareMatrixQ[Am]]:=(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given a list of adjacency matrices  Amats and a new adjacency matrix Am, it returns True if there is a graph isomorphic to Am in glist, False otherwise. It stops checking at the first occurence of isimorphic graph.*)If[Length[Am[[1]]]==0,Return[True],Return[AnyTrue[Amats,IsomorphicGraphQ[AdjacencyGraph[Am],AdjacencyGraph[#]]&]]
];

(* Catch-all Pattern *)
ECGrav`GraphIsContained[args___]:=(Message[ECGrav`GraphIsContained::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicGraphs*)


(* Primary Pattern *)
SetAttributes[ECGrav`ChooseNonIsomorphicGraphs,Flat];
ECGrav`ChooseNonIsomorphicGraphs[li__List/;Length[{li}]>1&&GraphQ[{li}[[1,1]]]]:=
(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given an arbitrary number of lists of graphs, it merges them into one list of graphs where no two graphs belonging to different input lists are isomorphic. This method does not check if the graphs WITHIN a given list are non-isomorphic.*)
Block[{result={li}[[1]],length=Length[{li}]},
Do[
result=Join[result,Select[{li}[[n]],!ECGrav`GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicGraphs[li_List/;GraphQ[li[[1]]]]:=
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
ECGrav`ChooseNonIsomorphicGraphs[li_List/;GraphQ[li[[1,1]]]]:=
(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given a list of lists of graphs, it merges them into one list where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)Block[{result=li[[1]],length=Length[li]},
Do[
result=Join[result,Select[li[[n]],!ECGrav`GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicGraphs[li__List/;Length[{li}]>1&&SquareMatrixQ[{li}[[1,1]]]]:=(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given an arbitrary number of lists of adjacency matrices of graphs, it merges them into one list of graphs where no two graphs belonging to different input lists are isomorphic. This method does not check if the graphs WITHIN a given list are non-isomorphic.*)Block[{result={li}[[1]],length=Length[{li}]},
Do[
result=Join[result,Select[{li}[[n]],!ECGrav`GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicGraphs[li_List/;MatrixQ[li[[1]]]]:=

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
ECGrav`ChooseNonIsomorphicGraphs[li_List/;SquareMatrixQ[li[[1,1]]]]:=(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given a list of lists of adjacency matrices, it merges them into one list where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)
Block[{result=li[[1]],length=Length[li]},
Do[
result=Join[result,Select[li[[n]],!ECGrav`GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Catch-all Pattern *)
ECGrav`ChooseNonIsomorphicGraphs[args___]:=(Message[ECGrav`ChooseNonIsomorphicGraphs::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Clique Complex Isomorphism Classes*)


(* ::Item::Closed:: *)
(*IsContainedClqComp*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`IsContainedClqComp[cmpxesLst_List/;Depth[cmpxesLst]==4,cmpx_List/;(Depth[cmpx]==3||cmpx=={})]:=
(****************************)
(*Last updated: 09/13/2023  *)
(****************************)(*Given a list of clique complexes cmpxesLst (each given as a list of maximal cliques)  and a single other clique complex cmpx, it returns True if there is a graph isomorphic to clq in clqsLst, False otherwise. It stops checking at the first occurence of isimorphic graph.*)
If[Length[cmpx]==0||cmpx=={{}},Return[True],AnyTrue[cmpxesLst,IsomorphicGraphQ[ECGrav`GraphFromCliques[cmpx],ECGrav`GraphFromCliques[#]]&]
];

(* Catch-all Pattern *)
ECGrav`IsContainedClqComp[args___]:=(Message[ECGrav`IsContainedClqComp::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicClqComplexes*)


(* Primary Pattern *)
SetAttributes[ECGrav`ChooseNonIsomorphicClqComplexes,Flat];
ECGrav`ChooseNonIsomorphicClqComplexes[li__List/;(Depth[{li}]==5&&Length[{li}]>1)]:=
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
result=Join[result,Select[{li}[[n]],!ECGrav`IsContainedClqComp[result,#]&]];
,{n,2,length}];

result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicClqComplexes[li_List/;(Depth[li]==4&&Length[{li}]==1)]:=
(*(*****************************)
(* Last Updated: 09/13/2023  *)
(*****************************)*)

(*given a single list of clique complexes (each given as a list of maximal cliques),
 it generates a new list composed of non-isomorphic clique complexes.*)
Block[{result={},remaining=li},


result=Reap[
While[Length[remaining]>0,
Sow[First[remaining]];
remaining=Select[remaining,IsomorphicGraphQ[ECGrav`GraphFromCliques[First[remaining]],ECGrav`GraphFromCliques[#]]==False&];
];
][[2,1]];

result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicClqComplexes[li_List/;(Depth[{li}]==6&&Length[li]>1)]:=
(*(*****************************)
(* Last Updated: 09/30/2025  *)
(*****************************)*)

(*Given a list of lists of graph maximal cliques, it merges them into one list 
where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)Block[{result=li[[1]],length=Length[li]},


If[length==1,Return[result]];

Do[
result=Join[result,Select[li[[n]],!ECGrav`IsContainedClqComp[result,#]&]];
,{n,2,length}];

result
];

(* Catch-all Pattern *)
ECGrav`ChooseNonIsomorphicClqComplexes[args___]:=(Message[ECGrav`ChooseNonIsomorphicClqComplexes::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Simplicial Complex Isomorphism Classes*)


(* ::Item::Closed:: *)
(*IsomorphicSimplicialComplexQ*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`IsomorphicSimplicialComplexQ[c1_List,c2_List]:=
(*
(****************************)
(*Last updated: 03/05/2024  *)
(****************************)
(*Checks whether two pure simplicial complexes given as lists of facets are isomorphic or not(i.e. whether or not there is a bijection of the vertex sets that is also a bijection of the facets.)
It does so by brute force enumeration of all isomorphism between their respective underlying graphs and testing if any of the isomorphisms are bijections of the facets*)
*)
With[{isomorphisms=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[c1],ECGrav`GraphFromCliques[c2],All]]},
If[c1=={}&&c2=={},Return[True]];
If[Length[isomorphisms]==0,Return[False]];

If[Length[Select[isomorphisms, Complement[Sort/@(c1/.#),Sort/@c2]=={}&
]]>0,Return[True]];

False
];

(* Catch-all Pattern *)
ECGrav`IsomorphicSimplicialComplexQ[args___]:=(Message[ECGrav`IsomorphicSimplicialComplexQ::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*IsContainedSimplicialComp*)


(* Primary Pattern *)
ECGrav`IsContainedSimplicialComp[cmpxesLst_List/;Depth[cmpxesLst]==4,cmpx_List/;(Depth[cmpx]==3||cmpx=={})]:=
(*
(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)
(*Given a list of pure simplicial complexes cmpxesLst (each given as a list of facets of equal size)  and a single other pure simplicial complex cmpx, it returns True if there is a complex isomorphic to cmpx in cmpxesLst, False otherwise. It stops checking at the first occurence of isomorphic complex.*)
If[Length[cmpx]==0||cmpx=={{}},Return[True],AnyTrue[cmpxesLst,ECGrav`IsomorphicSimplicialComplexQ[cmpx,#]&]
];

(* Catch-all Pattern *)
ECGrav`IsContainedSimplicialComp[args___]:=(Message[ECGrav`IsContainedSimplicialComp::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicSimplicialComplexes*)


(* Primary Pattern *)
SetAttributes[ECGrav`ChooseNonIsomorphicSimplicialComplexes, Flat];
ECGrav`ChooseNonIsomorphicSimplicialComplexes[li__List/;(Depth[{li}]==5&&Length[{li}]>1)]:=
(*
(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)
(*Given an arbitrary number of lists of facets of simplicial complexes, it merges them into one list of graphs where no two simplexes belonging to different input lists are isomorphic. This method does not check if the complexes WITHIN a given list are non-isomorphic.*)Block[{result={li}[[1]],length=Length[{li}]},


Do[
result=Join[result,Select[{li}[[n]],!ECGrav`IsContainedSimplicialComp[result,#]&]];
,{n,2,length}];

result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicSimplicialComplexes[li_List/;(Depth[li]==4&&Length[{li}]==1)]:=(*(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)

(*given a single list of simplicial complex facet lists, it generates a new list composed of non-isomorphic facet lists.*)Block[{result={},remaining=li},

result=Reap[
While[Length[remaining]>0,
Sow[First[remaining]];
remaining=Select[remaining[[2;;-1]],ECGrav`IsomorphicSimplicialComplexQ[First[remaining],#]==False&];
];
][[2,1]];

result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicSimplicialComplexes[li_List/;(Depth[li]==5)]:=
(*(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)

(*Given a list of lists of pure simplicial complexes, it merges them into one list where no two graphs originally in different sublists are isomorphic. This method does not check whether or not graphs WITHIN a sublist are non-isomorphic.*)
Block[{result=li[[1]],length=Length[li]},


If[length==1,Return[result]];

Do[
result=Join[result,Select[li[[n]],!ECGrav`IsContainedSimplicialComp[result,#]&]];
,{n,2,length}];

result
];

(* Catch-all Pattern *)
ECGrav`ChooseNonIsomorphicSimplicialComplexes[args___]:=(Message[ECGrav`ChooseNonIsomorphicSimplicialComplexes::argerr, args];
$Failed);


(* ::Section::Closed:: *)
(*Automorphism Groups and Orders*)


(* ::Subsection::Closed:: *)
(*Simplicial Complex Automorphism and Facet Automorphism*)


(* ::Item::Closed:: *)
(*SimplicialComplexAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`SimplicialComplexAutomorphismGroupOrderConn[facetsLst_List]:=

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

If[Length[facetsLst]==1,Return[Length[facetsLst[[1]]]!]];


g=ECGrav`GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Return[]
];

VertexFacetDegree[v_Integer]:=Length[Select[facetsLst,MemberQ[#,v]&]];


leafsLst=Table[Select[i,VertexFacetDegree[#]==1&],
{i,facetsLst}];
 

reducedFacetsParentsAsn=<|Table[With[{leafs=leafsLst[[i]]},
If[Length[leafs]>=1,Sort[Join[{leafs[[1]]},Complement[facetsLst[[i]],leafs]]],facetsLst[[i]]]->facetsLst[[i]]],{i,1,Length[facetsLst]}]|>;


reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[Keys[reducedFacetsParentsAsn]]];


relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[Keys[reducedFacetsParentsAsn]],reducedGraph][[1]]];


relabeledreducedFacetsParentsAsn=<|Table[Sort[(i/.relabelingRule)]->reducedFacetsParentsAsn[[Key[i]]],{i,Keys[reducedFacetsParentsAsn]}]|>;


canonicalReducedFacets=Keys[relabeledreducedFacetsParentsAsn];


reducedAutG=GraphAutomorphismGroup[reducedGraph];


thereAreMissingCanonicalReducedFacets=Complement[FindClique[reducedGraph,\[Infinity],All],canonicalReducedFacets,SameTest->(Complement[#1,#2]=={}&)];


isPermutationOfFacets[x_]:=AllTrue[canonicalReducedFacets,With[{permval=Sort[PermutationReplace[#,x]]},
((MemberQ[canonicalReducedFacets,permval])&&(Length[Lookup[relabeledreducedFacetsParentsAsn,Key[#],{}]]==Length[Lookup[relabeledreducedFacetsParentsAsn,Key[permval],{}]]))]&];


If[thereAreMissingCanonicalReducedFacets!={},reducedAutG=PermutationGroup[Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]]
];


GroupOrder[reducedAutG]*Product[Length[i]!,{i,leafsLst}]

];

(* Catch-all Pattern *)
ECGrav`SimplicialComplexAutomorphismGroupOrderConn[args___]:=(Message[ECGrav`SimplicialComplexAutomorphismGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*SimplicialComplexAutomorphismGroupOrder*)


(* Primary Pattern *)
ECGrav`SimplicialComplexAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 2/22/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Given a general simplicial complex (not necessarily clique complex and not necessarily pure) that may or may not be connected, it computes the order of the automorphism group. It uses SimplicialComplexAutomorphismGroupOrderConn to compute the automorphism group orders of each component and then combines them. *)
(*Note: It only works by first relabeling the complex into canonically labeled complex*)
*)
With[{components=Tally[ECGrav`ConnectedComplexComponents[facetsLst],ECGrav`IsomorphicSimplicialComplexQ[#1,#2]&]},


Product[(ECGrav`SimplicialComplexAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
ECGrav`SimplicialComplexAutomorphismGroupOrder[args___]:=(Message[ECGrav`SimplicialComplexAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroup*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`PureComplexAutomorphismGroup[facetsLst_List]:=
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

g=CanonicalGraph[ECGrav`GraphFromCliques[facetsLst]];


(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[facetsLst],g][[1]]];


If[Length[facetsLst]==1,Return[{SymmetricGroup[Length[facetsLst[[1]]]],Reverse/@relabelingRule}]];

facets=Sort/@(facetsLst/.relabelingRule);


autGroup=GraphAutomorphismGroup[g];


doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingfacets,MemberQ[facets,Sort[PermutationReplace[#,x]]]&];

missingfacets=Complement[Join@@(Subsets[#,{Length[facets[[1]]]}]&/@FindClique[g,\[Infinity],All]),facets];


If[missingfacets!={},autGroup=PermutationGroup[Select[GroupElements[autGroup],doesntMixupMissingAndNonMissing[#]&]
];
];


{autGroup,(Reverse/@relabelingRule)}

];

(* Catch-all Pattern *)
ECGrav`PureComplexAutomorphismGroup[args___]:=(Message[ECGrav`PureComplexAutomorphismGroup::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`PureComplexAutomorphismGroupOrderConn[facetsLst_List]:=

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


If[Length[facetsLst]==1,Return[Length[facetsLst[[1]]]!]];

If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Return[]];

g=ECGrav`GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Return[]
];

VertexFacetDegree[v_Integer]:=Length[Select[facetsLst,MemberQ[#,v]&]];


leafsLst=Table[Select[i,VertexFacetDegree[#]==1&],
{i,facetsLst}];
 

reducedFacets=Table[With[{leafs=Select[i,VertexFacetDegree[#]==1&]},
If[Length[leafs]>=1,Sort[Join[{leafs[[1]]},Complement[i,leafs]]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];


thereAreMissingCanonicalReducedFacets=Complement[FindClique[reducedGraph,\[Infinity],All],canonicalReducedFacets,SameTest->(Complement[#1,#2]=={}&)];


isPermutationOfFacets[x_]:=AllTrue[canonicalReducedFacets,MemberQ[canonicalReducedFacets,Sort[PermutationReplace[#,x]]]&];


If[thereAreMissingCanonicalReducedFacets!={},reducedAutGorder=Length[Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]],reducedAutGorder=GroupOrder[reducedAutG]
];


reducedAutGorder*Product[Length[i]!,{i,leafsLst}]

];

(* Catch-all Pattern *)
ECGrav`PureComplexAutomorphismGroupOrderConn[args___]:=(Message[ECGrav`PureComplexAutomorphismGroupOrderConn, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroupOrder*)


(* Primary Pattern *)
ECGrav`PureComplexAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 1/19/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Given a pure simplicial complex (not necessarily clique complex) that may or may not be connected, it computes the automorphism group order. It uses PureComplexAutomorphismGroupOrderConn to compute the automorphism group orders of each component and then combines them. *)
(*Note: It only works by first relabeling the complex into canonically labeled complex, so the output group will have different labeleing convention than the input unless the input is canonically labeled*)
*)
With[{components=Tally[ECGrav`ConnectedComplexComponents[facetsLst],ECGrav`IsomorphicSimplicialComplexQ[#1,#2]&]},


Product[(ECGrav`PureComplexAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
ECGrav`PureComplexAutomorphismGroupOrder[args___]:=(Message[ECGrav`PureComplexAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetStabilizerGroupOrder*)


(* Primary Pattern *)
ECGrav`PureComplexFacetStabilizerGroupOrder[facetsLst_List]:=
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
ECGrav`PureComplexFacetStabilizerGroupOrder[args___]:=(Message[ECGrav`PureComplexFacetStabilizerGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`PureComplexFacetAutomorphismGroupOrderConn[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 2.0 *) *) 
(*   (* Note:  *) *)
(****************************************)(*computes the facet-automorphism group order of a simpli clique complex whose maximal cliques are given by facetsLst. It does not check whether the complex is a clique complex or not, it assumes it is. It returns Null if it is not connected.*)
*)
Module[{purity,g,reducedFacets,reducedGraph, reducedAutG,canonicalReducedFacets, missingCanonicalReducedFacets,doesntMixupMissingAndNonMissing,reducedFacetStabilizerGroup},

If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Return[]];

If[Length[facetsLst]==1,Return[1]];

g=ECGrav`GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Return[]
];

purity=Length[facetsLst[[1]]];

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];


missingCanonicalReducedFacets=Complement[Join@@(Subsets[#,{purity}]&/@FindClique[reducedGraph,\[Infinity],All]),canonicalReducedFacets];


doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingCanonicalReducedFacets,MemberQ[canonicalReducedFacets,Sort[PermutationReplace[#,x]]]&];


If[missingCanonicalReducedFacets!={},reducedAutG=PermutationGroup[Select[GroupElements[reducedAutG],doesntMixupMissingAndNonMissing[#]&]];
];


reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];


GroupOrder[reducedAutG]/GroupOrder[reducedFacetStabilizerGroup]

];

(* Catch-all Pattern *)
ECGrav`PureComplexFacetAutomorphismGroupOrderConn[args___]:=(Message[ECGrav`PureComplexFacetAutomorphismGroupOrderConn];$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroupOrder*)


(* Primary Pattern *)
ECGrav`PureComplexFacetAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  *) *)
(****************************************)
(*Given a pure simplicial complex (not necessarily a clique complex) as a list of facets, it computes the order of the facet automorphism group. It does not check whether the complex is pure or not, it assumes it is. *)
*)
With[{components=Tally[ECGrav`ConnectedComplexComponents[facetsLst],ECGrav`IsomorphicSimplicialComplexQ[#1,#2]&]},


Product[(ECGrav`PureComplexFacetAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
ECGrav`PureComplexFacetAutomorphismGroupOrder[args___]:=(Message[ECGrav`PureComplexFacetAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroup*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`PureComplexFacetAutomorphismGroup[facetsLst_List]:=
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

If[Length[facetsLst]==1,Return[SymmetricGroup[1]]];

g=CanonicalGraph[ECGrav`GraphFromCliques[facetsLst]];


(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[facetsLst],g][[1]]];


facets=Sort/@(facetsLst/.relabelingRule);


autGroup=GraphAutomorphismGroup[g];


doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingfacets,MemberQ[facets,Sort[PermutationReplace[#,x]]]&];

missingfacets=Complement[Join@@(Subsets[#,{Length[facets[[1]]]}]&/@FindClique[g,\[Infinity],All]),facets];


If[missingfacets!={},
	autGroup=PermutationGroup[Select[GroupElements[autGroup],doesntMixupMissingAndNonMissing[#]&]];
];


facetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[autGroup,#]]&/@facets)];


(*If[GroupOrder[facetStabilizerGroup]==1,Return[autGroup/.(Reverse/@relabelingRule)]];*)

If[GroupOrder[facetStabilizerGroup]==1,Return[autGroup]];

(PermutationGroup[DeleteDuplicates[RightCosetRepresentative[facetStabilizerGroup,#]&/@GroupElements[autGroup]]])/.(Reverse/@relabelingRule)

];

(* Catch-all Pattern *)
ECGrav`PureComplexFacetAutomorphismGroup[args___]:=(Message[ECGrav`PureComplexFacetAutomorphismGroup::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Clique Complex Automorphism and Facet Automorphism*)


(* ::Item::Closed:: *)
(*CliqueFacetStabilixerGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`CliqueFacetStabilizerGroupOrderConn[facetsLst_List]:=
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


If[Length[facetsLst]==1,Return[Length[facetsLst[[1]]]!]];

g=ECGrav`GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Return[]
];

purity=Length[facetsLst[[1]]];

leafsLst=Table[Select[i,VertexDegree[g,#]==purity-1&],
{i,facetsLst}];
 

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];


GroupOrder[reducedFacetStabilizerGroup]*Product[Length[i]!,{i,leafsLst}]

];

(* Catch-all Pattern *)
ECGrav`CliqueFacetStabilizerGroupOrderConn[args___]:=(Message[ECGrav`CliqueFacetStabilizerGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetStabilixerGroupOrder*)


(* Primary Pattern *)
ECGrav`CliqueFacetStabilizerGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  *) *)
(****************************************)
(*Given a clique complex as a list of facets, it computes the order of the facet stabilizer group. It does not check whether the complex is a clique complex or not, it assumes it is. *)
*)With[{components=Tally[ECGrav`ConnectedComplexComponents[facetsLst],IsomorphicGraphQ[ECGrav`GraphFromCliques[#1],ECGrav`GraphFromCliques[#2]]&]},


Product[(ECGrav`CliqueFacetStabilizerGroupOrderConn[i[[1]]]^(i[[2]])),{i,components}]

];

(* Catch-all Pattern *)
ECGrav`CliqueFacetStabilizerGroupOrder[args___]:=(Message[ECGrav`CliqueFacetStabilizerGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`CliqueFacetAutomorphismGroupOrderConn[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/09/2024. *) *)
(*   (* Version: 2.0 *) *) 
(*   (* Note:  *) *)
(****************************************)(*computes the facet-automorphism group order of a connected clique complex whose maximal cliques are given by facetsLst. It does not check whether the complex is a clique complex or not, it assumes it is. It returns Null if it is not connected.*)
*)
Module[{purity,g,reducedFacets,reducedGraph, reducedAutG,canonicalReducedFacets, reducedFacetStabilizerGroup},


If[Length[facetsLst]==1,Return[1]];

g=ECGrav`GraphFromCliques[facetsLst];


If[!ConnectedGraphQ[g],Return[]
];

purity=Length[facetsLst[[1]]];

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];


reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);


reducedAutG=GraphAutomorphismGroup[reducedGraph];

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];


GroupOrder[reducedAutG]/GroupOrder[reducedFacetStabilizerGroup]

];

(* Catch-all Pattern *)
ECGrav`CliqueFacetAutomorphismGroupOrderConn[args___]:=(Message[ECGrav`CliqueFacetAutomorphismGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroupOrder*)


(* Primary Pattern *)

ECGrav`CliqueFacetAutomorphismGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  *) *)
(****************************************)
(*Given a clique complex as a list of facets, it computes the order of the facet automorphism group. It does not check whether the complex is a clique complex or not, it assumes it is. *)
*)With[{components=Tally[ECGrav`ConnectedComplexComponents[facetsLst],IsomorphicGraphQ[ECGrav`GraphFromCliques[#1],ECGrav`GraphFromCliques[#2]]&]},


Product[(ECGrav`CliqueFacetAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
ECGrav`CliqueFacetAutomorphismGroupOrder[args___]:=(Message[ECGrav`CliqueFacetAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroup*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`CliqueFacetAutomorphismGroup[facetsLst_List]:=(*
(****************************************)
(*   (* Last updated 09/30/2025. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)
(*Given a clique complex as a list of facets, it computes the facet automorphism group. 
It does not check whether the complex is a clique complex or not, it assumes it is. *)
*)
Module[{facets,g,relabelingRule,autGroup,facetStabilizerGroup},

If[Length[facetsLst]==1,Return[SymmetricGroup[1]]];

g=CanonicalGraph[ECGrav`GraphFromCliques[facetsLst]];


(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[facetsLst],g][[1]]];


facets=Sort/@(facetsLst/.relabelingRule);


autGroup=GraphAutomorphismGroup[g];


facetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[autGroup,#]]&/@facets)];


If[GroupOrder[facetStabilizerGroup]==1,Return[autGroup]];


(PermutationGroup[DeleteDuplicates[RightCosetRepresentative[facetStabilizerGroup,#]&/@GroupElements[autGroup]]])/.(Reverse/@relabelingRule)

];

(* Catch-all Pattern *)
ECGrav`CliqueFacetAutomorphismGroup[args___]:=(Message[ECGrav`CliqueFacetAutomorphismGroup::argerr, args];
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
ECGrav`GenerateAllVertexLabeledPureSimplicialComplexes[{purity_Integer,facetOrder_Integer, nV_Integer}]:=
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
ECGrav`GenerateAllVertexLabeledPureSimplicialComplexes[{purity_Integer,facetOrder_Integer}]:=
With[{nmin=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]]},
	Join@@ParallelTable[
		ECGrav`GenerateAllVertexLabeledPureSimplicialComplexes[{purity,facetOrder,n}],
	{n,nmin,purity*facetOrder},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}]
]


(* Catch-all Pattern *)
ECGrav`GenerateAllVertexLabeledPureSimplicialComplexes[args___]:=(Message[ECGrav`GenerateAllVertexLabeledPureSimplicialComplexes::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Generate all facet-labeled pure simplicial complexes*)


(* ::Item::Closed:: *)
(*GenerateAllFacetLabeledPureSimplicialComplexes*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`GenerateAllFacetLabeledPureSimplicialComplexes[{purity_Integer,facetOrder_Integer}]:=
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
			ParallelMap[AddOneVertex[#]&,curSXLst,DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}]);
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
ECGrav`GenerateAllFacetLabeledPureSimplicialComplexes[args___]:=(Message[ECGrav`GenerateAllFacetLabeledPureSimplicialComplexes::argerr, args];
$Failed);


(* ::Section:: *)
(*Generate A Random Pure Simplicial Complex*)


(* ::Subsection::Closed:: *)
(*Random [purity*facetOrder]-labeled Pure Simplicial Complex labeled by *)


(* ::Item::Closed:: *)
(*RandomPQLabeledPureSimplicialComplex*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`RandomPQLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
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

ECGrav`RandomPQLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer},numSamples_Integer]:=
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
	
	,{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}]
	];

result
];

(* Catch-all Pattern *)
ECGrav`RandomPQLabeledPureSimplicialComplex[args___]:=(Message[ECGrav`RandomPQLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random Vertex labeled Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomVertexLabeledPureSimplicialComplex*)


(* Primary Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}, numSamples_Integer]:=
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
Module[{pickRandomComposition,buildOneComplex},

(*A subroutine to construct a composition given the total numbe of fvertices*)
pickRandomComposition[]:=
Module[{curVCount=0,newkTable=Table[0,{facetOrder}],newkWeights},

	Do[
		curVCount=nV-Total[newkTable];
		
		Which[q>2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k]-KroneckerDelta[k,0]*(q-1))
					,{k,0,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]],
			q==2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k])
					,{k,1,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[1,purity]]
			];
	,{q,facetOrder,2,-1}];

	newkTable[[1]]=purity;

	newkTable
	
	];

buildOneComplex[]:=
Module[{facets={},composition,allVertices,coveredVertices={},uncoveredVertices,newVertices,newfacet={}},

	(* Pick a random composition *)
	composition=pickRandomComposition[];

	(* Enforce the condition that the composition be an allowed composition *)
	While[AnyTrue[Range[2,facetOrder],Binomial[Sum[composition[[i]],{i,1,#}],purity]<#&],
		composition=pickRandomComposition[]
	];

	(*Construct facets*)
	allVertices=Range[nV];

	coveredVertices={};
	uncoveredVertices=Complement[allVertices,coveredVertices];

	Do[

		newVertices=RandomSample[uncoveredVertices,composition[[k]]];
		newfacet=Sort[Join[RandomSample[coveredVertices,purity-composition[[k]]],newVertices]];

		If[newVertices=={},

			While[MemberQ[facets,newfacet],

				newfacet=Sort[RandomSample[coveredVertices,purity]];
			];
		];

		uncoveredVertices=Complement[uncoveredVertices,newVertices];
		coveredVertices=Join[coveredVertices,newVertices];

		facets=Join[facets,{newfacet}];

	,{k,1,facetOrder}];


	facets

];


	(*Table[buildOneComplex[],{numSamples}]*)

	If[numSamples>10^4,
		ParallelTable[buildOneComplex[],{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}],
		Table[buildOneComplex[],{numSamples}]]
];


(* Overload Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}, numSamples_Integer]:=
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
Module[{nmin,svTable,pickRandomComposition,buildOneComplex},


(*Set nmin*)
nmin=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]];


(*table of sv(p,q,n)*)
svTable=Table[ECGrav`NumPureComplexes[purity,facetOrder,n]*1.0,{n,nmin,purity*facetOrder}];

(*Normalize svTable by the geometric mean for easier computation*)
svTable=svTable/GeometricMean[svTable];

(*A subroutine to construct a composition given the total numbe of fvertices*)
pickRandomComposition[locTotNumVertices_Integer]:=
Module[{curVCount=0,newkTable=Table[0,{facetOrder}],newkWeights},

	Do[
		curVCount=locTotNumVertices-Total[newkTable];
		
		Which[q>2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k]-KroneckerDelta[k,0]*(q-1))
					,{k,0,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]],
			q==2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k])
					,{k,1,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[1,purity]]
			];
	,{q,facetOrder,2,-1}];

	newkTable[[1]]=purity;

	newkTable
	
	];

buildOneComplex[]:=
Module[{facets={},numTotVertices,composition,allVertices,coveredVertices={},uncoveredVertices,newVertices,newfacet={}},

	(*Pick number of vertices*)
	numTotVertices=RandomChoice[svTable->Range[nmin,purity*facetOrder]];

	(* Pick a random composition *)
	composition=pickRandomComposition[numTotVertices];

	(* Enforce the condition that the composition be an allowed composition *)
	While[AnyTrue[Range[2,facetOrder],Binomial[Sum[composition[[i]],{i,1,#}],purity]<#&],
		composition=pickRandomComposition[numTotVertices]
	];

	(*Construct facets*)
	allVertices=Range[numTotVertices];

	coveredVertices={};
	uncoveredVertices=Complement[allVertices,coveredVertices];

	Do[

		newVertices=RandomSample[uncoveredVertices,composition[[k]]];
		newfacet=Sort[Join[RandomSample[coveredVertices,purity-composition[[k]]],newVertices]];

		If[newVertices=={},

			While[MemberQ[facets,newfacet],

				newfacet=Sort[RandomSample[coveredVertices,purity]];
			];
		];

		uncoveredVertices=Complement[uncoveredVertices,newVertices];
		coveredVertices=Join[coveredVertices,newVertices];

		facets=Join[facets,{newfacet}];

	,{k,1,facetOrder}];


	facets

];

If[numSamples>10^4,ParallelTable[buildOneComplex[],{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}],Table[buildOneComplex[],{numSamples}]]

];


(* Overload Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=First[ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder}, 1]];


(* Overload Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}]:=First[ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder, nV}, 1]];


(* Catch-all Pattern *)
ECGrav`RandomVertexLabeledPureSimplicialComplex[args___]:=(Message[ECGrav`RandomVertexLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random Uniform Unlabeled Pure Simplicial Complex from rejection on vertex labeled*)


(* ::Item::Closed:: *)
(*RandomUniformUnlabeledPureSimplicialComplex*)


(* Primary Pattern *)

ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}, numSamples_Integer]:=
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
RandomUniformUnlabeledPureSimplicialComplex[{2,4,6},10] generates 10, 2-pure random simplicial 
complex with 4 edges and 6 vertices.     *)
*)
*)
Module[{pickRandomComposition,buildOneComplex},


(*A subroutine to construct a composition given the total numbe of fvertices*)
pickRandomComposition[]:=
Module[{curVCount=0,newkTable=Table[0,{facetOrder}],newkWeights},

	Do[
		curVCount=nV-Total[newkTable];
		
		Which[q>2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k]-KroneckerDelta[k,0]*(q-1))
					,{k,0,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]],
			q==2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k])
					,{k,1,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[1,purity]]
			];
	,{q,facetOrder,2,-1}];

	newkTable[[1]]=purity;

	newkTable
	
	];

buildOneComplex[]:=
Module[{facets={},composition,allVertices,coveredVertices={},uncoveredVertices,newVertices,newfacet={},weight},

	allVertices=Range[nV];
	While[facets=={},
		
		(* Pick a random composition *)
		composition=pickRandomComposition[];

		(* Enforce the condition that the composition be an allowed composition *)
		While[AnyTrue[Range[2,facetOrder],Binomial[Sum[composition[[i]],{i,1,#}],purity]<#&],
			composition=pickRandomComposition[]
		];

		(*Construct facets*)

		coveredVertices={};
		uncoveredVertices=Complement[allVertices,coveredVertices];

		Do[

			newVertices=RandomSample[uncoveredVertices,composition[[k]]];
			newfacet=Sort[Join[RandomSample[coveredVertices,purity-composition[[k]]],newVertices]];

			If[newVertices=={},

				While[MemberQ[facets,newfacet],

					newfacet=Sort[RandomSample[coveredVertices,purity]];
				];
			];

			uncoveredVertices=Complement[uncoveredVertices,newVertices];
			coveredVertices=Union[coveredVertices,newVertices];

			facets=Join[facets,{newfacet}];

		,{k,1,facetOrder}];
		
		weight=1.0*ECGrav`PureComplexAutomorphismGroupOrder[facets]/nV!;
		
		If[RandomReal[]>weight,facets={}];

	];
	
	facets
];

If[numSamples>20,
	ParallelTable[buildOneComplex[],{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}],
		Table[buildOneComplex[],{numSamples}]]
];


(* Overload Pattern *)

ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}, numSamples_Integer]:=
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
RandomUniformUnlabeledPureSimplicialComplex[{2,4},10] generates 10, 2-pure random 
unlabeled simplicial complex with 4 edges.     *)
*)
*)
Module[{nmin,svTable,pickRandomComposition,buildOneComplex},


(*Set nmin*)
nmin=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]];


(*table of sv(p,q,n)*)
svTable=Table[ECGrav`NumPureComplexes[purity,facetOrder,n]*1.0,{n,nmin,purity*facetOrder}];

(*Normalize svTable by the geometric mean for easier computation*)
svTable=svTable/GeometricMean[svTable];

(*A subroutine to construct a composition given the total numbe of fvertices*)
pickRandomComposition[locTotNumVertices_Integer]:=
Module[{curVCount=0,newkTable=Table[0,{facetOrder}],newkWeights},

	Do[
		curVCount=locTotNumVertices-Total[newkTable];
		
		Which[q>2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k]-KroneckerDelta[k,0]*(q-1))
					,{k,0,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]],
			q==2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k])
					,{k,1,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[1,purity]]
			];
	,{q,facetOrder,2,-1}];

	newkTable[[1]]=purity;

	newkTable
	
	];

buildOneComplex[]:=
Module[{numTotVertices,facets={},composition,allVertices,coveredVertices={},uncoveredVertices,newVertices,newfacet={},weight},

	While[facets=={},

		(*Pick number of vertices*)
		numTotVertices=RandomChoice[svTable->Range[nmin,purity*facetOrder]];

		allVertices=Range[numTotVertices];
		
		(* Pick a random composition *)
		composition=pickRandomComposition[numTotVertices];

		(* Enforce the condition that the composition be an allowed composition *)
		While[AnyTrue[Range[2,facetOrder],Binomial[Sum[composition[[i]],{i,1,#}],purity]<#&],
			composition=pickRandomComposition[numTotVertices]
		];

		(*Construct facets*)

		coveredVertices={};
		uncoveredVertices=Complement[allVertices,coveredVertices];

		Do[

			newVertices=RandomSample[uncoveredVertices,composition[[k]]];
			newfacet=Sort[Join[RandomSample[coveredVertices,purity-composition[[k]]],newVertices]];

			If[newVertices=={},

				While[MemberQ[facets,newfacet],

					newfacet=Sort[RandomSample[coveredVertices,purity]];
				];
			];

			uncoveredVertices=Complement[uncoveredVertices,newVertices];
			coveredVertices=Union[coveredVertices,newVertices];

			facets=Join[facets,{newfacet}];

		,{k,1,facetOrder}];
		
		weight=1.0*ECGrav`PureComplexAutomorphismGroupOrder[facets]/(numTotVertices!);
		
		If[RandomReal[]>weight,facets={}];

	];
	
	facets
];

If[numSamples>20,ParallelTable[buildOneComplex[],{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}],Table[buildOneComplex[],{numSamples}]]

];


(* Overload Pattern *)

ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer,nV_Integer}]:=First[ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{purity,facetOrder,nV}, 1]];


(* Overload Pattern *)

ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=First[ECGrav`RandomUniformUnlabeledPureSimplicialComplex[{purity,facetOrder}, 1]];


(* Catch-all Pattern *)
ECGrav`RandomUniformUnlabeledPureSimplicialComplex[args___]:=(Message[ECGrav`RandomUniformUnlabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random uniform Facet-labeled Pure Simplicial Complex from rejection on vertex labeled*)


(* ::Item::Closed:: *)
(*RandomUniformFacetLabeledPureSimplicialComplex*)


(* Primary Pattern *)

ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}, numSamples_Integer]:=
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
Module[{pickRandomComposition,buildOneComplex},


(*A subroutine to construct a composition given the total number of fvertices*)
pickRandomComposition[]:=
Module[{curVCount=0,newkTable=Table[0,{facetOrder}],newkWeights},

	Do[
		curVCount=nV-Total[newkTable];
		
		Which[q>2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k]-KroneckerDelta[k,0]*(q-1))
					,{k,0,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]],
			q==2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k])
					,{k,1,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[1,purity]]
			];
	,{q,facetOrder,2,-1}];

	newkTable[[1]]=purity;

	newkTable
	
	];

buildOneComplex[]:=
Module[{facets={},composition,allVertices,coveredVertices={},uncoveredVertices,newVertices,newfacet={},acceptanceProb},

	allVertices=Range[nV];
	While[facets=={},
		
		(* Pick a random composition *)
		composition=pickRandomComposition[];

		(* Enforce the condition that the composition be an allowed composition *)
		While[AnyTrue[Range[2,facetOrder],Binomial[Sum[composition[[i]],{i,1,#}],purity]<#&],
			composition=pickRandomComposition[]
		];

		(*Construct facets*)

		coveredVertices={};
		uncoveredVertices=Complement[allVertices,coveredVertices];

		Do[

			newVertices=RandomSample[uncoveredVertices,composition[[k]]];
			newfacet=Sort[Join[RandomSample[coveredVertices,purity-composition[[k]]],newVertices]];

			If[newVertices=={},

				While[MemberQ[facets,newfacet],

					newfacet=Sort[RandomSample[coveredVertices,purity]];
				];
			];

			uncoveredVertices=Complement[uncoveredVertices,newVertices];
			coveredVertices=Union[coveredVertices,newVertices];

			facets=Join[facets,{newfacet}];

		,{k,1,facetOrder}];
		
		(*weight=1.0*(facetOrder!/nV!)*
				(ECGrav`PureComplexAutomorphismGroupOrder[facets]/
					ECGrav`PureComplexFacetAutomorphismGroupOrder[facets]);*)
		(*weight=1.0*(facetOrder!/nV!)*
				(ECGrav`PureComplexFacetStabilizerGroupOrder[facets]);*)
				
		acceptanceProb=Min[1.0,1.0*((Min[facetOrder,nV]!)/(nV!))*
					(ECGrav`PureComplexFacetStabilizerGroupOrder[facets])];
								
		If[RandomReal[]<=acceptanceProb,Return[facets],facets={}];
		

	];
	
	facets
];

If[numSamples>20,
	ParallelTable[buildOneComplex[],{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}],
		Table[buildOneComplex[],{numSamples}]]
];


(* Overload Pattern *)

ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}, numSamples_Integer]:=
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
Module[{nmin,svTable,pickRandomComposition,buildOneComplex},


(*Set nmin*)
nmin=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]];


(*table of sv(p,q,n)*)
svTable=Table[ECGrav`NumPureComplexes[purity,facetOrder,n]*1.0,{n,nmin,purity*facetOrder}];

(*Normalize svTable by the geometric mean for easier computation*)
svTable=svTable/GeometricMean[svTable];

(*A subroutine to construct a composition given the total numbe of fvertices*)
pickRandomComposition[locTotNumVertices_Integer]:=
Module[{curVCount=0,newkTable=Table[0,{facetOrder}],newkWeights},

	Do[
		curVCount=locTotNumVertices-Total[newkTable];
		
		Which[q>2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k]-KroneckerDelta[k,0]*(q-1))
					,{k,0,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]],
			q==2,
				newkWeights=
					Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]
						*(Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k])
					,{k,1,purity}];
				newkTable[[q]]=RandomChoice[newkWeights->Range[1,purity]]
			];
	,{q,facetOrder,2,-1}];

	newkTable[[1]]=purity;

	newkTable
	
	];

buildOneComplex[]:=
Module[{numTotVertices,facets={},composition,allVertices,coveredVertices={},uncoveredVertices,newVertices,newfacet={},acceptanceProb},

	While[facets=={},

		(*Pick number of vertices*)
		numTotVertices=RandomChoice[svTable->Range[nmin,purity*facetOrder]];

		allVertices=Range[numTotVertices];
		
		(* Pick a random composition *)
		composition=pickRandomComposition[numTotVertices];

		(* Enforce the condition that the composition be an allowed composition *)
		While[AnyTrue[Range[2,facetOrder],Binomial[Sum[composition[[i]],{i,1,#}],purity]<#&],
			composition=pickRandomComposition[numTotVertices]
		];

		(*Construct facets*)

		coveredVertices={};
		uncoveredVertices=Complement[allVertices,coveredVertices];

		Do[

			newVertices=RandomSample[uncoveredVertices,composition[[k]]];
			newfacet=Sort[Join[RandomSample[coveredVertices,purity-composition[[k]]],newVertices]];

			If[newVertices=={},

				While[MemberQ[facets,newfacet],

					newfacet=Sort[RandomSample[coveredVertices,purity]];
				];
			];

			uncoveredVertices=Complement[uncoveredVertices,newVertices];
			coveredVertices=Union[coveredVertices,newVertices];

			facets=Join[facets,{newfacet}];

		,{k,1,facetOrder}];
		
		(*weight=1.0*(facetOrder!/(numTotVertices!))*
				(ECGrav`PureComplexAutomorphismGroupOrder[facets]/
					ECGrav`PureComplexFacetAutomorphismGroupOrder[facets]);*)
					
		(*weight=1.0*(facetOrder!/(numTotVertices!))*
					(ECGrav`PureComplexFacetStabilizerGroupOrder[facets]);*)
		
		(*acceptanceProb=Min[1.0,1.0*(facetOrder/(numTotVertices!))*
					(ECGrav`PureComplexFacetStabilizerGroupOrder[facets])];*)
		
		acceptanceProb=Min[1.0,1.0*((Min[facetOrder,nmin]!)/(numTotVertices!))*
					(ECGrav`PureComplexFacetStabilizerGroupOrder[facets])];
								
		If[RandomReal[]<=acceptanceProb,Return[facets],facets={}];

	];
	
	facets
];

If[numSamples>20,ParallelTable[buildOneComplex[],{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}],Table[buildOneComplex[],{numSamples}]]

];


(* Overload Pattern *)

ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer,nV_Integer}]:=First[ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{purity,facetOrder,nV}, 1]];


(* Overload Pattern *)

ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=First[ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[{purity,facetOrder}, 1]];


(* Catch-all Pattern *)
ECGrav`RandomUniformFacetLabeledPureSimplicialComplex[args___]:=(Message[ECGrav`RandomUniformFacetLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Random facet labeled Pure Simplicial Complex (not uniform)*)


(* ::Item::Closed:: *)
(*RandomFacetLabeledPureSimplicialComplex*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`RandomFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
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
ECGrav`RandomFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer},numSamples_Integer]:=
If[numSamples<10^4,
	Table[ECGrav`RandomFacetLabeledPureSimplicialComplex[{purity,facetOrder}],{numSamples}],
	ParallelTable[ECGrav`RandomFacetLabeledPureSimplicialComplex[{purity,facetOrder}],{numSamples},DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`"}]
];

(* Catch-all Pattern *)
ECGrav`RandomFacetLabeledPureSimplicialComplex[args___]:=(Message[ECGrav`RandomFacetLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection:: *)
(*MCMC Random Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexIterMCMCSweep*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`RandomPureSimplicialComplexMCMCSweep[seedComplex_List,labelingChoise_Integer,NN_Integer]:=
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
				*ECGrav`PureComplexFacetStabilizerGroupOrder[seedComplex]*(facetOrder!/(data[[Key["vertexCount"]]]!)),
		labelingChoise==2,
			(1.0/Binomial[purity*facetOrder,data[[Key["vertexCount"]]]])
				*ECGrav`PureComplexAutomorphismGroupOrder[seedComplex]/((data[[Key["vertexCount"]]])!)];

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
				*ECGrav`PureComplexFacetStabilizerGroupOrder[newComplex]*(facetOrder!/newVCount!),
		labelingChoise==2,
			newComplexWeight=(1.0/Binomial[purity*facetOrder,newVCount])
				*ECGrav`PureComplexAutomorphismGroupOrder[newComplex]/(newVCount!)];

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

Remove[purity,facetOrder,step];

data
]


(* Overload Pattern *)

ECGrav`RandomPureSimplicialComplexMCMCSweep[seedComplex_List,HoldNumberOfVerticesFixed_?BooleanQ,labelingChoise_Integer,NN_Integer]:=
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

If[HoldNumberOfVerticesFixed==False,Return[ECGrav`RandomPureSimplicialComplexMCMCSweep[seedComplex,labelingChoise,NN]]];

minNv=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]];

allVertexLabels=Union@@seedComplex;

data=<|"complex"->Sort[Sort/@seedComplex],"edgeCount"->Length[Union@@(Subsets[#,{2}]&/@(Sort/@seedComplex))],"weight"->0.0,"energy"->0.0|>;

data[[Key["weight"]]]=
	Which[
		labelingChoise==0,
			1.0,
		labelingChoise==1,
			ECGrav`PureComplexFacetStabilizerGroupOrder[seedComplex]*(1.0*(facetOrder!)/(nV!)),
		labelingChoise==2,
			1.0*ECGrav`PureComplexAutomorphismGroupOrder[seedComplex]]/(nV!);

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
			newComplexWeight=ECGrav`PureComplexFacetStabilizerGroupOrder[newComplex]*(1.0*(facetOrder!)/(nV!)),
		labelingChoise==2,
			newComplexWeight=1.0*ECGrav`PureComplexAutomorphismGroupOrder[newComplex]/(nV!)];

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

	Remove[purity,facetOrder,data,step];

	result
];


(* Catch-all Pattern *)
ECGrav`RandomPureSimplicialComplexMCMCSweep[args___]:=(Message[ECGrav`RandomPureSimplicialComplexMCMCSweep::argerr, args];
$Failed);


(* ::Item:: *)
(*RandomPureSimplicialComplexMCMCEquilibriate*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[seedComplex_List, labelingChoise_Integer]:=
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
		"minNumV"-><|"complex"->With[{minN=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]]},ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder,minN}]],
				"vertexCount"->Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]],"weight"->1.0,"energy" ->0.0|>,
		"random"-><|"complex"->ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder}],
				"vertexCount"->0,"weight"->1.0,"energy" ->0.0|>|>;

data[["random","vertexCount"]]=Length[DeleteDuplicates[Flatten[data[["random","complex"]]]]];


Which[
	labelingChoise==1,(*facet-labeled*)
		data[["state","weight"]]=ECGrav`PureComplexFacetStabilizerGroupOrder[seedComplex]*(facetOrder!/(data[["state","vertexCount"]]!));
		data[["maxNumV","weight"]]=((purity!)^facetOrder)*(facetOrder!/(purity*facetOrder)!);
		data[["minNumV","weight"]]=ECGrav`PureComplexFacetStabilizerGroupOrder[data[["minNumV","complex"]]]*(facetOrder!/data[["minNumV","vertexCount"]]!);
		data[["random","weight"]]=ECGrav`PureComplexFacetStabilizerGroupOrder[data[["random","complex"]]]*(facetOrder!/data[["random","vertexCount"]]!),
	labelingChoise==2,(*unlabeled*)	
		data[["state","weight"]]=ECGrav`PureComplexAutomorphismGroupOrder[seedComplex]/((data[["state","vertexCount"]])!);
		data[["maxNumV","weight"]]=(facetOrder!*(purity!)^facetOrder)/((purity*facetOrder)!);
		data[["minNumV","weight"]]=ECGrav`PureComplexAutomorphismGroupOrder[data[["minNumV","complex"]]]/(data[["minNumV","vertexCount"]]!);
		data[["random","weight"]]=ECGrav`PureComplexAutomorphismGroupOrder[data[["random","complex"]]]/(data[["random","vertexCount"]]!);
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

			sweepOutput=ECGrav`RandomPureSimplicialComplexMCMCSweep[data[[Key[replicaName],Key["complex"]]],
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

Remove[purity,facetOrder,data,sweepOutput,Entable,outWinLength,inWinLength,AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime,maxNumSweeps];

result

]


(* Overload Pattern *)

ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[seedComplex_List, HoldNumberOfVerticesFixed_?BooleanQ, labelingChoise_Integer]:=
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

If[HoldNumberOfVerticesFixed==False,
	Return[ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[seedComplex, labelingChoise]]
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
			ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,allowedFacetOrderPair[[1]],nV-purity*allowedFacetOrderPair[[2]]}];
		
		unusedVertices=Complement[Range[nV],Union@@highlyConnectedIsland];
		
		remainingFacets=Partition[unusedVertices,{purity}];
		Join[highlyConnectedIsland,remainingFacets]
	];

data=<|
	"state"-><|"complex"->Sort[Sort/@seedComplex],"edgeCount"->Length[Union@@(Subsets[#,{2}]&/@(Sort/@seedComplex))],"weight"->1.0,"energy" ->0.0|>,
	"leastEvenlyConnected"-><|"complex"->makeLeastEvenlyConnectedComplex[],"edgeCount"->0,"weight"->1.0,"energy" ->0.0|>,
	"random1"-><|"complex"->ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder,nV}],"edgeCount"->0,"weight"->1.0,"energy" ->0.0|>,
	"random2"-><|"complex"->ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder,nV}],"edgeCount"->0,"weight"->1.0,"energy" ->0.0|>|>;


data[["leastEvenlyConnected","edgeCount"]]=Length[Union@@(Subsets[#,{2}]&/@(Sort/@data[["leastEvenlyConnected","complex"]]))];
data[["random1","edgeCount"]]=Length[Union@@(Subsets[#,{2}]&/@(Sort/@data[["random1","complex"]]))];
data[["random2","edgeCount"]]=Length[Union@@(Subsets[#,{2}]&/@(Sort/@data[["random2","complex"]]))];

Which[

	labelingChoise==1,

	data[["state","weight"]]=ECGrav`PureComplexFacetStabilizerGroupOrder[seedComplex]*(1.0*(facetOrder!)/(nV!));
	data[["leastEvenlyConnected","weight"]]=ECGrav`PureComplexFacetStabilizerGroupOrder[data[["leastEvenlyConnected","complex"]]]*(1.0*(facetOrder!)/(nV!));
	data[["random1","weight"]]=ECGrav`PureComplexFacetStabilizerGroupOrder[data[["random1","complex"]]]*(1.0*(facetOrder!)/(nV!));
	data[["random2","weight"]]=ECGrav`PureComplexFacetStabilizerGroupOrder[data[["random2","complex"]]]*(1.0*(facetOrder!)/(nV!)),

	labelingChoise==2,

	data[["state","weight"]]=1.0*ECGrav`PureComplexAutomorphismGroupOrder[seedComplex]/(nV!);
	data[["leastEvenlyConnected","weight"]]=1.0*ECGrav`PureComplexAutomorphismGroupOrder[data[["leastEvenlyConnected","complex"]]]/(nV!);
	data[["random1","weight"]]=1.0*ECGrav`PureComplexAutomorphismGroupOrder[data[["random1","complex"]]]/(nV!);
	data[["random2","weight"]]=1.0*ECGrav`PureComplexAutomorphismGroupOrder[data[["random2","complex"]]]/(nV!)

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
			sweepOutput=ECGrav`RandomPureSimplicialComplexMCMCSweep[data[[Key[replicaName],Key["complex"]]],HoldNumberOfVerticesFixed,
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

Remove[purity,facetOrder,nV,makeLeastEvenlyConnectedComplex,data,sweepOutput,Entable,outWinLength,inWinLength,AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime,maxNumSweeps];

result

];


(* Catch-all Pattern *)
ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[args___]:=(Message[ECGrav`RandomPureSimplicialComplexMCMCEquilibriate::argerr, args];
$Failed);


(* ::Item:: *)
(*RandomPureSimplicialComplexMCMCCorrelationTime*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[seedComplex_List,eqlT_Integer, operators_/;MatchQ[operators,{__Function}],
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
		sweepOutput=ECGrav`RandomPureSimplicialComplexMCMCSweep[data[[Key["state"],Key["complex"]]],labelingChoise,1];
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

				Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the energy", observablesTable[[1,1]]],

			i==2,
				Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the magnetization", observablesTable[[2,1]]],

			i>2,
				Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::stuck, operators[[i-2]], observablesTable[[i-2,1]]]];
	Nothing],
{i,1,Length[observablesTable]}];

If[fluctuatingObservableIndices=={},
		Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::alldefault,
			DeleteDuplicates[#]&/@observablesTable,data[[Key["corrTValues"]]]],
	corrTable=Table[
		norm=ECGrav`CorrelationTime[0,observablesTable[[i]]];
		If[norm==0.0,norm=1.0];
			(*the time 0 correlation can be 0 sometimes*)
		Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,PrintTemporary["      computing corrT at t = ",t]];
			ECGrav`CorrelationTime[t,observablesTable[[i]]],
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

Remove[purity,facetOrder,sweepOutput,observablesTable,
	fluctuatingObservableIndices,norm,corrTable,tmaxVals,corrTValues,numsweeps];

data

];


(* Overload Pattern *)
ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[seedComplex_List,HoldNumberOfVerticesFixed_?BooleanQ,eqlT_Integer, operators_/;MatchQ[operators,{__Function}],
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

If[HoldNumberOfVerticesFixed==False,Return[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[seedComplex,eqlT, operators,labelingChoise]]];

data=<|"eqlT"->eqlT,"corrT"->2,"corrTValues"->Table[2,{Length[operators]+2}],
	"state"-><|"complex"->Sort[Sort/@seedComplex],"edgeCount"->Length[Union@@(Subsets[#,{2}]&/@(Sort/@seedComplex))],"weight"->1.0,"energy" ->0.0|>|>;

numsweeps=10*eqlT;

PrintTemporary["computing correlation time with numsweeps ",numsweeps];


observablesTable=
Transpose[
	Table[
		If[Mod[i,Ceiling[numsweeps/5.0]]==0,PrintTemporary[" sweepno ",i]];
		sweepOutput=ECGrav`RandomPureSimplicialComplexMCMCSweep[data[[Key["state"],Key["complex"]]],HoldNumberOfVerticesFixed,labelingChoise,1];
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
			Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the energy", observablesTable[[1,1]]],
			i==2,
			Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::stuck, "the magnetization", observablesTable[[2,1]]],
			i>2,
			Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::stuck, operators[[i-2]], observablesTable[[i-2,1]]]];
		Nothing],
{i,1,Length[observablesTable]}];

If[fluctuatingObservableIndices=={},
		Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::alldefault,
			DeleteDuplicates[#]&/@observablesTable,data[[Key["corrTValues"]]]],
	corrTable=
		Table[
			norm=ECGrav`CorrelationTime[0,observablesTable[[i]]];If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)
			Table[
				If[Mod[t,Ceiling[numsweeps/5.0]]==0,PrintTemporary["      computing corrT at t = ",t]];
				ECGrav`CorrelationTime[t,observablesTable[[i]]],
			{t,0,numsweeps-10}]/norm
		,{i,fluctuatingObservableIndices}];
	
	tmaxVals=Table[(FirstPosition[cT,_?(#<=0&),numsweeps-10])[[1]],{cT,corrTable}]; 
	(* A place to stop the integration for calculation of correlation time is when the autocorrelation value first becomes negative*)
	
	
	corrTValues=Table[Max[Ceiling[Sum[corrTable[[i,t]],{t,1,tmaxVals[[i]]}]],2],{i,1,Length[corrTable]}];
		
	data[[Key["corrT"]]]=Max[corrTValues];

	Do[data[[Key["corrTValues"],j]]=corrTValues[[First@Flatten[Position[fluctuatingObservableIndices,j]]]],{j,fluctuatingObservableIndices}];
];


Remove[purity,facetOrder,sweepOutput,observablesTable,
	fluctuatingObservableIndices,norm,corrTable,tmaxVals,corrTValues,numsweeps];

data

];


(* Catch-all Pattern *)
ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[args___]:=(Message[ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexMCMC*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`RandomPureSimplicialComplexMCMC[seedComplex_List, operators_/;MatchQ[operators,{__Function}],
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

data=ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[seedComplex,labelingChoise];

(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[data[["state","complex"]],data[["eqlT"]], operators,
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

		data[[Key["state"]]]=ECGrav`RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],labelingChoise,data[["corrT"]]];

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

Remove[Tempoutput,numsweeps,stopnum,printCase];

{measurements,data}

]


(* Overload Pattern *)
ECGrav`RandomPureSimplicialComplexMCMC[seedComplex_List, HoldNumberOfVerticesFixed_?BooleanQ, operators_/;MatchQ[operators,{__Function}],
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

If[HoldNumberOfVerticesFixed==False,Return[ECGrav`RandomPureSimplicialComplexMCMC[seedComplex, operators,labelingChoise, numSamples]]];

PrintTemporary[" Starting RandomPureSimplicialComplexMCMC "];

(*
(********************)
(*   Equilibriate   *)
(********************)
*)


data=ECGrav`RandomPureSimplicialComplexMCMCEquilibriate[seedComplex,HoldNumberOfVerticesFixed,labelingChoise];


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

	
Tempoutput=ECGrav`RandomPureSimplicialComplexMCMCCorrelationTime[data[["state","complex"]],HoldNumberOfVerticesFixed,data[["eqlT"]], operators,
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

			data[[Key["state"]]]=ECGrav`RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],
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

Remove[Tempoutput,numsweeps,stopnum,printCase];

{measurements,data}

];


(* Overload Pattern *)
ECGrav`RandomPureSimplicialComplexMCMC[priorRunInput_Association, operators_/;MatchQ[operators,{__Function}],
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
			ECGrav`RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],
				labelingChoise,data[["corrT"]]];

		If[Mod[numsweeps,1]==0,
			Sow[Join[Flatten[{numsweeps,data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["vertexCount"]]],Through[operators [data[[Key["state"],Key["complex"]]]]]}],{data[[Key["state"],Key["complex"]]]}]];
		];

	numsweeps++;

	]
][[2,1]];


Remove[Tempoutput,numsweeps,stopnum,printCase];

{measurements,data}

];


(* Overload Pattern *)
ECGrav`RandomPureSimplicialComplexMCMC[priorRunInput_Association,HoldNumberOfVerticesFixed_?BooleanQ, operators_/;MatchQ[operators,{__Function}],
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


If[HoldNumberOfVerticesFixed==False,Return[ECGrav`RandomPureSimplicialComplexMCMC[priorRunInput, operators,labelingChoise, numSamples]]];

data=priorRunInput;

numsweeps=1;
printCase=Floor[(numSamples*1.0)/5.0];

measurements=Reap[
	While[numsweeps<=numSamples,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		data[[Key["state"]]]=ECGrav`RandomPureSimplicialComplexMCMCSweep[data[["state","complex"]],
										HoldNumberOfVerticesFixed,labelingChoise,data[["corrT"]]
							];


		If[Mod[numsweeps,1]==0,
			Sow[Join[Flatten[{numsweeps,data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["edgeCount"]]],Through[operators [data[[Key["state"],Key["complex"]]]]]}],{data[[Key["state"],Key["complex"]]]}]];
		];

		numsweeps++;

	]
][[2,1]];

Remove[Tempoutput,numsweeps,stopnum,printCase];

{measurements,data}

];


(* ::Chapter:: *)
(*Generate A Random Pseudo Manifold*)


(* ::Subsection::Closed:: *)
(*Random Pseudo Manifold Through Successive Facet Addition*)


(* ::Item::Closed:: *)
(*RandomPseudoManifold*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`RandomUnlabeledPseudoManifold[{purity_Integer,facetOrder_Integer,connected_Integer:0}]:=
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

If[purity<2,
If[facetOrder==1,Return[facets]];
Return[]
];

If[facetOrder==0,
Return[]];

If[facetOrder==1,Return[facets]
];


Do[

{facets,pstingSites} =ECGrav`AddRandomUnlabeledFacetToPseudoManifold[facets,pstingSites,connected],{n,1,facetOrder-1}

];

{facets,pstingSites}

];

(* Catch-all Pattern *)
ECGrav`RandomUnlabeledPseudoManifold[args___]:=(Message[ECGrav`RandomUnlabeledPseudoManifold::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*Add  a random unlabeled facet to a pseudo-manifold*)


(* ::Item::Closed:: *)
(*AddRandomUnlabeledFacetToPseudoManifold*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`AddRandomUnlabeledFacetToPseudoManifold[facetsLst_List,apastingSites_List,connected_Integer:0]:=
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

If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Return[]];If[DeleteDuplicates[Length/@facetsLst][[1]]!=purity,
Return[]];


glen=Length[DeleteDuplicates[Flatten[facetsLst]]];

(*facetAutGroup=PureComplexFacetAutomorphismGroup[facets];*)

{autGroup,relabelingRule}={#1,Reverse/@#2}&@@ECGrav`PureComplexAutomorphismGroup[facetsLst];


g=CanonicalGraph[ECGrav`GraphFromCliques[facetsLst]];

(*relabelingRule=Normal[FindGraphIsomorphism[GraphFromCliques[facetsLst],g][[1]]];*)


facets=Sort[Sort/@(facetsLst/.relabelingRule)];

newPastingSites=Sort[Sort/@(apastingSites/.relabelingRule)];


newclqVset=Range[glen+1,glen+purity];(*vertex set for the new clique to be added*)


If[apastingSites=={},Return[{Join[facets,{newclqVset}],Subsets[newclqVset,{purity-1}]}]];
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
If[Length[candPsite]<2,Return[True]];
facetComplements=Complement[Subsets[vset,{purity-1}],candPsite];


NoneTrue[Tuples[{facets,facetComplements}],SubsetQ[#[[1]],#[[2]]]&]
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

newfacets=facets;


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

If[CheckMatchingFunctionness[psiteSet,bsiteSet]==False,(*Print["Can not paste. No function from psite to bsite. Returning 0. "];*)Return[0]
];

(*(** Check 2.  For each member of bsiteApexSet, check that its preimage, i.e, all vertices in psite mapped to the same vertex in bsiteApexSet have distance of greater than 2. **)*)

If[CheckVertexDistance[psiteSet,bsiteSet]==False,(*Print["Two or more vertices of separations less than 3 will be identified by this pasting. Returning 0 "];*)
Return[0]
];


newfacets=Sort[Sort/@(newfacets/.rules)];


AppendTo[newfacets,newclqVset];


newPastingSites=Complement[newPastingSites,psite];


newPastingSites=Sort[Sort/@(newPastingSites/.rules)];


newPastingSites=Join[newPastingSites,Complement[Subsets[newclqVset,{purity-1}],bsite]];


newPastingSites=Sort[Sort/@newPastingSites];


(*tgraph=GraphFromCliques[newfacets];*)


1

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

{Sort[Sort/@(newfacets/.relabelingRules)], Sort[Sort/@(newPastingSites/.relabelingRules)]}

];

(* Catch-all Pattern *)
ECGrav`AddRandomUnlabeledFacetToPseudoManifold[args___]:=(Message[ECGrav`AddRandomUnlabeledFacetToPseudoManifold::argerr, args];
$Failed);


(* ::Title:: *)
(*PureComplexes Protect and End*)


(* End private context *)
End[]

(* Protect exported symbols *)

Protect @@ Names["ECGrav`PureComplexes`*"];

EndPackage[]
