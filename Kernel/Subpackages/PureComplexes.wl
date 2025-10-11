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


(* ::Section:: *)
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
(*Print["vlst ",vertexlist," edgelist ",edgelist];*)
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
(*Print[" facetLabels ",facetLabels," cmlxAsn ",cmlxAsn];*)
Sort[
	Sort/@Table[
			With[{k=Select[facetLabeledVertices,MemberQ[#,q]&]},
				(*Print[" cmlxAsn ",cmlxAsn];*)
				(*Print["    q ",q," k ",k, " facet ",Keys[Select[cmlxAsn,MemberQ[k,#]&]]];*)
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


(* ::Subsection:: *)
(*Spheres and Balls, Links and Stars*)


(* ::Item:: *)
(*Sph*)


(* Primary Pattern *)
ECGrav`Sph[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer]:=
With[{size = Length[Amat],rowcolsToKeep=Flatten[Position[Amat[[i]],1]]},
	(*Print["in case Amat. Amat = ",Amat];*)
	If[i>size,Print[" Warning! attempting to find the sphere around a vertex that is not in the graph!"];
		Print["graph ",Graph[AdjacencyGraph[Amat],VertexLabels->"Name"]," node ",i];
	Return[{}]];
	If[rowcolsToKeep=={},Return[{}]];
	If[size==0,Return[{}]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep]
];

(* Overload Pattern *)
ECGrav`Sph[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer,r_Integer]:=
Module[{rowcolsToKeep,dm},
(*	Print["in case Amat. Amat = ",Amat];*)
	If[Length[Amat]==0,Return[{}]];
	dm=GraphDistanceMatrix[AdjacencyGraph[Amat]];
	(*Print[" Amat ",Amat, " dm ",dm];*)
	rowcolsToKeep=Flatten[Position[dm[[i]],r]];
	If[rowcolsToKeep=={},Return[{}]];
	(*Print["rowcolsToKeep ",rowcolsToKeep];*)
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep]
];

(* Overload Pattern *)
ECGrav`Sph[g_Graph,i_Integer]:=
Module[{},
	(*Print["in case graph. graph = ",g];*)
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Print["Warning! attempting to find the sphere around a vertex that is not in the graph!"];{}
	];
	
	Subgraph[g,AdjacencyList[g,i]]
];

(* Overload Pattern *)
ECGrav`Sph[g_Graph,i_Integer,r_Integer]:=
With[{sphVertices=Select[VertexList[g],GraphDistance[g,i,#]==r&]},
	Print["in case graph. graph = ",g];
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Print["Warning! attempting to find the sphere around a vertex that is not in the graph!"];{}
	];

	Subgraph[g,sphVertices]
];

(* Catch-all Pattern *)
ECGrav`Sph[args___]:=(Message[ECGrav`Sph::argerr, args];
$Failed);


(* ::Item:: *)
(*Bll*)


(* Primary Pattern *)

ECGrav`Bll[Amat_List/;(SymmetricMatrixQ[Amat]&&SubsetQ[{0,1},Sort[DeleteDuplicates[Flatten[Amat]]]]),i_Integer]:=
With[{size = Length[Amat],rowcolsToKeep=Flatten[Position[Amat[[i]],1]]},
	(*Print["in case Amat. Amat = ",Amat];*)
	If[i>size,Print[" Warning! attempting to find the sphere around a vertex that is not in the graph!"];
		Print["graph ",Graph[AdjacencyGraph[Amat],VertexLabels->"Name"]," node ",i];
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
	(*Print["in case Amat. Amat = ",Amat];*)
	If[Length[Amat]==0,Return[{}]];
	dm=GraphDistanceMatrix[AdjacencyGraph[Amat]];
	(*Print[" Amat ",Amat, " dm ",dm];*)
	rowcolsToKeep=Flatten[Position[dm[[1]],_?(#<=r&)]];
	(*Print["rowcolsToKeep ",rowcolsToKeep];*)
	If[rowcolsToKeep=={},Return[{}]];
	Part[Transpose[Amat[[rowcolsToKeep]]],rowcolsToKeep]
];

(* Overload Pattern *)
ECGrav`Bll[g_Graph,i_Integer]:=
(*Unit ball in a graph g at vertex i*)
With[{unitBallVertices=Union[{i},AdjacencyList[g,i]]},
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Print["Warning! attempting to find the sphere around a vertex that is not in the graph!"];{}
	];

Subgraph[g,unitBallVertices]

];

(* Overload Pattern *)
ECGrav`Bll[g_Graph,i_Integer,r_Integer]:=
(*Ball of radius r in a graph g at vertex i*)
With[{ballVertices=Select[VertexList[g],GraphDistance[g,i,#]<=r&]},
	If[VertexCount[g]==0,Return[{}]];
	If[MemberQ[VertexList[g],i]==False,Print["Warning! attempting to find the ball around a vertex that is not in the graph!"];{}
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
If[CompleteGraphQ[Subgraph[g,clq]]==False,Print["Error! the input ",clq," is not a clique in the graph"];Return[]];
Length[Intersection@@spheres]
];

(* Overload Pattern *)
ECGrav`HyperDeg[Amat_List/;(SymmetricMatrixQ[Amat]&&Sort[DeleteDuplicates[Flatten[Amat]]]=={0,1}),clq_List]:=
(*Computes the degree of the clique clq given as a list of vertices. 
Checks whether or not the input is a clique in the graph*)
With[{g=AdjacencyGraph[Amat]},
If[CompleteGraphQ[Subgraph[g,clq]]==False,Print["Error! the input ",clq," is not a clique in the graph"];Return[]];
ECGrav`HyperDeg[g,clq]
];

(* Overload Pattern *)
ECGrav`HyperDeg[facetsLst_List/;(Depth[facetsLst]==3&&Length[facetsLst]>=1&&
	Sort[DeleteDuplicates[Flatten[facetsLst]]]!={0,1}),clq_List]:=
(*Computes the hyperdegree of the face clq given as a list of vertices. 
Checks whether or not the input is a face in the graph*)
With[{lnk=ECGrav`Lnk[facetsLst]},
	If[NoneTrue[facetsLst,SubsetQ[#,clq]&],Print["Error! the input ",clq," 
	is not a face of the complex"];Return[]];
	Total[Length/@lnk]
];

(* Catch-all Pattern *)
ECGrav`HyperDeg[args___]:=(Message[ECGrav`HyperDeg::argerr, args];
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
ECGrav`KpathConnectedComponents[facets_List,facedim_Integer]:=
(*Given a list of facets of a complex and 
dimensionality facedim, this method computes the number of components connected by 
(facedim+1)-paths where those clqdim faces that form a component connected by a 
path of (facedim+1)-dimensional faces will form a component.  
It  uses the built in Mathematica KVertexConnectedComponents function on a 
graph constructed such that the vertices are all clqdim-faces and edges between 
them are 1 if two faces are contained in a bigger clique and 0 otherwise. *)

Module[{destinations,clqGraphAmat},

If[facedim==1,Return[ConnectedComponents[ECGrav`GraphFromCliques[facets]]]];

destinations=Union@@(Subsets[#,{facedim}]&/@facets);

clqGraphAmat = 
	Table[
		Table[
			If[Complement[i,j]=={},0,
				If[(Length[Union[i,j]]==facedim+1)&&AnyTrue[facets,SubsetQ[#,Union[i,j]]&],1,0]
			],{i,destinations}
		],
{j,destinations}];

(Union@@Part[destinations,#])&/@KVertexConnectedComponents[AdjacencyGraph[clqGraphAmat],1]

];

(* Overload Pattern *)
ECGrav`KpathConnectedComponents[g_Graph,clqdim_Integer]:=
With[{clqs=FindClique[g,\[Infinity],All]},
	ECGrav`KpathConnectedComponents[clqs,clqdim]
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
(*Print[" g ",g," k ",k," kpcomps ",kpcomps];*)
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


(* ::Item::Closed:: *)
(*EulerChi*)


(* Primary Pattern *)
ECGrav`EulerChi[facetsLst_List/;(Depth[facetsLst]==3&&Max[facetsLst]>1)]:=
(*Given the list of facets of a simplicial complex (not necessarily 
clique complex), it finds the Euler Characteristic by counting all simplices.*)
Module[{sortedfacetsLst=Sort[#]&/@facetsLst,maxDim=Max[Length/@facetsLst],altSignVector,fVector},
If[facetsLst=={{1}},Return[1]];
altSignVector=Table[(-1)^(i+1),{i,1,maxDim}];
fVector=Table[Length[DeleteDuplicates[Join@@(Subsets[#,{i}]&/@sortedfacetsLst)]],{i,1,maxDim}];
altSignVector . fVector
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
ECGrav`RmsPurity[AdjacencyGraph[amat]];

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

(* Overload Pattern *)
ECGrav`Branchedness[amat_List,n_Integer]:=ECGrav`Branchedness[AdjacencyGraph[amat],n];

(* Overload Pattern *)
ECGrav`Branchedness[amat_List]:=ECGrav`Branchedness[AdjacencyGraph[amat]];

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
(*Print["temp ",temp];*)
clqAssoc=<|Table[Length[i[[1]]]->i,{i,temp}]|>;
numClqSizes=Length[clqAssoc];
w=Max[Keys[clqAssoc]];


If[numClqSizes==1,Return[w-1]];

degassoc=<||>;
iassoc=<||>;
rassoc=<||>;
resultassoc=<||>;

ComputeDeg[clq1_List]:=Block[{spheres,deg},
(*If[CompleteGraphQ[Subgraph[g,clq1]]\[Equal]False,Print["Error! the input ",clq1," is not a clique in the graph"];Return[]];*)
deg=Lookup[degassoc,Key[clq1],-1];
If[deg>=0,Return[deg]];
spheres=Table[AdjacencyList[g,v],{v,clq1}];
(*Print["spheres",spheres]*);
deg=Length[Intersection@@spheres];
(*Print["  degree ",deg];*)
degassoc[clq1]=deg;
deg
];

ComputeI[clq2_List]:=Block[{spheres,sph,isolatedNodes,isolatedNodesCt},
(*If[CompleteGraphQ[Subgraph[g,clq2]]\[Equal]False,Print["Error! the input ",clq2," is not a clique in the graph"];Return[]];*)
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
(*Geometric Graphs*)


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
(*Print[" EdgeCount ",EdgeCount[g]," n ",n];*)
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
(*Print["interior points ",interiorPoints, "bdryPoints ",bdryPoints];*)
Return[Subgraph[g,bdryPoints]]
];

(* Catch-all Pattern *)
ECGrav`DGraphBoundary[args___]:=(Message[ECGrav`DGraphBoundary::argerr, args];
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
If[AnyTrue[set,#>=numLabels&],Print[" The set can not include a number >= ",numLabels];Return[]];
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
If[Binomial[numLabels,setSize]<=l,Print[" l = ",l," has to be betweeon 0 and numLabels choose setSize - 1 = ",Binomial[numLabels,setSize]-1," exiting "];
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
(*Print["In case 1"];*)
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
(*Print["Beginning with result ",result," remaining ",remaining];*)
(*Print["In case 2"];*)
Reap[
While[Length[remaining]>0,
(*Print["  result ",result," remaining ",remaining];*)
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
(*Print["In case 3"];*)
Do[
(*Print["  result ",result," n ",n," li[[n]] ",li[[n]]];*)
result=Join[result,Select[li[[n]],!ECGrav`GraphIsContained[result,#]&]];
,{n,2,length}];
result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicGraphs[li__List/;Length[{li}]>1&&SquareMatrixQ[{li}[[1,1]]]]:=(*(*****************************)
(* Last Updated: 07/26/2023  *)
(*****************************)*)

(*Given an arbitrary number of lists of adjacency matrices of graphs, it merges them into one list of graphs where no two graphs belonging to different input lists are isomorphic. This method does not check if the graphs WITHIN a given list are non-isomorphic.*)Block[{result={li}[[1]],length=Length[{li}]},
(*Print["In case 4"];*)
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
(*Print["Beginning with result ",result," remaining ",remaining];*)
(*Print["In case 5"];*)
Reap[
While[Length[remaining]>0,
(*Print["  result ",result," remaining ",remaining];*)
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
(*Print["In case 6"];*)
Do[
(*Print["  result ",result," n ",n," li[[n]] ",li[[n]]];*)
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

(*Print["In case 7"];*)
(*Print[" Depth[{li}] ",Depth[{li}]," Length[{li}] ",Length[{li}]];*)
(*Print[" result ",result, " length ",length];*)

Do[
(*Print[" n ",n, " {li}[[n]] ",{li}[[n]]];*)
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

(*Print["Beginning with result ",result," remaining ",remaining];*)
(*Print["In case 8"];*)
(*Print[" Depth[li]] ",Depth[li]," Length[{li}] ",Length[{li}]];*)

result=Reap[
While[Length[remaining]>0,
(*Print["  result ",result," remaining ",remaining];*)
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

(*Print["In case 9"];*)
(*Print[" Depth[{li}]] ",Depth[{li}]," Length[li] ",Length[li]];*)

If[length==1,Return[result]];

Do[
(*Print["  result ",result," n ",n," li[[n]] ",li[[n]]];*)
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
(*Print[" isomorphisms ",isomorphisms];*)
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

(*Print["In case 7"];*)
(*Print[" Depth[{li}] ",Depth[{li}]," Length[{li}] ",Length[{li}]];*)
(*Print[" result ",result, " length ",length];*)

Do[
(*Print[" n ",n, " {li}[[n]] ",{li}[[n]]];*)
result=Join[result,Select[{li}[[n]],!ECGrav`IsContainedSimplicialComp[result,#]&]];
,{n,2,length}];

result
];

(* Overload Pattern *)
ECGrav`ChooseNonIsomorphicSimplicialComplexes[li_List/;(Depth[li]==4&&Length[{li}]==1)]:=(*(*****************************)
(* Last Updated: 03/05/2024  *)
(*****************************)*)

(*given a single list of simplicial complex facet lists, it generates a new list composed of non-isomorphic facet lists.*)Block[{result={},remaining=li},
(*Print["Beginning with result ",result," remaining ",remaining];*)
(*Print["In case 8"];*)
(*Print[" Depth[li] ",Depth[li]," Length[{li}] ",Length[{li}]];*)

result=Reap[
While[Length[remaining]>0,
(*Print["  result ",result," remaining ",remaining];*)
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

(*Print["In case 9"];*)
(*Print[" Depth[li]] ",Depth[li]," Length[li] ",Length[li]];*)

If[length==1,Return[result]];

Do[
(*Print["  result ",result," n ",n," li[[n]] ",li[[n]]];*)
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

(*If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Print[" PureComplexFacetStabilizerGroupOrderConn requires a pure complex. Exiting"];Return[]];*)

g=ECGrav`GraphFromCliques[facetsLst];

(*Print["input facets ",facetsLst];*)

(*Print["input complex ",GraphPlot[GraphFromCliques[facetsLst],VertexLabels->Automatic,ImageSize->Large]];*)


If[!ConnectedGraphQ[g],Print[" PureComplexFacetStabilizerGroupOrderConn requires a connected complex. Exiting"];Return[]
];

VertexFacetDegree[v_Integer]:=Length[Select[facetsLst,MemberQ[#,v]&]];



leafsLst=Table[Select[i,VertexFacetDegree[#]==1&],
{i,facetsLst}];
 
(*Print["leafsLst ",leafsLst];*)

reducedFacetsParentsAsn=<|Table[With[{leafs=leafsLst[[i]]},
(*Print[leafs ];*)
If[Length[leafs]>=1,Sort[Join[{leafs[[1]]},Complement[facetsLst[[i]],leafs]]],facetsLst[[i]]]->facetsLst[[i]]],{i,1,Length[facetsLst]}]|>;


(*Print["reducedFacetsParentsAsn ",reducedFacetsParentsAsn];*)


reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[Keys[reducedFacetsParentsAsn]]];

(*Print[" reduced graph ",GraphPlot[reducedGraph,VertexLabels->Automatic,ImageSize->Medium]];*)

relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[Keys[reducedFacetsParentsAsn]],reducedGraph][[1]]];

(*Print[" relabelingRule ",relabelingRule];*)


relabeledreducedFacetsParentsAsn=<|Table[Sort[(i/.relabelingRule)]->reducedFacetsParentsAsn[[Key[i]]],{i,Keys[reducedFacetsParentsAsn]}]|>;

(*Print["reducedFacetsParentsAsn ",reducedFacetsParentsAsn];*)

(*Print["relabeledreducedFacetsParentsAsn ",relabeledreducedFacetsParentsAsn];*)

canonicalReducedFacets=Keys[relabeledreducedFacetsParentsAsn];



(*Print["canonicalReducedFacets ",canonicalReducedFacets," graph ",GraphPlot[reducedGraph,VertexLabels->Automatic,ImageSize->Medium]];*)

reducedAutG=GraphAutomorphismGroup[reducedGraph];

(*Print["Before adjusting for missing facets, reducedAutG ",reducedAutG, " order ",GroupOrder[reducedAutG]];*)

thereAreMissingCanonicalReducedFacets=Complement[FindClique[reducedGraph,\[Infinity],All],canonicalReducedFacets,SameTest->(Complement[#1,#2]=={}&)];

(*Print[" thereAreMissingCanonicalReducedFacets ",thereAreMissingCanonicalReducedFacets];*)

isPermutationOfFacets[x_]:=AllTrue[canonicalReducedFacets,With[{permval=Sort[PermutationReplace[#,x]]},
((MemberQ[canonicalReducedFacets,permval])&&(Length[Lookup[relabeledreducedFacetsParentsAsn,Key[#],{}]]==Length[Lookup[relabeledreducedFacetsParentsAsn,Key[permval],{}]]))]&];

(*Print[" thereAreMissingCanonicalReducedFacets ",thereAreMissingCanonicalReducedFacets, " ok group elements ",Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]];*)

If[thereAreMissingCanonicalReducedFacets!={},reducedAutG=PermutationGroup[Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]]
];

(*Print["After checking for facet permutation, reducedAutG ",reducedAutG, " order ",GroupOrder[reducedAutG]];*)


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

(*Print[" components ",ConnectedComplexComponents[facetsLst], " tally ",components," comp aut group order ",PureComplexAutomorphismGroupOrderConn[#[[1]]]&/@components];*)

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

(*Print["g ",GraphPlot[g,VertexLabels->Automatic]];*)




(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[facetsLst],g][[1]]];


Print[" relabelingRule ", relabelingRule, " inverseRelabelingRule ",Reverse/@relabelingRule];

If[Length[facetsLst]==1,Return[{SymmetricGroup[Length[facetsLst[[1]]]],Reverse/@relabelingRule}]];

facets=Sort/@(facetsLst/.relabelingRule);

(*Print["relabeled facets according to canonical labeling ",facets];*)

(*Print[" original graph ",GraphPlot[GraphFromCliques[facetsLst],VertexLabels\[Rule]Automatic]," relabeled graph ",GraphPlot[g,VertexLabels\[Rule]Automatic]];*)


autGroup=GraphAutomorphismGroup[g];

(*Print["g Aut group ",autGroup, " order ",GroupOrder[autGroup]];*)

doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingfacets,MemberQ[facets,Sort[PermutationReplace[#,x]]]&];

missingfacets=Complement[Join@@(Subsets[#,{Length[facets[[1]]]}]&/@FindClique[g,\[Infinity],All]),facets];

(*Print[" missingfacets ",missingfacets];*)

If[missingfacets!={},autGroup=PermutationGroup[Select[GroupElements[autGroup],doesntMixupMissingAndNonMissing[#]&]
];
];

(*Print["g Aut group ",autGroup, " order ",GroupOrder[autGroup]];*)


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

(*Print["In PureComplexAutomorphismGroupOrderConn, starting with complex ",facetsLst];*)

If[Length[facetsLst]==1,Return[Length[facetsLst[[1]]]!]];

If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Print[" PureComplexFacetStabilizerGroupOrderConn requires a pure complex. Exiting"];Return[]];

g=ECGrav`GraphFromCliques[facetsLst];

(*Print["g ",GraphPlot[g,VertexLabels->Automatic,ImageSize->Large]];*)

If[!ConnectedGraphQ[g],Print[" PureComplexFacetStabilizerGroupOrderConn requires a connected complex. Exiting"];Return[]
];

VertexFacetDegree[v_Integer]:=Length[Select[facetsLst,MemberQ[#,v]&]];


(*Print["facets ",facetsLst];*)


leafsLst=Table[Select[i,VertexFacetDegree[#]==1&],
{i,facetsLst}];
 
(*Print["leafsLst ",leafsLst];*)

reducedFacets=Table[With[{leafs=Select[i,VertexFacetDegree[#]==1&]},
(*Print[leafs ];*)
If[Length[leafs]>=1,Sort[Join[{leafs[[1]]},Complement[i,leafs]]],i]],{i,facetsLst}];

(*Print["reducedFacets ",reducedFacets," graph ",GraphPlot[GraphFromCliques[reducedFacets],VertexLabels->Automatic,ImageSize->Medium]];*)

reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);

(*Print["canonicalReducedFacets ",canonicalReducedFacets," graph ",GraphPlot[reducedGraph,VertexLabels->Automatic,ImageSize->Medium]];*)

reducedAutG=GraphAutomorphismGroup[reducedGraph];

(*Print["Before adjusting for missing facets, reducedAutG ",reducedAutG, " order ",GroupOrder[reducedAutG]];*)

thereAreMissingCanonicalReducedFacets=Complement[FindClique[reducedGraph,\[Infinity],All],canonicalReducedFacets,SameTest->(Complement[#1,#2]=={}&)];


isPermutationOfFacets[x_]:=AllTrue[canonicalReducedFacets,MemberQ[canonicalReducedFacets,Sort[PermutationReplace[#,x]]]&];

(*Print[" thereAreMissingCanonicalReducedFacets ",thereAreMissingCanonicalReducedFacets, " ok group elements ",Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]];*)

If[thereAreMissingCanonicalReducedFacets!={},reducedAutGorder=Length[Select[GroupElements[reducedAutG],isPermutationOfFacets[#]&]],reducedAutGorder=GroupOrder[reducedAutG]
];

(*Print["After adjusting for missing facets, reducedAutGorder ",reducedAutGorder];*)


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

(*Print[" components ",ConnectedComplexComponents[facetsLst], " tally ",components," comp aut group order ",PureComplexAutomorphismGroupOrderConn[#[[1]]]&/@components];*)

Product[(ECGrav`PureComplexAutomorphismGroupOrderConn[i[[1]]]^(i[[2]]))*(i[[2]]!),{i,components}]

];

(* Catch-all Pattern *)
ECGrav`PureComplexAutomorphismGroupOrder[args___]:=(Message[ECGrav`PureComplexAutomorphismGroupOrder::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetStabilizerGroupOrderConn*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`PureComplexFacetStabilizerGroupOrderConn[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Note:  *) *)
(****************************************)(*computes the facet stabilizer group order of a connected clique complex whose maximal cliques are given by facetsLst. It does not check whether the complex is a clique complex or not, it assumes it is. It returns Null if it is not connected.*)
*)
Module[{purity,g,leafsLst,reducedFacets,reducedGraph, reducedAutG,canonicalReducedFacets, missingCanonicalReducedFacets,doesntMixupMissingAndNonMissing,reducedFacetStabilizerGroup},

If[Length[facetsLst]==1,Return[Length[facetsLst[[1]]]!]];

If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Print[" PureComplexFacetStabilizerGroupOrderConn requires a pure complex. Exiting"];Return[]];

g=ECGrav`GraphFromCliques[facetsLst];

(*Print["g ",GraphPlot[g,VertexLabels->Automatic,ImageSize->Large]];*)

If[!ConnectedGraphQ[g],Print[" PureComplexFacetStabilizerGroupOrderConn requires a connected complex. Exiting"];Return[]
];

purity=Length[facetsLst[[1]]];

(*Print["facets ",facetsLst];

Print[" all max cliques ",FindClique[g,\[Infinity],All]];

Print[" missingFacets ",Complement[Sort/@(Join@@(Subsets[#,{purity}]&/@FindClique[g,\[Infinity],All])),facetsLst]];*)



leafsLst=Table[Select[i,VertexDegree[g,#]==purity-1&],
{i,facetsLst}];
 
(*Print["leafsLst ",leafsLst];*)

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
(*Print[leafs ];*)
If[Length[leafs]>=1,Sort[Join[{leafs[[1]]},Complement[i,leafs]]],i]],{i,facetsLst}];

(*Print["reducedFacets ",reducedFacets," graph ",GraphPlot[GraphFromCliques[reducedFacets],VertexLabels->Automatic,ImageSize->Medium]];*)

reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);

(*Print["canonicalReducedFacets ",canonicalReducedFacets," graph ",GraphPlot[reducedGraph,VertexLabels->Automatic,ImageSize->Medium]];*)

reducedAutG=GraphAutomorphismGroup[reducedGraph];

(*Print["Before adjusting for missing facets, reducedAutG ",reducedAutG, " order ",GroupOrder[reducedAutG]];*)

missingCanonicalReducedFacets=Complement[Join@@(Subsets[#,{purity}]&/@FindClique[reducedGraph,\[Infinity],All]),canonicalReducedFacets];

(*Print[" missingCanonicalReducedFacets ",missingCanonicalReducedFacets];*)

doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingCanonicalReducedFacets,MemberQ[canonicalReducedFacets,Sort[PermutationReplace[#,x]]]&];


If[missingCanonicalReducedFacets!={},reducedAutG=PermutationGroup[Select[GroupElements[reducedAutG],doesntMixupMissingAndNonMissing[#]&]];
];

(*Print["After adjusting for missing facets, reducedAutG ",reducedAutG, " order ",GroupOrder[reducedAutG]];*)

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];


(*Print["reducedFacetStabilizerGroup ",reducedFacetStabilizerGroup, " order ",GroupOrder[reducedFacetStabilizerGroup]];*)

(*Print[" Product[Length[i]!,{i,leafsLst}] ",Product[Length[i]!,{i,leafsLst}]];*)

GroupOrder[reducedFacetStabilizerGroup]*Product[Length[i]!,{i,leafsLst}]

];

(* Catch-all Pattern *)
ECGrav`PureComplexFacetStabilizerGroupOrderConn[args___]:=(Message[ECGrav`PureComplexFacetStabilizerGroupOrderConn::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*PureComplexFacetStabilizerGroupOrder*)


(* Primary Pattern *)
ECGrav`PureComplexFacetStabilizerGroupOrder[facetsLst_List]:=
(*
(****************************************)
(*   (* Last updated 02/13/2024. *) *)
(*   (* Version: 1.0 *) *) 
(*   (* Notes:  *) *)
(****************************************)
(*Given a clique complex as a list of facets, it computes the order of the facet stabilizer group. It does not check whether the complex is pure or not, it assumes it is. *)
*)With[{components=Tally[ECGrav`ConnectedComplexComponents[facetsLst],ECGrav`IsomorphicSimplicialComplexQ[#1,#2]&]},

(*Print[" components ",ConnectedComplexComponents[facetsLst], " tally ",components];*)

(*Print[" comp aut group ",CliqueFacetAutomorphismGroupOrderConnV2[#[[1]]]&/@components];*)

Product[(ECGrav`PureComplexFacetStabilizerGroupOrderConn[i[[1]]]^(i[[2]])),{i,components}]

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

If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Print[" FacetAutomorphismGroupOrderConnV1 requires a pure complex. Exiting"];Return[]];

If[Length[facetsLst]==1,Return[1]];

g=ECGrav`GraphFromCliques[facetsLst];

(*Print["g ",GraphPlot[g,VertexLabels->Automatic,ImageSize->Large]];*)

If[!ConnectedGraphQ[g],Print[" FacetAutomorphismGroupOrderConnV1 requires a connected complex. Exiting"];Return[]
];

purity=Length[facetsLst[[1]]];

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
(*Print[leafs ];*)
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];

(*Print["reducedFacets ",reducedFacets," graph ",GraphPlot[GraphFromCliques[reducedFacets],VertexLabels->Automatic,ImageSize->Medium]];*)

reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);

(*Print["canonicalReducedFacets ",canonicalReducedFacets," graph ",GraphPlot3D[reducedGraph,VertexLabels->Automatic,ImageSize->Medium]];*)

reducedAutG=GraphAutomorphismGroup[reducedGraph];


missingCanonicalReducedFacets=Complement[Join@@(Subsets[#,{purity}]&/@FindClique[reducedGraph,\[Infinity],All]),canonicalReducedFacets];

(*Print[" missingCanonicalReducedFacets ",missingCanonicalReducedFacets];*)

doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingCanonicalReducedFacets,MemberQ[canonicalReducedFacets,Sort[PermutationReplace[#,x]]]&];


If[missingCanonicalReducedFacets!={},reducedAutG=PermutationGroup[Select[GroupElements[reducedAutG],doesntMixupMissingAndNonMissing[#]&]];
];

(*Print["reducedAutG ",reducedAutG, " order ",GroupOrder[reducedAutG]];*)

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];

(*Print[ " reducedFacetStabilizerGroup ",reducedFacetStabilizerGroup, " order ",GroupOrder[reducedFacetStabilizerGroup]];*)

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

(*Print[" components ",ConnectedComplexComponents[facetsLst], " tally ",components];*)

(*Print[" comp aut group ",PureComplexFacetAutomorphismGroupOrderConn[#[[1]]]&/@components];*)

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

(*Print["g ",GraphPlot[g,VertexLabels->Automatic]];*)




(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[facetsLst],g][[1]]];


Print[" relabelingRule ", relabelingRule, " inverseRelabelingRule ",Reverse/@relabelingRule];

facets=Sort/@(facetsLst/.relabelingRule);

Print["relabeled facets according to canonical labeling ",facets];

Print[" original graph ",GraphPlot[ECGrav`GraphFromCliques[facetsLst],VertexLabels->Automatic]," relabeled graph ",GraphPlot[g,VertexLabels->Automatic]];


autGroup=GraphAutomorphismGroup[g];

(*Print["g Aut group ",autGroup, " order ",GroupOrder[autGroup]];*)

doesntMixupMissingAndNonMissing[x_]:=NoneTrue[missingfacets,MemberQ[facets,Sort[PermutationReplace[#,x]]]&];

missingfacets=Complement[Join@@(Subsets[#,{Length[facets[[1]]]}]&/@FindClique[g,\[Infinity],All]),facets];

Print[" missingfacets ",missingfacets];

If[missingfacets!={},autGroup=PermutationGroup[Select[GroupElements[autGroup],doesntMixupMissingAndNonMissing[#]&]
];
];

(*Print["g Aut group ",autGroup, " order ",GroupOrder[autGroup]];*)



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

(*If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Print[" CliqueAutomorphismGroupOrderConn requires a pure complex. Exiting"];Return[]];*)

If[Length[facetsLst]==1,Return[Length[facetsLst[[1]]]!]];

g=ECGrav`GraphFromCliques[facetsLst];

(*Print["g ",GraphPlot3D[g,VertexLabels->Automatic,ImageSize->Large]];*)

If[!ConnectedGraphQ[g],Print[" CliqueAutomorphismGroupOrderConn requires a connected complex. Exiting"];Return[]
];

purity=Length[facetsLst[[1]]];

leafsLst=Table[Select[i,VertexDegree[g,#]==purity-1&],
{i,facetsLst}];
 
(*Print["leafsLst ",leafsLst];*)

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
(*Print[leafs ];*)
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];

(*Print["reducedFacets ",reducedFacets," graph ",GraphPlot[GraphFromCliques[reducedFacets],VertexLabels->Automatic,ImageSize->Medium]];*)

reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);

(*Print["canonicalReducedFacets ",canonicalReducedFacets," graph ",GraphPlot3D[reducedGraph,VertexLabels->Automatic,ImageSize->Medium]];*)

reducedAutG=GraphAutomorphismGroup[reducedGraph];

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];

(*Print[" reducedAutG ",reducedAutG, " order ",reducedAutG, " reducedFacetStabilizerGroup ",reducedFacetStabilizerGroup, " order ",GroupOrder[reducedFacetStabilizerGroup]];*)

(*Print["Product[Length[i]!,{i,leafsLst}] ",Product[Length[i]!,{i,leafsLst}]];*)

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

(*Print[" components ",ConnectedComplexComponents[facetsLst], " tally ",components];*)

(*Print[" comp aut group ",CliqueFacetStabilizerGroupOrderConn[#[[1]]]&/@components];*)

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

(*If[Length[DeleteDuplicates[Length/@facetsLst]]!=1,Print[" CliqueAutomorphismGroupOrderConn requires a pure complex. Exiting"];Return[]];*)

If[Length[facetsLst]==1,Return[1]];

g=ECGrav`GraphFromCliques[facetsLst];

(*Print["g ",GraphPlot[g,VertexLabels->Automatic,ImageSize->Large]];*)

If[!ConnectedGraphQ[g],Print[" CliqueAutomorphismGroupOrderConn requires a connected complex. Exiting"];Return[]
];

purity=Length[facetsLst[[1]]];

reducedFacets=Table[With[{leafs=Select[i,VertexDegree[g,#]==purity-1&]},
(*Print[leafs ];*)
If[Length[leafs]>=1,Join[{leafs[[1]]},Complement[i,leafs]],i]],{i,facetsLst}];

(*Print["reducedFacets ",reducedFacets," graph ",GraphPlot[GraphFromCliques[reducedFacets],VertexLabels->Automatic,ImageSize->Medium]];*)

reducedGraph=CanonicalGraph[ECGrav`GraphFromCliques[reducedFacets]];

canonicalReducedFacets=Sort/@(reducedFacets/.Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[reducedFacets],reducedGraph][[1]]]);

(*Print["canonicalReducedFacets ",canonicalReducedFacets," graph ",GraphPlot[reducedGraph,VertexLabels->Automatic,ImageSize->Medium]];*)

reducedAutG=GraphAutomorphismGroup[reducedGraph];

reducedFacetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[reducedAutG,#]]&/@canonicalReducedFacets)];

(*Print[" reducedAutG ",reducedAutG, " order ",reducedAutG, " reducedFacetStabilizerGroup ",reducedFacetStabilizerGroup, " order ",GroupOrder[reducedFacetStabilizerGroup]];*)

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

(*Print[" components ",ConnectedComplexComponents[facetsLst], " tally ",components," comp aut group ",CliqueFacetAutomorphismGroupOrderConnV2[#[[1]]]&/@components];*)

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

(*Print["g ",GraphPlot[g,VertexLabels->Automatic]];*)




(*The facets have to be labeled from 1 to glen so that the action of the graph automorphism will be well defined.*)

relabelingRule=Normal[FindGraphIsomorphism[ECGrav`GraphFromCliques[facetsLst],g][[1]]];


(*Print[" relabelingRule ", relabelingRule, " inverseRelabelingRule ",Reverse/@relabelingRule];*)

facets=Sort/@(facetsLst/.relabelingRule);

(*Print["relabeled facets according to canonical labeling ",facets];*)

(*Print[" original graph ",GraphPlot[GraphFromCliques[facetsLst],VertexLabels\[Rule]Automatic]," relabeled graph ",GraphPlot[g,VertexLabels\[Rule]Automatic]];*)


autGroup=GraphAutomorphismGroup[g];

(*Print["g Aut group ",autGroup, " order ",GroupOrder[autGroup]];*)


facetStabilizerGroup=PermutationGroup[Intersection@@(GroupElements[GroupSetwiseStabilizer[autGroup,#]]&/@facets)];


If[GroupOrder[facetStabilizerGroup]==1,Return[autGroup]];


(PermutationGroup[DeleteDuplicates[RightCosetRepresentative[facetStabilizerGroup,#]&/@GroupElements[autGroup]]])/.(Reverse/@relabelingRule)

];

(* Catch-all Pattern *)
ECGrav`CliqueFacetAutomorphismGroup[args___]:=(Message[ECGrav`CliqueFacetAutomorphismGroup::argerr, args];
$Failed);


(* ::Chapter:: *)
(*Generate Pure Simplicial Complexes*)


(* ::Section:: *)
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


(*Print["facets ",facets];*)

zeroSeqLength[ind_Integer,cmpsn_List]:=(*counts the number of continuous sequences of 
zeros following the position ind in the composition of n given as a list cmpsn *)
Module[{count=0},
	Do[
		(*Print["j ",j," lst[[j]] ",lst[[j]]];*)
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

(*Print["partitions ",partitions];*)

allowedPartitionPermutations=Join@@Table[Select[Permutations[i],(#[[1]]==purity&&capacityCheckOk[#,purity,facetOrder])&],{i,partitions}];

(*Print["allowedPartitionPermutations ",allowedPartitionPermutations];*)

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

(*Print[ " curlst ",curlst," ptn ",ptn," facetsLst ",facetsLst, " lastpos ",lastpos];*)

Catch[

	If[curCt==clqCt,Throw[curlst],

		(*Print["curCt ",curCt," clqCt ",clqCt, " curVerts ",curVerts];*)

		cands=Select[facetsLst[[lastpos+1;;-1]],
				Length[Complement[#,curVerts]]==ptn[[curCt+1]]&];
				(*all possible new facets that increase the number of vertices by the 
				value of the partition at the next position to curCt*)

		(*Print["cands ", cands];*)


		If[cands!={},
			addOneClique[Append[curlst,#],ptn,facetsLst]&/@cands(*recursive call*)
		,Nothing]

	]
]

];

result=Table[addOneClique[{},j,facets],{j,allowedPartitionPermutations}];

(*Print["result ",result];*)

result=Apply[Join,result,{0,facetOrder-1}];
(*Print["result ",result];*)

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

	(*Print["          Adding one vertex to ",curLstOfVertices];*)

	vertexDegree=<|Table[i->Length[Select[curLstOfVertices,SubsetQ[#,i]&]],{i,facetSubsets}]|>;

	(*Print["               vertexDegree ",vertexDegree];*)

	unAvailableVertices=Table[If[(Length[i]==1&&vertexDegree[[Key[i]]]==purity)||(Length[i]>1&&vertexDegree[[Key[i]]]==purity-1),i,Nothing],{i,facetSubsets}];

	(*Print["               unAvailableVertices ",unAvailableVertices];*)

	availableVertices=Fold[Function[{x,y},Select[x,!SubsetQ[#,y]&]],facetSubsets,unAvailableVertices];

	(*Print["               availableVertices ",availableVertices];*)

	availableVertices=Select[availableVertices,Order[#,curLstOfVertices[[-1]]]>=0&];

	(*Print["               availableVertices ",availableVertices];*)

	(*Print["               Add One returning ",If[availableVertices=={},Nothing,Join[curLstOfVertices,{#}]&/@availableVertices]];*)

	If[availableVertices=={},Nothing,Join[curLstOfVertices,{#}]&/@availableVertices]


];

(*Print["Adding one ", AddOneVertex[{{2,3},{1,2},{1,2}}]];*)

BuildComplex[]:=
	Module[{finishedSXLst={},curSXLst={}, nextSxLst={#}&/@facetSubsets[[facetOrder;;-1]],iterCount=0},

	(*Print[" In BuildComplex , curSXLst ",curSXLst, " nextSxLst ",nextSxLst];*)

	While[nextSxLst!={}&&iterCount<=purity*facetOrder,
		iterCount++;
		curSXLst=nextSxLst;
		(*Print["     In BuildComplex while loop, iterCount ",iterCount];*)

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

(*Print[" complexLst ", complexLst];*)

GetFacetsFromVList[cmpxVlst_List]:=With[{cmlxAsn=<|Table[i->cmpxVlst[[i]],{i,1,Length[cmpxVlst]}]|>},
	Table[
		With[{k=Select[cmpxVlst,MemberQ[#,q]&]},
			(*Print[" cmlxAsn ",cmlxAsn];*)
			(*Print["    q ",q," k ",k, " facet ",Keys[Select[cmlxAsn,MemberQ[k,#]&]]];*)
			Keys[Select[cmlxAsn,MemberQ[k,#]&]]
		]
	,{q,1,facetOrder}]
];

result=Reverse[GetFacetsFromVList[#]]&/@complexLst;

(*Print[" Result ",result];*)

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

(*Print[" facets ",facets," labelsToChoseFrom ",labelsToChoseFrom];*)

Do[
newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

(*Print["newfacet ",newfacet];
Print["memberQ ",MemberQ[facets,newfacet]];*)

iterNum=0;

While[MemberQ[facets,newfacet]&&iterNum<100,
iterNum++;
newfacet=Sort[RandomSample[labelsToChoseFrom,purity]];

(*Print["i ",i," iterNum ",iterNum," facets[[i]] ",facets[[i]]," newfacet ",newfacet," memberQ ",MemberQ[facets,newfacet]]*)
];

(*Print[" i ",i,"labelsToChoseFrom[[1;;(i+1)*purity]]",labelsToChoseFrom[[1;;(i+1)*purity]], " newfacet ",newfacet];*)

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

(*Print[" facets ",facets," labelsToChoseFrom ",labelsToChoseFrom];*)

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

			(*Print["i ",i," iterNum ",iterNum," facets[[i]] ",facets[[i]]," newfacet ",newfacet," memberQ ",MemberQ[facets,newfacet]]*)
			];

			(*Print[" i ",i,"labelsToChoseFrom[[1;;(i+1)*purity]]",labelsToChoseFrom[[1;;(i+1)*purity]], " newfacet ",newfacet];*)

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

				(*Print["i ",i," iterNum ",iterNum," facets[[i]] ",facets[[i]]," newfacet ",newfacet," memberQ ",MemberQ[facets,newfacet]]*)
			];

			(*Print[" i ",i,"labelsToChoseFrom[[1;;(i+1)*purity]]",labelsToChoseFrom[[1;;(i+1)*purity]], " newfacet ",newfacet];*)

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


(* :Code Section: *)

(* Primary Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 10/03/2025. *) *)
(*   (* Note: A totally new approach from 
			previous versions based onrecursive picking of the number of new vertices 
			in the last facet. Very efficient. Gives the correct distribution for 
			vertex labeled pure complexes. *) 
(****************************************)
(*A code to generate a random sample of (possibly disconnected)vertex labeled pure 
simplicial complexes of a given purity, facet order, and number of vertices 
through iterative addition of facets. Outputs a list numSamples of them. Input is a list of three numbers, p,q, 
numSamples, where p is the purity number, q is the facet order, n is the number of vertices
and numSamples is the number of samples required. E.g. 
RandomVertexLabeledPureSimplicialComplex[{2,4,6},10] generates 10, 2-pure random simplicial 
complex with 4 edges and 6 vertices.     *)
(*Algorithm:,  
1. It simply randomly chooses q p-combinations (the facets) successively one at a time by randomly 
	sampling [n]. If at the kth stage the random facet picked has already been picked previously, it tries again., 
3. After generating q facets, if their union doesn't cover [n], it goes back to step 1 and repeats.*)
*)
*)
Module[{buildOneComplex},


buildOneComplex[]:=
	Module[{facets={},allVertices,newkTable,newkWeights,curVCount=0,coveredVertices={},
		uncoveredVertices,newVertices,newfacet={}},


	(*Construct facets*)
	allVertices=Range[nV];

	coveredVertices={};
	uncoveredVertices=Complement[allVertices,coveredVertices];
	curVCount=nV;

	(*Print[""];
	Print[""];
	Print[""];
	Print[""];
	Print[" In buildOneComplex starting  nV ",nV," curVCount ",curVCount," allVertices, ",allVertices, " coveredVertices ",coveredVertices," uncoveredVertices ",uncoveredVertices];*)

	(*Print[""];
	Print[""];*)


	(*Print["    In while loop, before entering do loop, facets ",facets, " covering ",Sort[DeleteDuplicates[Flatten[facets]]]," availableVertices ",availableVertices];*)

	(*Build a partition of n*)
	newkTable=Table[0,{facetOrder}];

	(*Print[" starting with newkTable ",newkTable];*)

	Do[

		curVCount=nV-Total[newkTable];

		(*Print["    In first do loop, at q ",q, " curVCount ",curVCount];*)

		newkWeights=Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]*Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k],{k,0,purity}];

		(*Print["        newkWeights ",newkWeights];*)

		newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]];

		(*Print["        newkTable after picking newk ",newkTable];*)

	,{q,facetOrder,2,-1}
	];

	newkTable[[1]]=purity;

	(*Print[" final newkTable ",newkTable];*)

	Do[

		(*Print[" k ",k," newkTable[[k]] ",newkTable[[k]], " availableVertices ",availableVertices," usedVertices ",usedVertices];*)

		newVertices=RandomSample[uncoveredVertices,newkTable[[k]]];
		newfacet=Sort[Join[RandomSample[coveredVertices,purity-newkTable[[k]]],newVertices]];

		(*Print[" before while "," newVertices ",newVertices," newkTable[[k]] ",newkTable[[k]], " newfacet ",newfacet];*)

		If[newVertices=={},

			While[MemberQ[facets,newfacet],

			newfacet=Sort[RandomSample[coveredVertices,purity]];
			];

			(*Print["iterNum ",iterNum," newVertices ",newVertices," newk ",newk];*)

		];

		uncoveredVertices=Complement[uncoveredVertices,newVertices];
		coveredVertices=Join[coveredVertices,newVertices];

		(*Print["i ",i," newk ",newk," newVertices ",newVertices," newfacet ",newfacet];*)

		facets=Join[facets,{newfacet}];

		(*Print[" i ",i," facets ",facets," vertices ",vertices," vCount ",vCount, " curMinVcount ",curMinVcount];*)

	,{k,1,facetOrder}];


	(*For diagnostics*)
	(*newkTable*)

	facets



	];

	(*Table[buildOneComplex[],{numSamples}]*)

	If[numSamples>10^4,
		ParallelTable[buildOneComplex[],{numSamples}],
		Table[buildOneComplex[],{numSamples}]]
];

(* Overload Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}, numSamples_Integer]:=
(*
(****************************************)
(*   (* Last updated 10/03/2025. *) *)
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
RandomPureSimplicialComplexSimp[2,4,10] generates 10, 2-pure random simplicial 
complex with 4 edges.     *)
(*Algorithm:, 
1. Pick the number of vertices n between n=nmin and n=p*q weighted by sv(p,q,n)., 2. Once numVertices n is picked, it simply randomly chooses q p-combinations (the facets) successively one at a time by randomly sampling [n]. If at the kth stage the random facet picked has already been picked previously, it tries again., 
3. After generating q facets, if their union doesn't cover [n], it goes back to step 1 and repeats.*)
*)
*)
Module[{nmin,svTable,buildOneComplex},


(*Set nmin*)
nmin=Catch[Do[If[Binomial[n,purity]>=facetOrder,Throw[n]],{n,purity,purity*facetOrder}]];

(*Print["nmin ",nmin];*)


(*table of sv(p,q,n)*)
svTable=Table[ECGrav`NumPureComplexes[purity,facetOrder,n]*1.0,{n,nmin,purity*facetOrder}];

(*Print[" Before normalizing svTable ",svTable];*)
(*Normalize svTable by the geometric mean for easier computation*)
svTable=svTable/GeometricMean[svTable];
(*Print[" After normalizing svTable ",svTable];*)


buildOneComplex[]:=Module[{facets={},numTotVertices,allVertices,newkTable,newkWeights,curVCount=0,coveredVertices={},uncoveredVertices,newVertices,newfacet={}},


(*Pick number of vertices*)
numTotVertices=RandomChoice[svTable->Range[nmin,purity*facetOrder]];

(*For Diagnostics*)
(*numTotVertices= 7;*)
(*numTotVertices= purity*facetOrder;*)


(*Print["numTotVertices picked  ",numTotVertices];*)

(*Construct facets*)
allVertices=Range[numTotVertices];

coveredVertices={};
uncoveredVertices=Complement[allVertices,coveredVertices];
curVCount=numTotVertices;

(*Print[""];
Print[""];
Print[""];
Print[""];
Print[" In buildOneComplex starting  numTotVertices ",numTotVertices," curVCount ",curVCount," allVertices, ",allVertices, " coveredVertices ",coveredVertices," uncoveredVertices ",uncoveredVertices];*)



(*Print[""];
Print[""];*)


(*facets=Join[{},{Sort[RandomSample[availableVertices,purity]]}];*)
(*availableVertices=Select[availableVertices,facetDegree[#,facets]<maxFacetDegree&];*)

(*Print["    In while loop, before entering do loop, facets ",facets, " covering ",Sort[DeleteDuplicates[Flatten[facets]]]," availableVertices ",availableVertices];*)

(*Build a partition of n*)
newkTable=Table[0,{facetOrder}];

(*Print[" starting with newkTable ",newkTable];*)

Do[

curVCount=numTotVertices-Total[newkTable];

(*Print["    In first do loop, at q ",q, " curVCount ",curVCount];*)

newkWeights=Table[ECGrav`NumPureComplexes[purity,q-1,curVCount-k]*Binomial[curVCount,curVCount-k]*Binomial[curVCount-k,purity-k],{k,0,purity}];

(*Print["        newkWeights ",newkWeights];*)

newkTable[[q]]=RandomChoice[newkWeights->Range[0,purity]];

(*Print["        newkTable after picking newk ",newkTable];*)

,{q,facetOrder,2,-1}];

newkTable[[1]]=purity;

(*Print[" final newkTable ",newkTable];*)

Do[

(*Print[" k ",k," newkTable[[k]] ",newkTable[[k]], " availableVertices ",availableVertices," usedVertices ",usedVertices];*)

newVertices=RandomSample[uncoveredVertices,newkTable[[k]]];
newfacet=Sort[Join[RandomSample[coveredVertices,purity-newkTable[[k]]],newVertices]];

(*Print[" before while "," newVertices ",newVertices," newkTable[[k]] ",newkTable[[k]], " newfacet ",newfacet];*)

If[newVertices=={},

While[MemberQ[facets,newfacet],

newfacet=Sort[RandomSample[coveredVertices,purity]];
];

(*Print["iterNum ",iterNum," newVertices ",newVertices," newk ",newk];*)

];

uncoveredVertices=Complement[uncoveredVertices,newVertices];
coveredVertices=Join[coveredVertices,newVertices];

(*Print["i ",i," newk ",newk," newVertices ",newVertices," newfacet ",newfacet];*)

facets=Join[facets,{newfacet}];

(*Print[" i ",i," facets ",facets," vertices ",vertices," vCount ",vCount, " curMinVcount ",curMinVcount];*)

,{k,1,facetOrder}];


(*For diagnostics*)
(*newkTable*)

facets



];
(*Table[buildOneComplex[],{numSamples}]*)

If[numSamples>10^4,ParallelTable[buildOneComplex[],{numSamples}],Table[buildOneComplex[],{numSamples}]]
];

(* Overload Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=First[ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder}, 1]];


(* Overload Pattern *)

ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer, nV_Integer}]:=First[ECGrav`RandomVertexLabeledPureSimplicialComplex[{purity,facetOrder, nV}, 1]];


(* Catch-all Pattern *)
ECGrav`RandomVertexLabeledPureSimplicialComplex[args___]:=(Message[ECGrav`RandomVertexLabeledPureSimplicialComplex::argerr, args];
$Failed);


(* ::Subsection:: *)
(*Random facet labeled Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomFacetLabeledPureSimplicialComplex*)


(* :Code Section *)

(* Primary Pattern *)
ECGrav`RandomFacetLabeledPureSimplicialComplex[{purity_Integer,facetOrder_Integer}]:=
(*
(****************************************)
(*   (* Last updated 08/15/2024. *) *)
(*   (* Version: 3 *)*)
(*   (* Note: faster than version 1. Stores degrees in an association .*) *) 
(****************************************)
(*A code to generate a random pure facet labeled(possibly disconnected) simplicial 
complex of a given purity and clique order through iterative addition of vertices 
labeled by the facets. Input is a list of two numbers, {a, b} where a is the purity 
number and b is the number of cliques. E.g. 
ECGrav`RandomFacetLabeledPureSimplicialComplex[{2,4}] generates a 2-pure facet-labeled 
random simplicial complex with 4 edges. It successively choses a random vertex from a 
decreasing list of available vertex space.  *)

*)
Module[{complex={},result,curAvailableVertices=Subsets[Range[facetOrder]][[2;;-1]],curUnAvailableVertices={},vertexDegree,newVertex={}},

(*randomSubset[lst_List]:=Map[If[RandomReal[]<=1/2,#,##&[]]&,lst];*)

vertexDegree=<|Table[i->0,{i,curAvailableVertices}]|>;

While[curAvailableVertices!={},

(*Print[""];
Print["In while loop"];
Print["Before updating, complex ",complex," vertexDegree, ",vertexDegree, " curAvailableVertices ",curAvailableVertices, " curUnAvailableVertices ",curUnAvailableVertices];*)

newVertex=RandomChoice[curAvailableVertices];

Do[vertexDegree[[Key[i]]]++,{i,Subsets[newVertex][[2;;-1]]}];
complex=Join[complex,{newVertex}];

curUnAvailableVertices=Union[curUnAvailableVertices,Table[If[(Length[i]==1&&vertexDegree[[Key[i]]]==purity)||(Length[i]>1&&vertexDegree[[Key[i]]]==purity-1),i,Nothing],{i,curAvailableVertices}]];

curAvailableVertices=Fold[Function[{x,y},Select[x,!SubsetQ[#,y]&]],curAvailableVertices,curUnAvailableVertices];
(*Print[""];
Print["After updating, newVertex ",newVertex," complex ",complex];
Print[" vertex degrees ",vertexDegree];

Print[" curUnAvailableVertices ",curUnAvailableVertices," curAvailableVertices ",curAvailableVertices];*)

];

(*Print[" After while loop, complex ",complex];*)

complex=Sort[complex];
(*Print[" After while loop, complex ",complex];*)

(*Join@@(Map[LexicographicSort[#]&,
GatherBy[ReverseSortBy[Table[Keys[Select[cmlxAsn,MemberQ[#,v]&]],{v,vlist}],Length],Length],{1}])
*)
result=With[{cmlxAsn=<|Table[i->complex[[i]],{i,1,Length[complex]}]|>},
Table[With[{k=Select[complex,MemberQ[#,q]&]},
(*Print[" cmlxAsn ",cmlxAsn];*)
(*Print["    q ",q," k ",k, " facet ",Keys[Select[cmlxAsn,MemberQ[k,#]&]]];*)
Keys[Select[cmlxAsn,MemberQ[k,#]&]]
],{q,1,facetOrder}]
];

(*Print[" result ",result];*)

(*result=Table[With[{k=Select[complex,MemberQ[#,v]&]},

(*Print["    v ",v," k ",k, " position ",Flatten[SubsetPosition[complex,k]]];*)

Flatten[SubsetPosition[complex,k]]],{v,1,facetOrder}];*)

(*Print[" result ",result];*)

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


(* ::Chapter:: *)
(*Generate A Random Pseudo Manifold*)


(* ::Subsection:: *)
(*Random Pseudo Manifold Through Successive Facet Addition*)


(* ::Subsection:: *)
(*Add  a random unlabeled facet to a pseudo-manifold*)


(* ::Title:: *)
(*PureComplexes Protect and End*)


(* End private context *)
End[]

(* Protect exported symbols *)

Protect @@ Names["ECGrav`PureComplexes`*"];

EndPackage[]
