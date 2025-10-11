(* ::Package:: *)

(* ::Input:: *)
(*(*:Name: MCSims - Monte Carlo Simulations for Emergent Combinatorial Gravity*)*)
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
(*Begin MCSims Package*)


BeginPackage["ECGrav`MCSims`"]


(* ::Title:: *)
(*MCSims Public Functions*)


(* ::Title:: *)
(*MCSims Private*)


Begin["`Private`"] (* Begin private context *)


(* ::Chapter:: *)
(*Helper Functions*)


(* ::Section::Closed:: *)
(*Aggregating Data*)


(* :Code Section: *)


(* ::Item::Closed:: *)
(*CorrelationTime*)


(* Primary Pattern*)
ECGrav`CorrelationTime[t_Integer,tbl_List]:=
(*Given a list of values (e.g. magnetization), it computes the correlation at time t.*)
With[{s1=(1/(Length[tbl]-t))*Sum[tbl[[i]]*tbl[[i+t]],{i,1,Length[tbl]-t}],s2=(1/(Length[tbl]-t)^2)*(Sum[tbl[[i]],{i,1,Length[tbl]-t}])*(Sum[tbl[[i+t]],{i,1,Length[tbl]-t}])},
(*
Print["s1 ",s1," s2 ",s2];
Print["s1-s2 ",s1-s2];*)
s1-s2
];

(* Catch-all Pattern *)
ECGrav`CorrelationTime[args___]:=(Message[ECGrav`CorrelationTime::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*EmpCorrelationTime*)


(* Primary Pattern *)
ECGrav`EmpCorrelationTime[t_Integer,tbl_List]:=
(*Computes the Empirical autocorrelation time based on the lecture by Evertz, 2020*)
Block[{tmax=Length[tbl],xt,yt,num, denom},
xt=(1/(tmax-t))*Sum[tbl[[i]],{i,1,tmax-t}];
yt=(1/(tmax-t))*Sum[tbl[[i]],{i,1+t,tmax}];
num=Sum[(tbl[[i]]-xt)(tbl[[i+t]]-yt),{i,1,tmax-t}];
denom=Sqrt[(Sum[(tbl[[i]]-xt)^2,{i,1,tmax-t}])*(Sum[(tbl[[i+t]]-yt)^2,{i,1,tmax-t}])];
(*Print["tmax ",tmax];
Print["s1 ",s1," s2 ",s2];
Print["s1-s2 ",s1-s2];*)
num/denom
];

(* Catch-all Pattern *)
ECGrav`EmpCorrelationTime[args___]:=(Message[ECGrav`EmpCorrelationTime::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ErrorBootstrap*)


(* Primary Pattern *)
ECGrav`ErrorBootstrap[formula_,data_List]:=
(*Computes the uncertainty in the value of formula[data] computed using the 
bootstrap or resampling method. Section 3.4.3 of Newman & Barkema*)
With[{n=Length[data]},
StandardDeviation[ParallelTable[formula[RandomChoice[data,n]],{200},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]]
];

(* Catch-all Pattern *)
ECGrav`ErrorBootstrap[args___]:=(Message[ECGrav`ErrorBootstrap::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*SpecificHeat*)


(* Primary Pattern *)
ECGrav`SpecificHeat[etable_List,NN_Integer,beta_Real]:=
(*specific heat per site given number of sites NN, inverse temperature beta, and a 
table of energies*)
(beta^2/NN)*Variance[etable];

(* Catch-all Pattern *)
ECGrav`SpecificHeat[args___]:=(Message[ECGrav`SpecificHeat::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*DSpecificHeat*)


(* Primary Pattern *)
ECGrav`DSpecificHeat[etable_List,NN_Integer,beta_Real]:=
(*Gives the derivative of specific heat per site wrt temperature given number of sites NN, 
inverse temperature beta, and a table of energies*)
With[{meanE=Mean[etable],meanEsq = Mean[etable^2],meanEcube=Mean[etable^3]},
(-beta^3/NN)*(2.0(meanEsq-meanE^2)(1+beta*meanE)-beta*(meanEcube-meanE*meanEsq ))
];

(* Catch-all Pattern *)
ECGrav`DSpecificHeat[args___]:=(Message[ECGrav`DSpecificHeat::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Susceptibility*)


(* Primary Pattern *)
ECGrav`Susceptibility[obsTable_List,NN_Integer,beta_Real]:=
(*susceptibility per site given number of sites NN, inverse temperature beta, 
and a table of magnetizations*)
(beta*NN)*Variance[obsTable];

(* Catch-all Pattern *)
ECGrav`Susceptibility[args___]:=(Message[ECGrav`Susceptibility::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*LogSumExp*)


(* Primary Pattern *)
(*Computes the log of the sum of exponentials of the list*)
ECGrav`LogSumExp[lst_List]:=With[{mn=Mean[lst]},
mn+Log[Total[Exp[lst-mn]]]
];

(* Catch-all Pattern *)
ECGrav`LogSumExp[args___]:=(Message[ECGrav`LogSumExp::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ComputeMinusBetaTimesFreeEnergy*)


ECGrav`ComputeMinusBetaTimesFreeEnergy[dat_Association]:=
(*Computes the value of -beta*free energy at the values of the inverse temperature beta
given as the keys in the association dat. 

Inputs is:,
1. dat - an association of inverse temperature and the corresponding 
   list of energies measured at that beta., 
   e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>

Returns an association with the inverse temperatures (betas) as keys and 
-(beta)(free energy) as the values.
*)
Module[{output,numTemps=Length[dat],betas=Keys[dat],
ntab=Table[Length[dat[[Key[i]]]],{i,Keys[dat]}],logntab,curFs,nextFs,iterNum,del,deltab},
logntab=Log[ntab*1.0];


Print[" Computing -beta*freeenergy for betas ",betas];
(*Print[" data ",dat];*)
(*Print[" numTemps ", numTemps, " ntable ",ntab," logntab ",logntab];*)

curFs=Table[1.0,{k,1,numTemps}];

nextFs=curFs;


(*Print["curFs ",curFs];*)

iterNum=0;
del=100.0;

(* Collecting delta For diagnostics *)
deltab=Reap[
While[del>10.0^(-5)&&iterNum<30000,
iterNum++;

nextFs=ParallelTable[
	ECGrav`LogSumExp[
		Flatten[
			Table[
				Table[
				-ECGrav`LogSumExp[
					Table[
						logntab[[j]]-curFs[[j]]+(betas[[k]]-betas[[j]])*(dat[[Key[betas[[i]]]]][[s]])
					,{j,1,numTemps}]]
				,{s,1,ntab[[i]]}]
			,{i,1,numTemps}]]
		]
,{k,1,numTemps},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}];

del=Sqrt[Sum[(1.0-curFs[[k]]/nextFs[[k]])^2,{k,1,numTemps}]];

(*If[Mod[iterNum,20]\[Equal]0,
Print["iterNum ",iterNum, " del ",del];
Print[" curFs ",curFs];
Print[" nextFs ",nextFs];
];*)


Sow[del];(* For diagnostics *)
If[Mod[iterNum,500]==0,PrintTemporary["iterNum ",iterNum, " del ",del]];
curFs=nextFs-(Mean[nextFs]);
];
];


(*Print[" curFs ",curFs,"deltab",deltab];*)
(*Print[" curFs ",curFs];*)

output=<|Table[betas[[i]]->curFs[[i]],{i,1,Length[betas]}]|>;

Remove[numTemps,betas,
ntab,logntab,curFs,nextFs,iterNum,del,deltab];

output

];

(* Catch-all Pattern *)
ECGrav`ComputeMinusBetaTimesFreeEnergy[args___]:=(Message[ECGrav`ComputeMinusBetaTimesFreeEnergy::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*NegativeBetaTimesFreeEnergy*)


(* Primary Pattern *)

ECGrav`NegativeBetaTimesFreeEnergy[bf_Real,minusBetaF_Association,energyMeasurements_Association]:=

(*Computes the value of -beta*free energy at the value of the inverse temperature bf. 
It requires two inputs:,
1. minusBetaF - an association of inverse temperatures and the corresponding value 
   of -beta*free energy, e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,
2. energyMeasurements - an association of inverse temperature and the corresponding 
   list of energies measured at that beta., 
   e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta values which are the keys for both associations have to be equal as 
sets! Also, the lengths of the lists of energy measurements have to be equal for all 
betas.
*)
With[{betas=Keys[minusBetaF],obsLength=Length[energyMeasurements[[1]]]},
(*Print[" betas ",betas];
Print[" minusBetaF ",minusBetaF];
Print[" measurements ",measurements];
Print[" obsLength ",obsLength];*)

	ECGrav`LogSumExp[
		Flatten[
			Table[
				Table[
					-ECGrav`LogSumExp[
						Table[
							Log[obsLength]-minusBetaF[[Key[betas[[j]]]]]+(bf-betas[[j]])*(energyMeasurements[[Key[betas[[i]]]]][[s]])
						,{j,1,Length[betas]}
						]
					]
				,{s,1,obsLength}]
			,{i,1,Length[betas]}]
		]
	]
];

(* Catch-all Pattern *)
ECGrav`NegativeBetaTimesFreeEnergy[args___]:=(Message[ECGrav`NegativeBetaTimesFreeEnergy::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*InternalEnergy*)


(* Primary Pattern *)
ECGrav`InternalEnergy[bf_Real,minusBetaF_Association,energyMeasurements_Association]:=
(*Computes the value of the interval energy at the value of the inverse temperature bf.
It requires two inputs:,
1. minusBetaF - an association of inverse temperatures and the corresponding value of -beta*free energy, e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>
2. energyMeasurements - an association of inverse temperature and the corresponding list of energies measured at that beta., e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta values which are the keys for both associations have to be equal as sets! Also, the lengths of the lists of energy measurements have to be equal for all betas.
*)
With[{betas=Keys[minusBetaF],obsLength=Length[energyMeasurements[[1]]]},
(Exp[-ECGrav`NegativeBetaTimesFreeEnergy[bf,minusBetaF,energyMeasurements]])*
Sum[
Sum[(energyMeasurements[[Key[betas[[i]]]]][[s]])/
(Sum[
obsLength*Exp[(bf-betas[[j]])*energyMeasurements[[Key[betas[[i]]]]][[s]]]*Exp[-minusBetaF[[Key[betas[[j]]]]]]
,{j,1,Length[betas]}])
,{s,1,obsLength}]
,{i,1,Length[betas]}]
];

(* Catch-all Pattern *)
ECGrav`InternalEnergy[args___]:=(Message[ECGrav`InternalEnergy::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*CvOverT*)


(* Primary Pattern *)
ECGrav`CvOverT[bf_Real,minusBetaF_Association,energyMeasurements_Association]:=

(*Cv/T*)

(*Computes the value of specific heat divided by temperature (Cv/T) at the value of the inverse temperature bf. Note, it is not computing Cv/T per site! 
It requires two inputs:,
1. minusBetaF - an association of inverse temperatures and the corresponding value of -beta*free energy, e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>
2. energyMeasurements - an association of inverse temperature and the corresponding list of energies measured at that beta., e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta values which are the keys for both associations have to be equal as sets! Also, the lengths of the lists of energy measurements have to be equal for all betas.
*)
With[{betas=Keys[minusBetaF],obsLength=Length[energyMeasurements[[1]]]},
bf^2(*/(vCount*vCount)**)(
Exp[-ECGrav`NegativeBetaTimesFreeEnergy[bf,minusBetaF,energyMeasurements]]*
Sum[
Sum[(energyMeasurements[[Key[betas[[i]]]]][[s]])^2/
(Sum[
obsLength*Exp[-minusBetaF[[Key[betas[[j]]]]]]*Exp[(bf-betas[[j]])*(energyMeasurements[[Key[betas[[i]]]]][[s]])]
,{j,1,Length[betas]}])
,{s,1,obsLength}]
,{i,1,Length[betas]}]-(ECGrav`InternalEnergy[bf,minusBetaF,energyMeasurements])^2
)
];

(* Catch-all Pattern *)
ECGrav`CvOverT[args___]:=(Message[ECGrav`CvOverT::argerr, args];
$Failed);


(* ::Chapter:: *)
(*Hamiltonians*)


(* ::Section::Closed:: *)
(*Graph Hamiltonians*)


(* ::Item::Closed:: *)
(*HIsing*)


(* Primary Pattern *)
ECGrav`HIsing[Am_List,J_Real,L_Real]:=
(*J/2 sum_{i!=j}A^2_{ij} + L/2 sum_{i!=j}A_{ij}*)
With[{Amsq=Am . Am, Fm=Table[1,{n,Length[Am]},{m,Length[Am]}]},
(J/2)*(Tr[Amsq . Fm]-Tr[Amsq])+(L/2)Tr[Am . Fm]
];

(* Catch-all Pattern *)
ECGrav`HIsing[args___]:=(Message[ECGrav`HIsing::argerr, args];
$Failed);


(* Primary Pattern *)
ECGrav`delHIsing[Am_List,J_,L_,a_Integer,b_Integer]:=
(* computes HIsing[Amnew] - HIsing[Am] where Amnew if found by toggling Am at 
row a and col b. HIsing = J/2 sum_{i!=j}A^2_{ij} + L/2 sum_{i!=j}A_{ij} *)
With[{n=Length[Am],togab=If[Am[[a,b]]==0,1,-1]},
togab*(J*(Total[Am[[a]]]+Total[Am[[b]]]-2Am[[a,b]])+L)
];

(* Catch-all Pattern *)
ECGrav`delHIsing[args___]:=(Message[ECGrav`delHIsing::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HWeightedFaceCounts*)


(* Primary Pattern *)
ECGrav`HWeightedFaceCounts[Am_List,J1_Real,J2_Real,J3_Real,J4_Real, J5_Real]:=
(*J1*(# of vertices) + J2*(number of edges)+... + J5*(number of 5-cliques).*)
Block[{n=Length[Am],num2cliques,num3cliques,num4cliques,num5cliques},
(*clqs=FindClique[g,\[Infinity],All];*)
num2cliques=Tr[Am . Am];
num3cliques=Tr[MatrixPower[Am,3]]/6;
(*Print["num3cliques ",num3cliques];*)
num4cliques=Sum[Product[Am[[a[[1]],a[[2]]]],{a,Subsets[i,{2}]}],{i,Subsets[Range[n],{4}]}];
(*Print["num4cliques ",num4cliques];*)
num5cliques=Sum[Product[Am[[a[[1]],a[[2]]]],{a,Subsets[i,{2}]}],{i,Subsets[Range[n],{5}]}];
J1*n+J2*num2cliques+J3*num3cliques+J4*num4cliques+J5*num5cliques
];

(* Catch-all Pattern *)
ECGrav`HWeightedFaceCounts[args___]:=(Message[ECGrav`HWeightedFaceCounts::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HEdgeDeg*)


(* Primary Pattern *)
ECGrav`HEdgeDeg[Am_List,J_Real,D1_Real,D2_Real]:=
(*J/(2N)*(Tr(A^4-A(2D1K + 2D2I)A)+D1^2*n(n-1)).*)
With[{n=Length[Am],Amsq=Am . Am},
(*Print["Lt "MatrixForm[2*D1*(Table[1,{i,n},{j,n}]-IdentityMatrix[n])+2*D2*IdentityMatrix[n]]];*)
(J/(2*n))*((Tr[Amsq . Amsq]
   -Tr[Amsq . (2*D1*(Table[1,{i,n},{j,n}]-IdentityMatrix[n])+2*D2*IdentityMatrix[n])]))
];

(* Catch-all Pattern *)
ECGrav`HEdgeDeg[args___]:=(Message[ECGrav`HEdgeDeg::argerr, args];
$Failed);


(* Primary Pattern *)
ECGrav`delHEdgeDeg[Am_List,J_Real, D1_Real, D2_Real,Amsq_List,a_Integer,b_Integer]:=
(*If you toggle edge {a,b} in Am to get Amnew, and Amsq = Am.Am, then this function 
gives the energy difference delta = HEdgeDeg[Amnew]-HEdgeDeg[Am].*)
With[{n=Length[Am],togab=If[Am[[a,b]]==0,1,-1]},
J*(togab-2*D2+2*(1-D1*togab)*(Amsq[[a,a]]+Amsq[[b,b]])+4*togab*(Amsq[[a]] . Am[[b]]-(D2-D1)Am[[a,b]]))/(n)
];

(* Overload Pattern *)
ECGrav`delHEdgeDeg[Am_List,J_Real, D1_Real, D2_Real,a_Integer,b_Integer]:=
(*If you toggle edge {a,b} in Am to get Amnew, and Amsq = Am.Am, then this function 
gives the energy difference delta = HEdgeDeg[Amnew]-HEdgeDeg[Am].*)
With[{Amsq=Am . Am,n=Length[Am],togab=If[Am[[a,b]]==0,1,-1]},
J*(togab-2*D2+2*(1-D1*togab)*(Amsq[[a,a]]+Amsq[[b,b]])+4*togab*(Amsq[[a]] . Am[[b]]-(D2-D1)Am[[a,b]]))/(n)
];

(* Catch-all Pattern *)
ECGrav`delHEdgeDeg[args___]:=(Message[ECGrav`delHEdgeDeg::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HLaplacian*)


(* Primary Pattern *)
ECGrav`HLaplacian[Amat_List,J_Real]:=
(*Sum_{ij} deg(i)L_{ij}deg(j), where L is the graph Laplacian, L = Diag[deg] - Amat*)
With[{degv={Total[#]}&/@Amat,Lmat=DiagonalMatrix[Total/@Amat]-Amat},
(*Print[" degv ",degv," Lmat ",Lmat];*)
J*(Transpose[degv] . Lmat . degv)[[1,1]]];

(* Catch-all Pattern *)
ECGrav`HLaplacian[args___]:=(Message[ECGrav`HLaplacian::argerr, args];
$Failed);


(* ::Chapter:: *)
(*Ground State Searches*)


(* ::Section::Closed:: *)
(*Gradient Descent (GD)*)


(* ::Item::Closed:: *)
(*GradDescent*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`GradDescent[seedAmat_List,delH_,cutoff_Integer]:=
(*This program runs the gradient descent algorithm.,
Inputs are:, 
1. seedAmat - a seed graph as an adjacency matrix, 
2. delH - a formula for delta E (when one edge is flipped),  
3. cutoff -  number of sweeps,

Outputs an adjacency matrix. *)

Block[{ nn=Length[seedAmat],Amcur,edgeList,deltaEtable,updateAmsq,step,stepagain,numsteps,printCase},

Amcur=seedAmat;
edgeList=Subsets[Range[nn],{2}];
deltaEtable=<|Table[i->0,{i,edgeList}]|>;

step[]:=Block[{negatives,flipSpin},

Do[

deltaEtable[[Key[i]]]=delH[Amcur,i[[1]],i[[2]]];
(*Print["In step, calculating deltaE for i ",i, " deltaE is ",deltaEtable[i]];*)
,{i,edgeList}];

(*Print[" deltaEtable ",deltaEtable];*)

negatives=Select[deltaEtable,#<=0&];
(*Print["negatives ", negatives];*)

If[Length[negatives]>0,
flipSpin=Keys[TakeSmallest[negatives,1]][[1]];

(*Print["flipping spin ", flipSpin];*)

Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2],
Return[0](* unable to find a flip that lowers the energy*)
];

1 (*success*)
];

printCase=Floor[(cutoff*1.0)/5.0];

numsteps=0;
stepagain=1;

While[stepagain==1&&numsteps<cutoff,
If[Mod[numsteps,printCase]==0,Print["step number ",numsteps]];

(*Print["before step, graph", AdjacencyGraph[Amcur,ImageSize\[Rule]Tiny]];
Print["before step, eDifftable ", MatrixPlot[to2DArray[nn,deltaEtable]]];
*)

stepagain*=step[];

(*Print["After step, graph", AdjacencyGraph[Amcur,ImageSize\[Rule]Tiny]];
Print["After step, eDifftable ", MatrixPlot[to2DArray[nn,deltaEtable]]];*)

numsteps++;
];

Amcur
];

(* Catch-all Pattern *)
ECGrav`GradDescent[args___]:=(Message[ECGrav`GradDescent::argerr, args];
$Failed);


(* ::Section::Closed:: *)
(*Stochastic Gradient Descent (SGD)*)


(* ::Item::Closed:: *)
(*SGradDescent*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`SGradDescent[seedAmat_List, hamiltonian_,delH_,beta_Real,NN_Integer]:=

(*This program runs stochastic gradient descent with softmax at parameter beta.,
Inputs are:, 
1. seedAmat - a seed graph as an adjacency matrix, 
2. hamiltonian - a hamiltonian, 
3. delH - a formula for delta E (when one edge is flipped to expedite computation),  
4. beta - inverse temperature,
5. NN -  number of sweeps,

Outputs an association with three elements,  
1. the minimum energy visited throughout the search,
2. states with that minimum energy. If multiple states have degenerate minimum energy, 
they will all be included.,
3. The last state visited. *)


Block[{ nn=Length[seedAmat],Amcur,edgeList,deltaEtable,weightsTable,curE,minE, minStates,computeWeights,step,numsteps,printCase},


Amcur=seedAmat;
edgeList=Subsets[Range[nn],{2}];
deltaEtable=<|Table[i->0,{i,edgeList}]|>;
weightsTable=deltaEtable;

(*Print["deltaEtable ",deltaEtable," weightsTable ",weightsTable];*)

curE = hamiltonian[Amcur];
minE = curE;
minStates={Amcur};


computeWeights[bta_Real]:=Module[{},
Do[weightsTable[i]=Exp[-bta*deltaEtable[i]],{i,edgeList}];
weightsTable=weightsTable/Total[weightsTable];
];


step[]:=Block[{negatives,flipSpin},


Do[

deltaEtable[[Key[i]]]=delH[Amcur,i[[1]],i[[2]]];
(*Print["In step, calculating deltaE for i ",i, " deltaE is ",deltaEtable[i]];*)
,{i,edgeList}];

computeWeights[beta];
flipSpin=RandomChoice[Values[weightsTable]->Keys[weightsTable]];

(*Print["Flipping ",flipSpin, " which gives delEval ",deltaEtable[flipSpin]];*)


Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2];


curE=curE+deltaEtable[flipSpin];
If[curE<minE,
(*Print["updating minE and minAm b/c: curE = ", curE, " minE = ",minE];*)
(*Print["direct calculations: H[Amcur]  = ", hamiltonian[Amcur,hparams], " H[minAm] = ",hamiltonian[minAm,hparams] ];*)
minE=curE;
minStates={Amcur},
If[curE==minE,
minStates=Join[minStates,{Amcur}];
];
];

1 (*success*)
];

printCase=Floor[(NN*1.0)/5.0];

numsteps=0;

While[numsteps<NN,
If[Mod[numsteps,printCase]==0,Print["step number ",numsteps]];


step[];

numsteps++;
];

<|"minE"->minE,"minEstates"->minStates,"LastState"->Amcur|>
];

(* Overload Pattern *)
ECGrav`SGradDescent[seedAmat_List, hamiltonian_,beta_Real,NN_Integer]:=

(*This program runs stochastic gradient descent with softmax at parameter beta. It is 
an overload that doesn't require the formula for delH (E(newgraph) - E(curgraph) when 
an edge is toggled),
Inputs are:, 
1. seedAmat - a seed graph as an adjacency matrix, 
2. hamiltonian - a hamiltonian,  
3. beta - inverse temperature,
4. NN -  number of sweeps,

Outputs an association with three elements,  
1. the minimum energy visited throughout the search,
2. states with that minimum energy. If multiple states have degenerate minimum energy, 
they will all be included.,
3. The last state visited. *)


Block[{ nn=Length[seedAmat],Amcur,Amnew, edgeList,deltaEtable,weightsTable,
curE,minE, minStates,computeWeights,step,numsteps,printCase},


Amcur=seedAmat;
edgeList=Subsets[Range[nn],{2}];
deltaEtable=<|Table[i->0,{i,edgeList}]|>;
weightsTable=deltaEtable;

(*Print["deltaEtable ",deltaEtable," weightsTable ",weightsTable];*)

curE = hamiltonian[Amcur];
minE = curE;
minStates={Amcur};
computeWeights[bta_Real]:=Module[{},
Do[weightsTable[i]=Exp[-bta*deltaEtable[i]],{i,edgeList}];
weightsTable=weightsTable/Total[weightsTable];
];


step[]:=Block[{negatives,flipSpin},


Do[

Amnew=Amcur;
Amnew[[i[[1]],i[[2]]]]=Amnew[[i[[2]],i[[1]]]]=Mod[Amnew[[i[[1]],i[[2]]]]+1,2];

deltaEtable[[Key[i]]]= hamiltonian[Amnew]-curE;
(*Print["In step, calculating deltaE for i ",i, " deltaE is ",deltaEtable[i]];*)
,{i,edgeList}];

computeWeights[beta];
flipSpin=RandomChoice[Values[weightsTable]->Keys[weightsTable]];

(*Print["Flipping ",flipSpin, " which gives delEval ",deltaEtable[flipSpin]];*)


Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2];


curE=curE+deltaEtable[flipSpin];
If[curE<minE,
(*Print["updating minE and minAm b/c: curE = ", curE, " minE = ",minE];*)
(*Print["direct calculations: H[Amcur]  = ", hamiltonian[Amcur,hparams], " H[minAm] = ",hamiltonian[minAm,hparams] ];*)
minE=curE;
minStates={Amcur},
If[curE==minE,
minStates=Join[minStates,{Amcur}];
];
];

1 (*success*)
];

printCase=Floor[(NN*1.0)/5.0];

numsteps=0;

While[numsteps<NN,
If[Mod[numsteps,printCase]==0,Print["step number ",numsteps]];


step[];

numsteps++;
];

<|"minE"->minE,"minEstates"->minStates,"LastState"->Amcur|>
];

(* Catch-all Pattern *)
ECGrav`SGradDescent[args___]:=(Message[ECGrav`SGradDescent::argerr, args];
$Failed);


(* ::Section::Closed:: *)
(*Simulated Annealing (SA)*)


(* ::Item::Closed:: *)
(*SimulatedAnnealing*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`SimulatedAnnealing[seedAmat_List, hamiltonian_,delH_,betai_Real,betaf_Real,roundLength_Integer,NN_Integer]:=                     
(*************************************)
(***  Last updated on: 01/15/2025  ***)
(*************************************)
(*Version: 1*)
(*Notes: ,
1. 01/07/23 update updated the sweep.,
2. 02/06/23 update added storage of first excited states.,
3. 03/10/23 added precisions.,
4. 01/15/2025 added overload and fixed a bug.,
*)
(*This program runs the simulated annealing ground state search from initial to final inverse temperatures given by betai and betaf.,
Inputs are:, 
1. seedAmat - a seed graph as an adjacency matrix, 
2. hamiltonian - a hamiltonian, 
3. delH - a formula for delta E (when one edge is flipped to expedite computation),  
4. betai - starting low inverse temperature,
5. betaf - the final high inverse temperature,
6. roundLength - how many rounds to do MC steps at each temperature
7. NN -  number of sweeps,

Outputs an association with five elements,
1. the minimum energy visited throughout the search,
2. the second minimum energy,
3. the "ground states", the states with that minimum energy. If multiple states have degenerate minimum energy, they will all be included.,
4. the "first excited states", i.e., the states correspondign to the second minimum energy,
5. The last state visited. 


*)

Module[{ nn=Length[seedAmat],result,maxNumOfSavedStates=400,(*precision=100,*)precision=MachinePrecision,Amcur,edgeList,minE,excitedEnergy,curE,minStates,excitedStates,beta,rate,step,numsteps,printCase},


Amcur=seedAmat;
edgeList=Subsets[Range[nn],{2}];

curE = SetPrecision[hamiltonian[Amcur],precision];
minE = curE;
excitedEnergy=N[minE+SetPrecision[1000.0,precision],precision];
minStates={Amcur};
excitedStates={Amcur};

beta = betai;
rate = SetPrecision[Exp[Log[betai/betaf]*roundLength/NN],precision];

Print["rate ",rate]; 



(*Print["state ", state];
Print["finalstate ", finalState];*)

step[bta_]:=
Block[{delE,expdelE,flipSpin,accept},

	flipSpin=RandomChoice[edgeList];
	delE = N[delH[Amcur,flipSpin[[1]],flipSpin[[2]]],precision];

	accept = 0;
	If[delE<=0,accept = 1,
		expdelE = N[Exp[-SetPrecision[bta,precision]*delE],precision];
		If[RandomReal[]<=expdelE,accept =1];
	];

If[accept ==0, Return[0]];

If[accept==1,

(*Print["flipSpin ",flipSpin," Amcur ",Amcur//MatrixForm];*)

Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2];

(*Print["flipSpin ",flipSpin," Amnext ",Amcur//MatrixForm];*)

(*Print[" Amnext ",Amcur//MatrixForm];*)

curE=curE+delE;

(*Print["Before Which:, curE ",curE," h[curstate] ",hamiltonian[Amcur]];
Print["  minE ",minE, " h[gstate] ",hamiltonian[minStates[[1]]]];
Print[ "  excitedEnergy ",excitedEnergy," h[excitedStates] ",hamiltonian[excitedStates[[1]]]];*)

Which[
	curE<minE,
		{excitedEnergy=minE; excitedStates=minStates; minE=curE; minStates={Amcur};},
	curE==minE&&Length[minStates]<maxNumOfSavedStates,
		{minStates=Join[minStates,{Amcur}];},
	minE<curE<excitedEnergy,
		{excitedEnergy=curE; excitedStates={Amcur};},
	curE==excitedEnergy&&Length[excitedStates]<maxNumOfSavedStates,
		{excitedStates=Join[excitedStates,{Amcur}];}
];

(*Print["After Which:, curE ",curE," h[curstate] ",hamiltonian[Amcur]];
Print["  minE ",minE, " h[gstate] ",hamiltonian[minStates[[1]]]];
Print[ "  excitedEnergy ",excitedEnergy," h[excitedStates] ",hamiltonian[excitedStates[[1]]]];*)

];

1 (*success*)
];

printCase=Floor[(NN*1.0)/5.0];


numsteps=0;

While[numsteps<NN,

If[Mod[numsteps,printCase]==0,Print["sweep number ",numsteps]];

(*Print["before step , deltaEtable ", deltaEtable];
Print["before step , state ", state];*)

Do[step[beta],{roundLength*nn*nn}];
(*Do[step[beta],{1}];*)
beta = N[betai*rate^(-Floor[numsteps/roundLength]),precision];

(*Print["after step , deltaEtable ", deltaEtable];
Print["after step , state ", state];*)

numsteps++;
];

result=<|"minE"->minE,"excitedEnergy"->excitedEnergy,"minEstates"->minStates,"excitedStates"->excitedStates,"LastState"->Amcur|>;

Remove[ nn,maxNumOfSavedStates,Amcur,edgeList,minE,excitedEnergy,curE,minStates,excitedStates,beta,rate,step,numsteps,printCase];

result
];

(* Overload Pattern *)
ECGrav`SimulatedAnnealing[seedAmat_List, hamiltonian_,betai_Real,betaf_Real,roundLength_Integer,NN_Integer]:=                     

(*************************************)
(***  Last updated on: 01/15/2025  ***)
(*************************************)
(*Notes: an overload for the case when delh is not given,

*)
(*This program runs the simulated annealing ground state search from initial to final inverse temperatures given by betai and betaf.,
Inputs are:, 
1. seedAmat - a seed graph as an adjacency matrix, 
2. hamiltonian - a hamiltonian,   
3. betai - starting low inverse temperature,
4. betaf - the final high inverse temperature,
5. roundLength - Integer, related to how many rounds of MC steps to do at each 
   temperaturenumber before advancing to the next temperature,
6. NN -  number of sweeps,

Outputs an association with five elements,  list with two elements associations,
1. the minimum energy visited throughout the search,
2. the second minimum energy,
3. the "ground states", the states with that minimum energy. If multiple states have degenerate minimum energy, they will all be included.,
4. the "first excited states", i.e., the states correspondign to the second minimum energy,
5. The last state visited. 


*)


Module[{ nn=Length[seedAmat],result,maxNumOfSavedStates=400,(*precision=100,*)precision=MachinePrecision,Amcur,Amnext,edgeList,minE,excitedEnergy,curE,minStates,excitedStates,beta,rate,step,numsteps,printCase},


Amcur=seedAmat;
Amnext=seedAmat;
edgeList=Subsets[Range[nn],{2}];

curE = SetPrecision[hamiltonian[Amcur],precision];
minE = curE;
excitedEnergy=N[minE+SetPrecision[1000.0,precision],precision];
minStates={Amcur};
excitedStates={Amcur};

beta = betai;
rate = SetPrecision[Exp[Log[betai/betaf]*roundLength/NN],precision];

Print["rate ",rate]; 



(*Print["state ", state];
Print["finalstate ", finalState];*)

step[bta_Real]:=
Block[{delE,expdelE,flipSpin,accept},

flipSpin=RandomChoice[edgeList];
Amnext[[flipSpin[[1]],flipSpin[[2]]]]=Amnext[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amnext[[flipSpin[[1]],flipSpin[[2]]]]+1,2];

(*Print["flipSpin ",flipSpin," Amcur ",Amcur//MatrixForm," Amnext ",Amnext//MatrixForm];*)

delE = N[hamiltonian[Amnext],precision]-curE;
(*delE = N[hamiltonian[Amnext]-hamiltonian[Amcur],precision];*)

(*Print[" delE ",delE," h[Amcur] ",N[hamiltonian[Amcur],precision], " h[Amnext] ",
	N[hamiltonian[Amnext],precision], " h[Amnext] - h[Amcur] ",N[hamiltonian[Amnext]-hamiltonian[Amcur],precision]];*)
	
accept = 0;
If[delE<=0,accept = 1,
	expdelE = N[Exp[-SetPrecision[bta,precision]*delE],precision];
	If[RandomReal[]<=expdelE,accept =1];
];

If[accept==0,Amnext=Amcur;Return[0]];

If[accept==1,

Amcur=Amnext;

curE=curE+delE;

(*Print["Before Which:, curE ",curE," h[curstate] ",hamiltonian[Amcur]];
Print["  minE ",minE, " h[gstate] ",hamiltonian[minStates[[1]]]];
Print[ "  excitedEnergy ",excitedEnergy," h[excitedStates] ",hamiltonian[excitedStates[[1]]]];*)

Which[
	curE<minE,
		{excitedEnergy=minE;
		excitedStates=minStates;
		minE=curE;
		minStates={Amcur};
		},
	curE==minE&&Length[minStates]<maxNumOfSavedStates,
		{minStates=Join[minStates,{Amcur}];},
	minE<curE<excitedEnergy,
		{excitedEnergy=curE;
		excitedStates={Amcur};},
	curE==excitedEnergy&&Length[excitedStates]<maxNumOfSavedStates,
		{excitedStates=Join[excitedStates,{Amcur}];}
]

(*Print["After Which:, curE ",curE," h[Amcur] ",hamiltonian[Amcur]];
Print["  minE ",minE, " h[gstate] ",hamiltonian[minStates[[1]]]];
Print[ "  excitedEnergy ",excitedEnergy," h[excitedStates] ",hamiltonian[excitedStates[[1]]]];*)

];

1 (*success*)
];

printCase=Floor[(NN*1.0)/5.0];


numsteps=0;

While[numsteps<NN,

If[Mod[numsteps,printCase]==0,Print["sweep number ",numsteps]];

(*Print["before step , deltaEtable ", deltaEtable];
Print["before step , state ", state];*)

Do[step[beta],{roundLength*nn*nn}];
(*Do[step[beta],{1}];*)
beta = N[betai*rate^(-Floor[numsteps/roundLength]),precision];

(*Print["after step , deltaEtable ", deltaEtable];
Print["after step , state ", state];*)

numsteps++;
];

result=<|"minE"->minE,"excitedEnergy"->excitedEnergy,"minEstates"->minStates,"excitedStates"->excitedStates,"LastState"->Amcur|>;

Remove[ nn,maxNumOfSavedStates,Amcur,edgeList,minE,excitedEnergy,curE,minStates,excitedStates,beta,rate,step,numsteps,printCase];

result
];

(* Catch-all Pattern *)
ECGrav`SimulatedAnnealing[args___]:=(Message[ECGrav`SimulatedAnnealing::argerr, args];
$Failed);



(* ::Chapter::Closed:: *)
(*Metropolis MC*)


(* ::Section:: *)
(*Single Flip Metropolis*)


(* ::Item::Closed:: *)
(*GraphMetropolis*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`GraphMetropolis[seedAmat_List, beta_Real, hamiltonian_,observables_,maxSweep_Integer]:=

(*(*****************************)
(* Last Updated: 08/26/2025  *)
(*****************************)*)
(*Version : 1 *)
(*Notes: *)
(*Runs Metropolic MC simulation on the seed graph with Adjacency matrix seedAmat at inverse temperature beta.,
Inputs are:,
1. seedAmat= seed adjacency matrix,
2. beta = inverse temperature,
3. hamiltonian_[hparams__]= the hamiltonian which can depend on arbitrary parameters,
4. observables = a list of cuntions on graph adjacency matrices e.g. magnetization, average degree etc.,
5. maxSweep= the total number of sweeps where one sweep is N(N-1)/2 MC flip attempts.,
It outputs the last graph in the simulation and a list of the values of the energy, magnetization, and observables.*)
Module[{N=Length[seedAmat],maxEdgeCount,expDEVals,MCStep,MCSweep,numsweeps,Amat,AmatNew,printVal,result},
maxEdgeCount=N (N-1)/2;
printVal=Floor[maxSweep/5];

Amat = seedAmat;
AmatNew = Amat;
Print["seed graph",AdjacencyGraph[Amat]];

MCStep[]:=Module[{l,lp,m,i,j,deltaE,expBetaDE,accept},(*one MC flip attempt*)
Do[
AmatNew = Amat;
l=RandomInteger[{1,maxEdgeCount}];

lp=maxEdgeCount+1-l;
m=0;
While[lp>0,
m++;
lp-=N-m;
];
j=N+1-m;
i=l-(j-1)(j-2)/2;

AmatNew[[i,j]]=Mod[Amat[[i,j]]+1,2];AmatNew[[j,i]]=AmatNew[[i,j]];
deltaE=hamiltonian[AmatNew]-hamiltonian[Amat];
accept=False;
If[
deltaE<=0,accept=True,
expBetaDE=Exp[-beta*deltaE];
If[RandomReal[]<expBetaDE,accept=True]
];
If[accept,
Amat = AmatNew;
],
(*Print["after ",Amat//MatrixForm];*)
maxEdgeCount
];
];

MCSweep[]:=Module[{},
Do[MCStep[],{maxEdgeCount}]
];

numsweeps=0;

result=Reap[
Do[MCSweep[];
numsweeps++;
If[Mod[numsweeps,printVal]==0,Print["Now at sweep number ",numsweeps]];


Sow[Flatten[{numsweeps,hamiltonian[Amat],Through[observables [Amat]]}]]
,maxSweep]
][[2,1]];

{Amat,result}

];

(* Catch-all Pattern *)
ECGrav`GraphMetropolis[args___]:=(Message[ECGrav`GraphMetropolis::argerr, args];
$Failed);


(* ::Chapter:: *)
(*Parallel Tempering*)


(* ::Section:: *)
(*Parallel Tempering*)


(* ::Subsection::Closed:: *)
(*GraphSweepReplica*)


(* ::Item::Closed:: *)
(*GraphMetropolis Primary*)


ECGrav`GraphSweepReplica[seedGraph_List,beta_Real,hamiltonian_,delH_,NN_Integer,minEToBeat_Real,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 1/30/2025  ***)
(*************************************)

(*Notes :,
*1. Memory leak in the parallelization is fixed. ,
*2. This version can implement selection probability to make graphs unlabaled, but GraphAutomorphismGroup function itself has memory leak.,
*3. 1/30/2025 Updated enabeled one to chose whether the weep is done with labeled or unlabeled graphs *)

(* This function performs NN sweeps on a seed graph state.,
Inputs are:, 
1. a seed graph as an adjacency matrix, 
2. inverse temperature beta, 
3. a hamiltonian, 
4. a formula for delta E (when one edge is flipped to expedite computation),  
5. number of sweeps NN,
6. a value for energy such that the lowest energy states with energy lower than minEToBeat will be saved (for ground state search).,
7. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.   

Outputs a list with two associations,
1. the minimum energy visited throughout the sweep states with that energy. If multiple states have degenerate minimum energy, they will all be included. It saves only non-isomorphic graphs. ,
2. The second association is the temperature, final graph, energy, and magnetization at the end of the sweeps.
*)
Module[{result,vCount=Length[seedGraph],minE, minStates,maxGStateCount,expDelETable,data,step},
minE = minEToBeat;
minStates={};
maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)
expDelETable=<||>;
data=<|"graph"->seedGraph,"energy"->hamiltonian[seedGraph], "mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>;

step[]:=(*Performs one spin flip step*)
Block[{
curE=data[[Key["energy"]]],
curM=data[[Key["mag"]]],
row=RandomInteger[{1,vCount-1}],
col,
newAmat,
selectionProb,
delE,
expdelE,
accept
},

col=RandomInteger[{row+1,vCount}];

delE=delH[data[[Key["graph"]]],row,col];

(*Print["In step, calculated delta E for edge ", {row,col}, " delE ",delE];*)

(*Print["Before spin flip of edge ",{row,col}, " state is ",data[[Key["graph"]]]//MatrixForm];*)


newAmat=data[[Key["graph"]]];
newAmat[[row,col]]=newAmat[[col,row]]=Mod[newAmat[[row,col]]+1,2];

Which[UnlabeledVerticesYes==0,selectionProb=1,UnlabeledVerticesYes==1,
selectionProb=(*Automorphism group order of new graph over auromorphism group order of the old*)
GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[newAmat]]]/GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[data[[Key["graph"]]]]]],True,Print[" Error, in GraphSweepReplica, UnlabeledVerticesYes = ",UnlabeledVerticesYes," is invalid input. It has to be 0 for labeled or 1 for unlabeled "];Abort[];
];

(*Print["current Amat is ",data[[Key["graph"]]]//MatrixForm];
Print["row col :",row," " ,col," newAmat is ",newAmat//MatrixForm];
Print[" selectionProb ",selectionProb];*)

accept = 0;
expdelE = Lookup[expDelETable,delE];
If[MissingQ[expdelE],
expdelE=Exp[-delE*beta];expDelETable[[Key[delE]]]=expdelE];

(*If[delE<=0,accept = 1,

If[RandomReal[]<selectionProb*expdelE,accept =1],
accept = 0
]; (*This seems to work well and ca reduce the number of automorphism group order calculationsthough it doesn't have a sound logical basis*)*)


If[selectionProb*expdelE>=1,accept = 1,

If[RandomReal[]<selectionProb*expdelE,accept =1],
accept = 0
];

If[accept==1,

(*Print["In step, accepting spin flip at edge ",{row,col}, " of state ",data[[Key["graph"]]]//MatrixForm];*)


data[[Key["graph"]]]=newAmat;
data[[Key["energy"]]]+= delE;
data[[Key["mag"]]]+= (2.0/(vCount(vCount-1)))*(2*data[[Key["graph"]]][[row,col]]-1);
(*Note the spin has already been flipped.*)

(*Note, below since flip is accepted, the energy of the new state is curE + delE. We cannot reassign curE because it is passed by value i.e., is is immutable, hence the use of curE + delE *)

If[curE+delE<minE,minE=curE+delE;minStates={data[[Key["graph"]]]},
If[curE+delE==minE&&Length[minStates]<maxGStateCount,minStates=DeleteDuplicates[Union[minStates,{data[[Key["graph"]]]}]]]];

];
];

Do[step[];

(*Print["After step, graph is ",data[[Key["graph"]]]//MatrixForm];*)
,{i,NN*vCount (vCount-1)/2}];

(*Print["expDelEtable ",expDelETable];*)

result={<|"minEnergy"->minE,"minEstates"->minStates|>,data};

Remove[vCount,minE, minStates,maxGStateCount,expDelETable,data,step];


result

];


(* ::Item::Closed:: *)
(*GraphMetropolis Overload *)


ECGrav`GraphSweepReplica[seedGraph_List,beta_Real,hamiltonian_,NN_Integer,minEToBeat_Real,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 10/10/2025  ***)
(*************************************)

(*Notes :,
*1. Memory leak in the parallelization is fixed. ,
*2. This version can implement selection probability to make graphs unlabaled, but GraphAutomorphismGroup function itself has memory leak.,
*3. 1/30/2025 Updated enabeled one to chose whether the weep is done with labeled or unlabeled graphs 
*4. 10/10/2025 Updated created the overload to enable delH to not be specified
s*)

(* This function performs NN sweeps on a seed graph state.,
Inputs are:, 
1. a seed graph as an adjacency matrix, 
2. inverse temperature beta, 
3. a hamiltonian,  
4. number of sweeps NN,
5. a value for energy such that the lowest energy states with energy lower than minEToBeat will be saved (for ground state search).,
6. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.   

Outputs a list with two associations,
1. the minimum energy visited throughout the sweep states with that energy. If multiple states have degenerate minimum energy, they will all be included. It saves only non-isomorphic graphs. ,
2. The second association is the temperature, final graph, energy, and magnetization at the end of the sweeps.
*)
Module[{result,vCount=Length[seedGraph],minE, minStates,maxGStateCount,expDelETable,data,step},
minE = minEToBeat;
minStates={};
maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)
expDelETable=<||>;
data=<|"graph"->seedGraph,"energy"->hamiltonian[seedGraph], "mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>;

step[]:=(*Performs one spin flip step*)
Block[{
	curE=data[[Key["energy"]]],
	curM=data[[Key["mag"]]],
	row=RandomInteger[{1,vCount-1}],
	col,
	newAmat,
	selectionProb,
	delE,
	expdelE,
	accept},

col=RandomInteger[{row+1,vCount}];

newAmat=data[[Key["graph"]]];
newAmat[[row,col]]=newAmat[[col,row]]=Mod[newAmat[[row,col]]+1,2];

delE=hamiltonian[newAmat]-curE;

(*Print["In step, calculated delta E for edge ", {row,col}, " delE ",delE];*)

(*Print["Before spin flip of edge ",{row,col}, " state is ",data[[Key["graph"]]]//MatrixForm];*)

Which[UnlabeledVerticesYes==0,selectionProb=1,UnlabeledVerticesYes==1,
selectionProb=(*Automorphism group order of new graph over auromorphism group order of the old*)
GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[newAmat]]]/GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[data[[Key["graph"]]]]]],True,Print[" Error, in GraphSweepReplica, UnlabeledVerticesYes = ",UnlabeledVerticesYes," is invalid input. It has to be 0 for labeled or 1 for unlabeled "];Abort[];
];

(*Print["current Amat is ",data[[Key["graph"]]]//MatrixForm];
Print["row col :",row," " ,col," newAmat is ",newAmat//MatrixForm];
Print[" selectionProb ",selectionProb];*)

accept = 0;
expdelE = Lookup[expDelETable,delE];
If[MissingQ[expdelE],
expdelE=Exp[-delE*beta];expDelETable[[Key[delE]]]=expdelE];


If[selectionProb*expdelE>=1,accept = 1,

	If[RandomReal[]<selectionProb*expdelE,accept =1],
	accept = 0
];

If[accept==1,

(*Print["In step, accepting spin flip at edge ",{row,col}, " of state ",data[[Key["graph"]]]//MatrixForm];*)


data[[Key["graph"]]]=newAmat;
data[[Key["energy"]]]+= delE;
data[[Key["mag"]]]+= (2.0/(vCount(vCount-1)))*(2*data[[Key["graph"]]][[row,col]]-1);
(*Note the spin has already been flipped.*)

(*Note, below since flip is accepted, the energy of the new state is curE + delE. We cannot reassign curE because it is passed by value i.e., is is immutable, hence the use of curE + delE *)

If[curE+delE<minE,minE=curE+delE;minStates={data[[Key["graph"]]]},
If[curE+delE==minE&&Length[minStates]<maxGStateCount,minStates=DeleteDuplicates[Union[minStates,{data[[Key["graph"]]]}]]]];

];
];

Do[step[];

(*Print["After step, graph is ",data[[Key["graph"]]]//MatrixForm];*)
,{i,NN*vCount (vCount-1)/2}];

(*Print["expDelEtable ",expDelETable];*)

result={<|"minEnergy"->minE,"minEstates"->minStates|>,data};

Remove[vCount,minE, minStates,maxGStateCount,expDelETable,data,step];


result

];


(* ::Item::Closed:: *)
(*GraphMetropolis Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphSweepReplica[args___]:=(Message[ECGrav`GraphSweepReplica::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*GraphEquilibriate*)


(* ::Item::Closed:: *)
(*GraphEquilibriate Primary*)


ECGrav`GraphEquilibriate[seedGraph_List, beta_Real, hamiltonian_,delH_,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 1/30/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/22/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs*)

(*Equilibriate an input graph configuration to the temperature beta.,
Inputs are:,
1. seedGraph - adjacency matrix of the seed graph,
2. beta = inverse temperature,
3. a hamiltonian, 
4. a formula for delta E (when one edge is flipped to expedite computation),  
5. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Depends on the function GraphSweepReplicas., 

Outputs a list with two associations:, the first one is the lowest energy found  throughout the run together with all states with that lowest energy;, 
The second is the final equilibriated state
*)
Module[{result,vCount=Length[seedGraph],data,maxGStateCount,sweepOutput,Entable,outWinLength,inWinLength,AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime=20000,maxNumSweeps=25000},

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)
data=<|"minEnergy" ->hamiltonian[seedGraph],"minEstates"->{seedGraph},
"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph],"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>,"empty"-><|"graph"->Table[0,{i,vCount},{j,vCount}],"energy" ->0.0,"mag"->0.0|>,
"random"-><|"graph"->With[{rm=RandomInteger[1,{vCount,vCount}]},Mod[rm+Transpose[rm],2]],"energy" ->0.0,"mag"->0.0|>|>;

data[[Key["empty"]]][[Key["energy"]]]=hamiltonian[data[[Key["empty"]]][[Key["graph"]]]];
data[[Key["random"]]][[Key["energy"]]]=hamiltonian[data[[Key["random"]]][[Key["graph"]]]];
data[[Key["random"]]][[Key["mag"]]]=Total[Flatten[data[[Key["random"]]][[Key["graph"]]]]]*1.0/(vCount(vCount-1));


Print["Equilibriating at beta ",beta];
(*Print[" starting with data ",data];*)

numsweeps=0;
outWinLength =500 ;(* length of a table to store running energy values to test equilibriation*)

inWinLength=20;(* A segment to be averaged over at the beginning and the end of the table*)

Entable=Table[0.0,{m,3},{n,outWinLength}];(*A table to store energy values of the seed, random and empty tracks for tunning calculation of tests of equilibriation*)

(*Main calculation*)

AllEnMat=
(*AllEnMat is for storing all energy values for diagnostic*)
Reap[
While[numsweeps<maxNumSweeps,

numsweeps++;

Sow[{data[[Key["state"]]][[Key["energy"]]],data[[Key["empty"]]][[Key["energy"]]],data[[Key["random"]]][[Key["energy"]]]}];

Do[

(*Print["sweeping at ",replicaName];*)
(*Print["minEnergy ",data[[Key["minEnergy"]]]];*)

sweepOutput=ECGrav`GraphSweepReplica[data[[Key[replicaName]]][[Key["graph"]]],beta,hamiltonian,delH,1,data[[Key["minEnergy"]]],UnlabeledVerticesYes];

(*Print["In GraphEquilibriate, sweep output ",sweepOutput];*)

data[[Key[replicaName]]]=sweepOutput[[2]];

(*Update minimum energy states after the sweep*)
If[sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
data[[Key["minEnergy"]]]=sweepOutput[[1,Key["minEnergy"]]];
data[[Key["minEstates"]]]=sweepOutput[[1,Key["minEstates"]]],
If[sweepOutput[[1,Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<=maxGStateCount,
data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1,Key["minEstates"]]]]];
];
,{replicaName,{"state","empty","random"}}];

(*Print["after one sweep, replicas are ",data];*)


Entable[[All,1]]={data[[Key["state"]]][[Key["energy"]]],data[[Key["empty"]]][[Key["energy"]]],data[[Key["random"]]][[Key["energy"]]]};


(*Print["Entable ",Entable//MatrixForm];*)


If[numsweeps>outWinLength,

sqMeanEMat=Table[(Mean[Entable[[i]][[1;;inWinLength]]]-Mean[Entable[[j]][[(outWinLength-inWinLength);;outWinLength]]])^2,{i,3},{j,3}];

(*sqMeanEMat is a three by three matrix of the squared means of difference in energy within and across the tracks at the beginning and end of outWinLength.
At equilibrium, these 9 mumbers should be randomly distributed with mean of 0. Their fluctuation from 0 should be within the variance of the newer (late time) variance in energy for the equilibriation to exit. sqMeanPairwiseDiff is the mean of these 9 numbers*)

sqMeanPairwiseDiff=Mean[Flatten[sqMeanEMat]];

meanLateVar=Mean[Table[Variance[Entable[[i]][[1;;inWinLength]]],{i,3}]];(*Mean variance of the newer energies*)


(*Print["sqMeanEMat ",sqMeanEMat//MatrixForm," sqMeanPairwiseDiff ",sqMeanPairwiseDiff," meanLateVar ",meanLateVar];*)


If[Abs[sqMeanPairwiseDiff]<meanLateVar,
(*Print["exiting because differences in mean less than standard deviation"];*)
eqlTime=numsweeps-outWinLength+inWinLength;
Break[]
];(*Exit the while loop and go to the Do loop*)

If[Abs[Tr[sqMeanEMat]]<0.000001,
(*Print["exiting because stuck in a metastable state, sqMeanEMat is ",MatrixForm[sqMeanEMat]];*)
eqlTime=numsweeps-outWinLength+inWinLength;
Break[]
](*Exit the while loop and go to the Do loop*)

];

Entable=RotateRight[#,1]&/@Entable;
(*Shifts every entry to the right by one cyclically. The new data will be written on the first slot, so newer data is at the beginning.*)

];
][[2,1]];


(********************
*  For diagnostics *
********************)

(*Print["AllEnMat ",Transpose[AllEnMat]];*)
(*Print["numsweeps ",numsweeps];*)
(*Print["In IsingEquilibriate, eqlTime ",eqlTime, " Length[AllEnMat] ",Length[AllEnMat]];*)

Print[ListLinePlot[Transpose[AllEnMat[[1;;Min[Length[AllEnMat],2*eqlTime]]]],PlotRange->All,PlotLabel->"t from 1 to 2 times eqlT for beta "<>ToString[beta]]];

(*Print["data",data];*)

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,<|"beta"->beta,"eqlT"->eqlTime,"state"->data[[Key["state"]]]|>};

Remove[{vCount,data,maxGStateCount,sweepOutput,Entable,outWinLength,inWinLength,AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime,maxNumSweeps}];

(*Share[];*)

result

];


(* ::Item::Closed:: *)
(*GraphEquilibriate Overload*)


ECGrav`GraphEquilibriate[seedGraph_List, beta_Real, hamiltonian_,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 10/10/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/22/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or 
* unlabeled graphs
*3. 10/10/2025 Updated created the overload to enable delH to not be specified
*)

(*Equilibriate an input graph configuration to the temperature beta.,
Inputs are:,
1. seedGraph - adjacency matrix of the seed graph,
2. beta = inverse temperature,
3. a hamiltonian, 
4. a formula for delta E (when one edge is flipped to expedite computation),  
5. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Depends on the function GraphSweepReplicas., 

Outputs a list with two associations:, the first one is the lowest energy found  throughout the run together with all states with that lowest energy;, 
The second is the final equilibriated state
*)
Module[{result,vCount=Length[seedGraph],data,maxGStateCount,sweepOutput,Entable,outWinLength,inWinLength,AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime=20000,maxNumSweeps=25000},

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)
data=<|"minEnergy" ->hamiltonian[seedGraph],"minEstates"->{seedGraph},
"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph],"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>,"empty"-><|"graph"->Table[0,{i,vCount},{j,vCount}],"energy" ->0.0,"mag"->0.0|>,
"random"-><|"graph"->With[{rm=RandomInteger[1,{vCount,vCount}]},Mod[rm+Transpose[rm],2]],"energy" ->0.0,"mag"->0.0|>|>;

data[[Key["empty"]]][[Key["energy"]]]=hamiltonian[data[[Key["empty"]]][[Key["graph"]]]];
data[[Key["random"]]][[Key["energy"]]]=hamiltonian[data[[Key["random"]]][[Key["graph"]]]];
data[[Key["random"]]][[Key["mag"]]]=Total[Flatten[data[[Key["random"]]][[Key["graph"]]]]]*1.0/(vCount(vCount-1));


Print["Equilibriating at beta ",beta];
(*Print[" starting with data ",data];*)

numsweeps=0;
outWinLength =500 ;(* length of a table to store running energy values to test equilibriation*)

inWinLength=20;(* A segment to be averaged over at the beginning and the end of the table*)

Entable=Table[0.0,{m,3},{n,outWinLength}];(*A table to store energy values of the seed, random and empty tracks for tunning calculation of tests of equilibriation*)

(*Main calculation*)

AllEnMat=
(*AllEnMat is for storing all energy values for diagnostic*)
Reap[
While[numsweeps<maxNumSweeps,

numsweeps++;

Sow[{data[[Key["state"]]][[Key["energy"]]],data[[Key["empty"]]][[Key["energy"]]],data[[Key["random"]]][[Key["energy"]]]}];

Do[

(*Print["sweeping at ",replicaName];*)
(*Print["minEnergy ",data[[Key["minEnergy"]]]];*)

sweepOutput=ECGrav`GraphSweepReplica[data[[Key[replicaName]]][[Key["graph"]]],beta,hamiltonian,1,data[[Key["minEnergy"]]],UnlabeledVerticesYes];

(*Print["In GraphEquilibriate, sweep output ",sweepOutput];*)

data[[Key[replicaName]]]=sweepOutput[[2]];

(*Update minimum energy states after the sweep*)
If[sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
data[[Key["minEnergy"]]]=sweepOutput[[1,Key["minEnergy"]]];
data[[Key["minEstates"]]]=sweepOutput[[1,Key["minEstates"]]],
If[sweepOutput[[1,Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<=maxGStateCount,
data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1,Key["minEstates"]]]]];
];
,{replicaName,{"state","empty","random"}}];

(*Print["after one sweep, replicas are ",data];*)


Entable[[All,1]]={data[[Key["state"]]][[Key["energy"]]],data[[Key["empty"]]][[Key["energy"]]],data[[Key["random"]]][[Key["energy"]]]};


(*Print["Entable ",Entable//MatrixForm];*)


If[numsweeps>outWinLength,

sqMeanEMat=Table[(Mean[Entable[[i]][[1;;inWinLength]]]-Mean[Entable[[j]][[(outWinLength-inWinLength);;outWinLength]]])^2,{i,3},{j,3}];

(*sqMeanEMat is a three by three matrix of the squared means of difference in energy within and across the tracks at the beginning and end of outWinLength.
At equilibrium, these 9 mumbers should be randomly distributed with mean of 0. Their fluctuation from 0 should be within the variance of the newer (late time) variance in energy for the equilibriation to exit. sqMeanPairwiseDiff is the mean of these 9 numbers*)

sqMeanPairwiseDiff=Mean[Flatten[sqMeanEMat]];

meanLateVar=Mean[Table[Variance[Entable[[i]][[1;;inWinLength]]],{i,3}]];(*Mean variance of the newer energies*)


(*Print["sqMeanEMat ",sqMeanEMat//MatrixForm," sqMeanPairwiseDiff ",sqMeanPairwiseDiff," meanLateVar ",meanLateVar];*)


If[Abs[sqMeanPairwiseDiff]<meanLateVar,
(*Print["exiting because differences in mean less than standard deviation"];*)
eqlTime=numsweeps-outWinLength+inWinLength;
Break[]
];(*Exit the while loop and go to the Do loop*)

If[Abs[Tr[sqMeanEMat]]<0.000001,
(*Print["exiting because stuck in a metastable state, sqMeanEMat is ",MatrixForm[sqMeanEMat]];*)
eqlTime=numsweeps-outWinLength+inWinLength;
Break[]
](*Exit the while loop and go to the Do loop*)

];

Entable=RotateRight[#,1]&/@Entable;
(*Shifts every entry to the right by one cyclically. The new data will be written on the first slot, so newer data is at the beginning.*)

];
][[2,1]];


(********************
*  For diagnostics *
********************)

(*Print["AllEnMat ",Transpose[AllEnMat]];*)
(*Print["numsweeps ",numsweeps];*)
(*Print["In IsingEquilibriate, eqlTime ",eqlTime, " Length[AllEnMat] ",Length[AllEnMat]];*)

Print[ListLinePlot[Transpose[AllEnMat[[1;;Min[Length[AllEnMat],2*eqlTime]]]],PlotRange->All,PlotLabel->"t from 1 to 2 times eqlT for beta "<>ToString[beta]]];

(*Print["data",data];*)

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,<|"beta"->beta,"eqlT"->eqlTime,"state"->data[[Key["state"]]]|>};

Remove[{vCount,data,maxGStateCount,sweepOutput,Entable,outWinLength,inWinLength,AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime,maxNumSweeps}];

(*Share[];*)

result

];


(* ::Item::Closed:: *)
(*GraphEquilibriate Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphEquilibriate[args___]:=(Message[ECGrav`GraphEquilibriate::argerr, args];
$Failed);


(* ::Subsection:: *)
(*GraphComputeCorrelationTime*)


(* ::Item::Closed:: *)
(*GraphComputeCorrelationTime Primary*)


ECGrav`GraphComputeCorrelationTime[seedGraph_List,beta_Real,hamiltonian_,delH_,eqlT_Integer, minEToBeat_Real,EnergyOrMag_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 1/30/2025  ***)
(*************************************)

(*Notes: ,
*1. 01/16/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs*)

(*
Takes an equilibriated graph model with equilibriation time, graph configurations, beta, and minEToBeat and computes the correlation time., 

Inputs are:,
1. seedGraph = List, adjacency matrix of the seed graph,
2. beta = Real, inverse temperature,
3. hamiltonian = formula, the hamiltonian, 
4. delH = a formula for delta E (when one edge is flipped to expedite computation),  
5. eqlT = equilibriation time.,6. minEToBeat = energy value such that if a state (or states) are found with energy lower than this then the lowest such energy and configurations with that energy will be saved.,
7. EnergyOrMag = Integer, a variable to specify whether energy or magnetization is used to 
    compute correlation time. If EnergyOrMag = 0, energy is used; if EnergyOrMag = 1, 
    magnetization is used., 
8. UnlabeledVerticesYes = Integer:  0 means no selection probability to make the graphs unlabeled so graphs are labeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with two associations,
1. the minimum energy visited throughout the equilibriation and the states with that energy. If multiple states have degenerate minimum energy, they will all be included; but if multiple identical adjacency graphs are found, only one unique adjacency graph is kept. ,
2. The second association has the inverse temperature, equilibriation time, 
    correlation timeand the final state visited,   which itself is an association which includes the adjacency matrix, 
   magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.,

Depends on the functions: GraphSweepReplicas, CorrelationTime., 

 *)

Module[{result,vCount=Length[seedGraph],data,maxGStateCount,sweepOutput,EorMTable,corrTable,tmax,norm,numsweeps,EorM = "energy",empair,maxNumSweeps=5000},

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

Which[EnergyOrMag==0,EorM = "energy",EnergyOrMag==1,EorM = "mag"];


data=<|"minEnergy" ->minEToBeat,"minEstates"->{},"beta"->beta,"eqlT"->eqlT,"corrT"->2,"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph],"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>|>;


numsweeps=50*eqlT;

Print["computing correlation time at beta ",beta, " using Energy or Magnetization ",EorM, " numsweeps ",numsweeps];


empair=Reap[EorMTable=
Table[
If[Mod[i,Ceiling[numsweeps/5.0]]==0,Print[" sweepno ",i]];
sweepOutput=ECGrav`GraphSweepReplica[data[[Key["state"]]][[Key["graph"]]],beta,hamiltonian,delH,1,data[[Key["minEnergy"]]],UnlabeledVerticesYes];

(*Print["In compute Corr T, sweep output ",sweepOutput];*)

data[[Key["state"]]]=sweepOutput[[2]];

(*Update minimum energy states after the sweep*)
If[
sweepOutput[[1]][[Key["minEnergy"]]]<data[[Key["minEnergy"]]],
data[[Key["minEnergy"]]]=sweepOutput[[1]][[Key["minEnergy"]]];
data[[Key["minEstates"]]]=sweepOutput[[1]][[Key["minEstates"]]],
If[sweepOutput[[1]][[Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<maxGStateCount,
data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1]][[Key["minEstates"]]]]
];
];

(*Sow[{hamiltonian[data[[Key["state"]]][[Key["graph"]]]],
Total[Flatten[data[[Key["state"]]][[Key["graph"]]]]]*1.0/(vCount(vCount-1))}];*)

data[[Key["state"]]][[Key[EorM]]]
,{i,numsweeps}
];
];

(*For diagnostics*)
(*Print["empair ",empair[[2,1]]];*)

(*Table of energy/magnetizations*)

(*Print[" EorMTable ",EorMTable];*)

If[Length[DeleteDuplicates[EorMTable]]>1,

(*If the magnetizations have all equal value, then correlation time can not be computed, so this loop will be exited with the default corrT left at 2 *)

norm=ECGrav`CorrelationTime[0,EorMTable];


If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)

(*Print["norm ",norm];*)

corrTable=Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,
	Print["      computing corrT at t = ",t]];
	ECGrav`CorrelationTime[t,EorMTable],{t,0,numsweeps-10}
]/norm;


tmax=(FirstPosition[corrTable,_?(#<0&),numsweeps-10])[[1]]; (* A place to stop the integration for calculation of correlation time is when the autocorrelation value first becomes negative*)

(*Print["tmax ",tmax];
Print[" data ",data];
Print["corrTable ",corrTable];*)

Print["Correlation table plot"];
Print[ListLinePlot[corrTable[[1;;Min[Length[corrTable],4*tmax]]],PlotRange->Full,PlotLabel->"t vs auto correlation for beta "<>ToString[beta]]];

data[[Key["corrT"]]]=
Max[Ceiling[Sum[corrTable[[t]],{t,tmax}]],2]

];

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,data[[3;;All]]};

Remove[vCount,data,maxGStateCount,sweepOutput,EorMTable,corrTable,tmax,norm,numsweeps,EorM ,empair,maxNumSweeps];

result

];


(* ::Item::Closed:: *)
(*GraphComputeCorrelationTime Overload*)


ECGrav`GraphComputeCorrelationTime[seedGraph_List,beta_Real,hamiltonian_,eqlT_Integer, minEToBeat_Real,EnergyOrMag_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 1/30/2025  ***)
(*************************************)

(*Notes: ,
*1. 01/16/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or 
	unlabeled graphs
*3. 10/10/2025 Updated created the overload to enable delH to not be specified	
*)

(*
Takes an equilibriated graph model with equilibriation time, graph configurations, beta, and minEToBeat and computes the correlation time., 

Inputs are:,
1. seedGraph = List, adjacency matrix of the seed graph,
2. beta = Real, inverse temperature,
3. hamiltonian = formula, the hamiltonian, 
4. delH = a formula for delta E (when one edge is flipped to expedite computation),  
5. eqlT = equilibriation time.,6. minEToBeat = energy value such that if a state (or states) are found with energy lower than this then the lowest such energy and configurations with that energy will be saved.,
7. EnergyOrMag = Integer, a variable to specify whether energy or magnetization is used to 
    compute correlation time. If EnergyOrMag = 0, energy is used; if EnergyOrMag = 1, 
    magnetization is used., 
8. UnlabeledVerticesYes = Integer:  0 means no selection probability to make the graphs unlabeled so graphs are labeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with two associations,
1. the minimum energy visited throughout the equilibriation and the states with that energy. If multiple states have degenerate minimum energy, they will all be included; but if multiple identical adjacency graphs are found, only one unique adjacency graph is kept. ,
2. The second association has the inverse temperature, equilibriation time, 
    correlation timeand the final state visited,   which itself is an association which includes the adjacency matrix, 
   magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.,

Depends on the functions: GraphSweepReplicas, CorrelationTime., 

 *)

Module[{result,vCount=Length[seedGraph],data,maxGStateCount,sweepOutput,EorMTable,corrTable,tmax,norm,numsweeps,EorM = "energy",empair,maxNumSweeps=5000},

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

Which[EnergyOrMag==0,EorM = "energy",EnergyOrMag==1,EorM = "mag"];


data=<|"minEnergy" ->minEToBeat,"minEstates"->{},"beta"->beta,"eqlT"->eqlT,"corrT"->2,"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph],"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>|>;


numsweeps=50*eqlT;

Print["computing correlation time at beta ",beta, " using Energy or Magnetization ",EorM, " numsweeps ",numsweeps];


empair=Reap[EorMTable=
Table[
If[Mod[i,Ceiling[numsweeps/5.0]]==0,Print[" sweepno ",i]];
sweepOutput=ECGrav`GraphSweepReplica[data[[Key["state"]]][[Key["graph"]]],beta,hamiltonian,1,data[[Key["minEnergy"]]],UnlabeledVerticesYes];

(*Print["In compute Corr T, sweep output ",sweepOutput];*)

data[[Key["state"]]]=sweepOutput[[2]];

(*Update minimum energy states after the sweep*)
If[
sweepOutput[[1]][[Key["minEnergy"]]]<data[[Key["minEnergy"]]],
data[[Key["minEnergy"]]]=sweepOutput[[1]][[Key["minEnergy"]]];
data[[Key["minEstates"]]]=sweepOutput[[1]][[Key["minEstates"]]],
If[sweepOutput[[1]][[Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<maxGStateCount,
data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1]][[Key["minEstates"]]]]
];
];

(*Sow[{hamiltonian[data[[Key["state"]]][[Key["graph"]]]],
Total[Flatten[data[[Key["state"]]][[Key["graph"]]]]]*1.0/(vCount(vCount-1))}];*)

data[[Key["state"]]][[Key[EorM]]]
,{i,numsweeps}
];
];

(*For diagnostics*)
(*Print["empair ",empair[[2,1]]];*)

(*Table of energy/magnetizations*)

(*Print[" EorMTable ",EorMTable];*)

If[Length[DeleteDuplicates[EorMTable]]>1,

(*If the magnetizations have all equal value, then correlation time can not be computed, so this loop will be exited with the default corrT left at 2 *)

norm=ECGrav`CorrelationTime[0,EorMTable];


If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)

(*Print["norm ",norm];*)

corrTable=Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,
	Print["      computing corrT at t = ",t]];
	ECGrav`CorrelationTime[t,EorMTable],{t,0,numsweeps-10}
]/norm;


tmax=(FirstPosition[corrTable,_?(#<0&),numsweeps-10])[[1]]; (* A place to stop the integration for calculation of correlation time is when the autocorrelation value first becomes negative*)

(*Print["tmax ",tmax];
Print[" data ",data];
Print["corrTable ",corrTable];*)

Print["Correlation table plot"];
Print[ListLinePlot[corrTable[[1;;Min[Length[corrTable],4*tmax]]],PlotRange->Full,PlotLabel->"t vs auto correlation for beta "<>ToString[beta]]];

data[[Key["corrT"]]]=
Max[Ceiling[Sum[corrTable[[t]],{t,tmax}]],2]

];

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,data[[3;;All]]};

Remove[vCount,data,maxGStateCount,sweepOutput,EorMTable,corrTable,tmax,norm,numsweeps,EorM ,empair,maxNumSweeps];

result

];


(* ::Item::Closed:: *)
(*GraphComputeCorrelationTime Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphComputeCorrelationTime[args___]:=(Message[ECGrav`GraphComputeCorrelationTime::argerr, args];
$Failed);


(* ::Subsection:: *)
(*GraphMultiHistogram*)


(* ::Item::Closed:: *)
(*GraphMultiHistogram Primary*)


ECGrav`GraphMultiHistogram[seedGraph_List,betaLow_Real,betaHigh_Real,hamiltonian_,
	delH_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 4/01/2025  ***)
(*************************************)

(*Implements the Multiple Histogram Method for the graph models to get a smooth plot of the quantity obs as a function of inverse temperature ranging from betaLow to betaHigh.,

(*Notes: ,
*1. 01/20/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs.,
*3. 3/22/2025 Update - wrote an overload to enable inputing a table of beta values,
*4. 4/01/2025 Update - made computation of -beta*freeenergy parallel, *)

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. betaLow = the lower bound of the inverse temperature,
3. betaHigh = the upper bound of the inverse temperatire,
4. hamiltonian =  function that assigns graphs energy,
5. delH = function that gives the change in energy when a single edge is flipped,
6. obs = the observable quantity in question,
7. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
8. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with three entries:,
1. an association of the minimum energy found from the run and the states found having that energy.,
2. an association with temperatures as keys and values of negative*beta*free energy and values for the energies and observables at each temperature,
3. the replicas at the last step *)

Module[{result,vCount=Length[seedGraph],edgeCount,groundStates,maxGStateCount,replicas,Tempoutput,btTable,defaultRatio,EnergyOrMag=1 (*comuting correlation time using energy (0) or magnetization (1)*),measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,curRootSpecificHeat,curRatios,newBetas,computeBFs,energyAsn,mBetaF},

edgeCount=vCount (vCount-1)/2;

replicas=<||>;
chart=<||>;
groundStates=<|"minEnergy"->hamiltonian[seedGraph],"minEstates"->{seedGraph}|>;

(*Print["groundStates ",groundStates];*)

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

btTable={{betaLow,betaHigh}};

defaultRatio={Exp[Log[betaHigh/betaLow]/(vCount*1.0)],Exp[Log[betaLow/betaHigh]/(vCount*1.0)]};
(*A default value for increment of beta if the specific heat happens to be 0. In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective, hence the factor of vCount.*)

Print[" Starting Multihistogram with default ratio ",defaultRatio];


stopnum=1;

While[stopnum<500,

stopnum++;

(*Print["     In While loop btTable ",btTable];*)

(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,i,hamiltonian,delH,UnlabeledVerticesYes],{i,btTable[[-1]]},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["After equilibriating, Tempoutput ",Tempoutput];*)

(*Prepare replicas*)
replicas=Union[replicas,Tempoutput[[All,2]]];

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
	(Values[
	Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],
		Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];

(*Print["     After Euilibriating and updating, replicas ",replicas," groundStates ",groundStates ];*)


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,delH,locrepl[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes],{i,btTable[[-1]]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*Print["     After computing correlation time, Tempoutput is ",Tempoutput];*)

(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,btTable[[-1]]}];

(*Print[" After correlation time, updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
	(Values[
	Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After updating groundStates ",groundStates ];*)


(*
(****************************************)
(* Take NN measurements at the two just  *
* equilibriated temperatures to compute *
* specific heats *)
(****************************************)
*)

numsweeps=0;
measurements=Reap[
While[numsweeps<NN,

numsweeps++;
If[Mod[numsweeps,Ceiling[NN/5.0]]==0,Print[" sweepno ",numsweeps]];

(*Print["At the start of parallel table ",replicas];*)

(*Print["corrTs, ",Table[replicas[[Key[i],"corrT"]],{i,btTable[[-1]]}]];*)


Tempoutput=Association[With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
	ParallelTable[
	(*repNumSweeps=replicas[[Key[i]]][[Key["corrT"]]];*)
	(* Each replica will be swept corrT times *)

	(*Print["in replica ", i, " repNumSweeps ",repNumSweeps];*)

	<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,delH,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

	{i,btTable[[-1]]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];

(*Print["Before breaking : "];
Print["replicas ",replicas];
Print["gstates ",groundStates];*)

(*Update replicas*)
Do[replicas[[Key[i]]][[Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,btTable[[-1]]}];

(* Extract the minimum energy and corresponding states *)

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
	(Values[
	Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After sweep, replicas",replicas, " groundStates ",groundStates ];*)




If[Mod[numsweeps,1]==0,
Table[

Sow[Flatten[{numsweeps,i,replicas[[Key[i]]][[Key["state"]]][[Key["energy"]]],Through[obs[replicas[[Key[i]]][[Key["state"]]][[Key["graph"]]]]]}],i]
,{i,btTable[[-1]]}];
];

]][[2]];


(*Print[" btTable ",btTable ];
Print[" measurements ",measurements ];*)

chart=Union[chart,AssociationThread[btTable[[-1]],measurements[[1;;2]]]];

(*Print["chart",chart];*)

If[(btTable[[-1,1]]>btTable[[-1,2]]),Break[]];

(* sqrt of the specific heat at current temp*)
curRootSpecificHeat=Table[Sqrt[ECGrav`SpecificHeat[chart[[Key[i]]][[All,3]],1,i]],{i,btTable[[-1]]}];
(*Total specific heat, not specific heat per site!*)

(*Print[" curRootSpecificHeat ",curRootSpecificHeat];*)


curRatios=
Table[If[curRootSpecificHeat[[k]]!=0,(1.0+1.0/curRootSpecificHeat[[k]])^k,defaultRatio[[k]]],{k,{1,-1}}];
(*when k=1 it multiplies, and divides when k=2. Ensures that the first factor is greater than 1, and the second less than one.*)

(*Print[" curRatios ",curRatios];*)

newBetas={btTable[[-1,1]]*curRatios[[1]],btTable[[-1,2]]*curRatios[[2]]};

Print["     curRootSpecificHeat ",curRootSpecificHeat," curRatios ",curRatios];
Print[" newBetas ",newBetas];

AppendTo[btTable,newBetas];

Print["btTable ",btTable];

];

(*Print["btTable ",btTable];*)

btTable=Sort[Flatten[btTable]];

(*Print["btTable ",btTable];*)

(*Print[" chart ",chart];*)

(*Print[" replicas ",replicas];*)

(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each temp, i.e., (-beta*F). It is the
* same as computing the partition function, *   
* but better due to overflow issues   *
*        *)
(**********************************************)
*)


energyAsn=<|Table[i->chart[[Key[i]]][[All,3]],{i,Keys[chart]}]|>;

(*Print["energyAsn ",energyAsn];*)

(*mBetaF=computeBFs[energyAsn];*)

mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[energyAsn];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas};

Remove[vCount,edgeCount,groundStates,maxGStateCount,replicas,Tempoutput,btTable,defaultRatio,measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,curRootSpecificHeat,curRatios,newBetas,computeBFs,energyAsn,mBetaF];

result

];


(* ::Item::Closed:: *)
(*GraphMultiHistogram Overload 1 no delH*)


ECGrav`GraphMultiHistogram[seedGraph_List,betaLow_Real,betaHigh_Real,hamiltonian_,
	obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 10/10/2025  ***)
(*************************************)

(*Implements the Multiple Histogram Method for the graph models to get a smooth plot of the quantity obs as a function of inverse temperature ranging from betaLow to betaHigh.,

(*Notes: ,
*1. 01/20/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs.,
*3. 3/22/2025 Update - wrote an overload to enable inputing a table of beta values,
*4. 4/01/2025 Update - made computation of -beta*freeenergy parallel, 
*5. 10/10/2025 Updated created the overload to enable delH to not be specified
*)

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. betaLow = the lower bound of the inverse temperature,
3. betaHigh = the upper bound of the inverse temperatire,
4. hamiltonian =  function that assigns graphs energy,
5. obs = the observable quantity in question,
6. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
7. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with three entries:,
1. an association of the minimum energy found from the run and the states found having that energy.,
2. an association with temperatures as keys and values of negative*beta*free energy and values for the energies and observables at each temperature,
3. the replicas at the last step *)

Module[{result,vCount=Length[seedGraph],edgeCount,groundStates,maxGStateCount,replicas,Tempoutput,btTable,defaultRatio,EnergyOrMag=1 (*comuting correlation time using energy (0) or magnetization (1)*),measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,curRootSpecificHeat,curRatios,newBetas,computeBFs,energyAsn,mBetaF},

edgeCount=vCount (vCount-1)/2;

replicas=<||>;
chart=<||>;
groundStates=<|"minEnergy"->hamiltonian[seedGraph],"minEstates"->{seedGraph}|>;

(*Print["groundStates ",groundStates];*)

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

btTable={{betaLow,betaHigh}};

defaultRatio={Exp[Log[betaHigh/betaLow]/(vCount*1.0)],Exp[Log[betaLow/betaHigh]/(vCount*1.0)]};
(*A default value for increment of beta if the specific heat happens to be 0. In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective, hence the factor of vCount.*)

Print[" Starting Multihistogram with default ratio ",defaultRatio];


stopnum=1;

While[stopnum<500,

stopnum++;

(*Print["     In While loop btTable ",btTable];*)

(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,i,hamiltonian,UnlabeledVerticesYes],{i,btTable[[-1]]},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["After equilibriating, Tempoutput ",Tempoutput];*)

(*Prepare replicas*)
replicas=Union[replicas,Tempoutput[[All,2]]];

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
	(Values[
	Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],
		Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];

(*Print["     After Euilibriating and updating, replicas ",replicas," groundStates ",groundStates ];*)


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,locrepl[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes],{i,btTable[[-1]]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*Print["     After computing correlation time, Tempoutput is ",Tempoutput];*)

(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,btTable[[-1]]}];

(*Print[" After correlation time, updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
	(Values[
	Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After updating groundStates ",groundStates ];*)


(*
(****************************************)
(* Take NN measurements at the two just  *
* equilibriated temperatures to compute *
* specific heats *)
(****************************************)
*)

numsweeps=0;
measurements=Reap[
While[numsweeps<NN,

numsweeps++;
If[Mod[numsweeps,Ceiling[NN/5.0]]==0,Print[" sweepno ",numsweeps]];

(*Print["At the start of parallel table ",replicas];*)

(*Print["corrTs, ",Table[replicas[[Key[i],"corrT"]],{i,btTable[[-1]]}]];*)


Tempoutput=Association[With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
	ParallelTable[
	(*repNumSweeps=replicas[[Key[i]]][[Key["corrT"]]];*)
	(* Each replica will be swept corrT times *)

	(*Print["in replica ", i, " repNumSweeps ",repNumSweeps];*)

	<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

	{i,btTable[[-1]]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];

(*Print["Before breaking : "];
Print["replicas ",replicas];
Print["gstates ",groundStates];*)

(*Update replicas*)
Do[replicas[[Key[i]]][[Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,btTable[[-1]]}];

(* Extract the minimum energy and corresponding states *)

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
	(Values[
	Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After sweep, replicas",replicas, " groundStates ",groundStates ];*)




If[Mod[numsweeps,1]==0,
Table[

Sow[Flatten[{numsweeps,i,replicas[[Key[i]]][[Key["state"]]][[Key["energy"]]],Through[obs[replicas[[Key[i]]][[Key["state"]]][[Key["graph"]]]]]}],i]
,{i,btTable[[-1]]}];
];

]][[2]];


(*Print[" btTable ",btTable ];
Print[" measurements ",measurements ];*)

chart=Union[chart,AssociationThread[btTable[[-1]],measurements[[1;;2]]]];

(*Print["chart",chart];*)

If[(btTable[[-1,1]]>btTable[[-1,2]]),Break[]];

(* sqrt of the specific heat at current temp*)
curRootSpecificHeat=Table[Sqrt[ECGrav`SpecificHeat[chart[[Key[i]]][[All,3]],1,i]],{i,btTable[[-1]]}];
(*Total specific heat, not specific heat per site!*)

(*Print[" curRootSpecificHeat ",curRootSpecificHeat];*)


curRatios=
Table[If[curRootSpecificHeat[[k]]!=0,(1.0+1.0/curRootSpecificHeat[[k]])^k,defaultRatio[[k]]],{k,{1,-1}}];
(*when k=1 it multiplies, and divides when k=2. Ensures that the first factor is greater than 1, and the second less than one.*)

(*Print[" curRatios ",curRatios];*)

newBetas={btTable[[-1,1]]*curRatios[[1]],btTable[[-1,2]]*curRatios[[2]]};

Print["     curRootSpecificHeat ",curRootSpecificHeat," curRatios ",curRatios];
Print[" newBetas ",newBetas];

AppendTo[btTable,newBetas];

Print["btTable ",btTable];

];

(*Print["btTable ",btTable];*)

btTable=Sort[Flatten[btTable]];

(*Print["btTable ",btTable];*)

(*Print[" chart ",chart];*)

(*Print[" replicas ",replicas];*)

(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each temp, i.e., (-beta*F). It is the
* same as computing the partition function, *   
* but better due to overflow issues   *
*        *)
(**********************************************)
*)


energyAsn=<|Table[i->chart[[Key[i]]][[All,3]],{i,Keys[chart]}]|>;

(*Print["energyAsn ",energyAsn];*)

(*mBetaF=computeBFs[energyAsn];*)

mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[energyAsn];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas};

Remove[vCount,edgeCount,groundStates,maxGStateCount,replicas,Tempoutput,btTable,defaultRatio,measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,curRootSpecificHeat,curRatios,newBetas,computeBFs,energyAsn,mBetaF];

result

];


(* ::Item::Closed:: *)
(*GraphMultiHistogram Overload 2 betaTable*)


ECGrav`GraphMultiHistogram[seedGraph_List,btTable_List,hamiltonian_,delH_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 4/01/2025  ***)
(*************************************)

(*An overload which implements the Multiple Histogram Method for the graph models to get a smooth plot of the quantity obs as a function of inverse temperatures given by btTable.,

(*Notes: ,

*)

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. btTable = table of inverse temperatures to be sampled,
3. hamiltonian =  function that assigns graphs energy,
4. delH = function that gives the change in energy when a single edge is flipped,
5. obs = the observable quantity in question,
6. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
7. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with three entries:,
1. an association of the minimum energy found from the run and the states found having that energy.,
2. an association with temperatures as keys and values of negative*beta*free energy and values for the energies and observables at each temperature,
3. the replicas at the last step *)

Module[{result,groundStates,maxGStateCount,replicas,Tempoutput,EnergyOrMag=1 (*comuting correlation time using energy (0) or magnetization (1)*),measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,computeBFs,energyAsn,mBetaF},


replicas=<||>;
chart=<||>;
groundStates=<|"minEnergy"->hamiltonian[seedGraph],"minEstates"->{seedGraph}|>;

(*Print["groundStates ",groundStates];*)

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)


(*In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective.*)

Print[" Starting multihistogram with btTable ",btTable];




(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,i,hamiltonian,delH,UnlabeledVerticesYes],{i,btTable},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
];

(*Print["After equilibriating, Tempoutput ",Tempoutput];*)

(*Prepare replicas*)
replicas=Union[replicas,Tempoutput[[All,2]]];

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
groundStates[[Key["minEnergy"]]]=candminE;
groundStates[[Key["minEstates"]]]=
Union@@
(Values[
Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];

(*Print["     After Euilibriating and updating, replicas ",replicas," groundStates ",groundStates ];*)


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,delH,locrepl[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes],{i,btTable},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*Print["     After computing correlation time, Tempoutput is ",Tempoutput];*)

(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,btTable}];

(*Print[" After correlation time, updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
		(Values[
		Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After updating groundStates ",groundStates ];*)


(*
(****************************************)
(* Take NN measurements at the           *
* equilibriated temperatures to compute *
* specific heats *)
(****************************************)
*)

numsweeps=0;
measurements=Reap[
While[numsweeps<NN,

numsweeps++;
If[Mod[numsweeps,Ceiling[NN/5.0]]==0,Print[" sweepno ",numsweeps]];

(*Print["At the start of parallel table ",replicas];*)

(*Print["corrTs, ",Table[replicas[[Key[i],"corrT"]],{i,btTable[[-1]]}]];*)


Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[
		(*repNumSweeps=replicas[[Key[i]]][[Key["corrT"]]];*)
		(* Each replica will be swept corrT times *)

		(*Print["in replica ", i, " repNumSweeps ",repNumSweeps];*)

		<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,delH,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

		{i,btTable},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];

(*Print["Before breaking : "];
Print["replicas ",replicas];
Print["gstates ",groundStates];*)

(*Update replicas*)
Do[replicas[[Key[i]]][[Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,btTable}];

(* Extract the minimum energy and corresponding states *)

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
groundStates[[Key["minEnergy"]]]=candminE;
groundStates[[Key["minEstates"]]]=
Union@@
(Values[
Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
	groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After sweep, replicas",replicas, " groundStates ",groundStates ];*)




If[Mod[numsweeps,1]==0,
Table[

Sow[Flatten[{numsweeps,i,replicas[[Key[i]]][[Key["state"]]][[Key["energy"]]],Through[obs[replicas[[Key[i]]][[Key["state"]]][[Key["graph"]]]]]}],i]
,{i,btTable}];
];

]][[2]];


(*Print[" btTable ",btTable ];
Print[" measurements ",measurements ];*)

chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];

(*Print["chart",chart];*)

(*Print["btTable ",btTable];*)

(*Print[" replicas ",replicas];*)

(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each temp, i.e., (-beta*F). It is the
* same as computing the partition function, *   
* but better due to overfloq issues   *
*        *)
(**********************************************)
*)



energyAsn=<|Table[i->chart[[Key[i]]][[All,3]],{i,Keys[chart]}]|>;

(*Print["energyAsn ",energyAsn];*)

mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[energyAsn];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas};

Remove[groundStates,maxGStateCount,replicas,Tempoutput,measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,computeBFs,energyAsn,mBetaF];

result

];


(* ::Item::Closed:: *)
(*GraphMultiHistogram Overload 3 betaTable and no delH*)


ECGrav`GraphMultiHistogram[seedGraph_List,btTable_List,hamiltonian_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 10/10/2025  ***)
(*************************************)

(*An overload which implements the Multiple Histogram Method for the graph models 
to get a smooth plot of the quantity obs as a function of inverse temperatures given 
by btTable.,

(*Notes: ,
*1. 10/10/2025 Updated created the overload to enable delH to not be specified
*)

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. btTable = table of inverse temperatures to be sampled,
3. hamiltonian =  function that assigns graphs energy,
4. obs = the observable quantity in question,
5. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
6. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with three entries:,
1. an association of the minimum energy found from the run and the states found having that energy.,
2. an association with temperatures as keys and values of negative*beta*free energy and values for the energies and observables at each temperature,
3. the replicas at the last step *)

Module[{result,groundStates,maxGStateCount,replicas,Tempoutput,EnergyOrMag=1 (*comuting correlation time using energy (0) or magnetization (1)*),measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,computeBFs,energyAsn,mBetaF},


replicas=<||>;
chart=<||>;
groundStates=<|"minEnergy"->hamiltonian[seedGraph],"minEstates"->{seedGraph}|>;

(*Print["groundStates ",groundStates];*)

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)


(*In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective.*)

Print[" Starting multihistogram with btTable ",btTable];




(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,i,hamiltonian,UnlabeledVerticesYes],{i,btTable},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
];

(*Print["After equilibriating, Tempoutput ",Tempoutput];*)

(*Prepare replicas*)
replicas=Union[replicas,Tempoutput[[All,2]]];

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
groundStates[[Key["minEnergy"]]]=candminE;
groundStates[[Key["minEstates"]]]=
Union@@
(Values[
Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];

(*Print["     After Euilibriating and updating, replicas ",replicas," groundStates ",groundStates ];*)


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,locrepl[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes],{i,btTable},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*Print["     After computing correlation time, Tempoutput is ",Tempoutput];*)

(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,btTable}];

(*Print[" After correlation time, updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
		(Values[
		Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After updating groundStates ",groundStates ];*)


(*
(****************************************)
(* Take NN measurements at the           *
* equilibriated temperatures to compute *
* specific heats *)
(****************************************)
*)

numsweeps=0;
measurements=Reap[
While[numsweeps<NN,

numsweeps++;
If[Mod[numsweeps,Ceiling[NN/5.0]]==0,Print[" sweepno ",numsweeps]];

(*Print["At the start of parallel table ",replicas];*)

(*Print["corrTs, ",Table[replicas[[Key[i],"corrT"]],{i,btTable[[-1]]}]];*)


Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[
		(*repNumSweeps=replicas[[Key[i]]][[Key["corrT"]]];*)
		(* Each replica will be swept corrT times *)

		(*Print["in replica ", i, " repNumSweeps ",repNumSweeps];*)

		<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],i,hamiltonian,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

		{i,btTable},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];

(*Print["Before breaking : "];
Print["replicas ",replicas];
Print["gstates ",groundStates];*)

(*Update replicas*)
Do[replicas[[Key[i]]][[Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,btTable}];

(* Extract the minimum energy and corresponding states *)

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
groundStates[[Key["minEnergy"]]]=candminE;
groundStates[[Key["minEstates"]]]=
Union@@
(Values[
Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]),

If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
	groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@(Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]])
		];
	];
];


(*Print["     After sweep, replicas",replicas, " groundStates ",groundStates ];*)




If[Mod[numsweeps,1]==0,
Table[

Sow[Flatten[{numsweeps,i,replicas[[Key[i]]][[Key["state"]]][[Key["energy"]]],Through[obs[replicas[[Key[i]]][[Key["state"]]][[Key["graph"]]]]]}],i]
,{i,btTable}];
];

]][[2]];


(*Print[" btTable ",btTable ];
Print[" measurements ",measurements ];*)

chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];

(*Print["chart",chart];*)

(*Print["btTable ",btTable];*)

(*Print[" replicas ",replicas];*)

(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each temp, i.e., (-beta*F). It is the
* same as computing the partition function, *   
* but better due to overfloq issues   *
*        *)
(**********************************************)
*)



energyAsn=<|Table[i->chart[[Key[i]]][[All,3]],{i,Keys[chart]}]|>;

(*Print["energyAsn ",energyAsn];*)

mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[energyAsn];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas};

Remove[groundStates,maxGStateCount,replicas,Tempoutput,measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,computeBFs,energyAsn,mBetaF];

result

];


(* ::Item::Closed:: *)
(*GraphMultiHistogram Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphMultiHistogram[args___]:=(Message[ECGrav`GraphMultiHistogram::argerr, args];
$Failed);


(* ::Subsection:: *)
(*GraphCEITempSchedule*)


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Primary*)


ECGrav`GraphCEITempSchedule[seedGraph_List,betaLow_Real,betaHigh_Real,hamiltonian_,delH_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 08/08/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/18/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs.,
*3. 08/08/2025 Update - created an overload that accepts an inverse temperature table to be used in the GraphMultiHistogram instead of just betaLow and betaHigh *)
(*A program to find a temperature schedule for parallel tempering for a graph model with a given hamiltonian and input seed graph.,

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. betaLow = the lower bound of the inverse temperature,
3. betaHigh = the upper bound of the inverse temperature, the overload combines betaLow and betaHigh into a single input that is a list of temperatures.,
4. hamiltonian =  function that assigns graphs energy,
5. delH = function that gives the change in energy when a single edge is flipped,
6. obs = the observable quantity in question,
7. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
8. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with the following elements:,
1. an association of the minimum energy found from the run and the states found having that energy.,
2. temperature entropy pairs,i.e., the temperature schedule,
3. an association of temperatures with minus*beta*free energy and list of energy and obs values sampled during the multihistogram process.,
4. an association of the last state of the replicas (association of temperatures with equilibrium time, correlation time, energy, magnetization, and last graph configuration)

*)
Module[{result,vCount=Length[seedGraph],hist,minusbetaFTable,chart,
       betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,
       entropyVals,tempSchedule},


hist=ECGrav`GraphMultiHistogram[seedGraph,betaLow,betaHigh,hamiltonian,delH,obs,NN,UnlabeledVerticesYes];

(*Print[" In GraphCEITempSchedule, hist[[1]] ",hist[[1]]];*)

minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)


incrementVal=(1.0/betaLow-1.0/betaHigh)/(20*vCount);
Print["incrementVal ",incrementVal];

CvOverTtab=With[{mbf=minusbetaFTable,msrments=chart},
		ParallelTable[{i,ECGrav`CvOverT[1.0/i,mbf,msrments]},
			{i,1.0/betaHigh,1.0/betaLow,incrementVal},
			DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["CvOverT table 1 through 100 ",CvOverTtab];*)

Print["CvOverT plot "];
Print[ListLinePlot[CvOverTtab,PlotRange->All,PlotLabel->"Cv/T plot"]];


(*entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]-1}];*)

entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]}];

delS=entropyRange/vCount;

Print["entropyRange ",entropyRange, " delS ",delS ];

entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];

(*Print["entropyVals ",entropyVals];*)

tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];

Print[" tempSchedule in T (not beta) ",tempSchedule];

result={hist[[1]],tempSchedule,hist[[2]],hist[[3]]};

Remove[vCount,hist,minusbetaFTable,chart,betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule];

result
];


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Overload 1 no delH*)


ECGrav`GraphCEITempSchedule[seedGraph_List,betaLow_Real,betaHigh_Real,hamiltonian_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 10/10/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/18/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or 
	unlabeled graphs.,
*3. 08/08/2025 Update - created an overload that accepts an inverse temperature table 
	to be used in the GraphMultiHistogram instead of just betaLow and betaHigh, 
*4. 10/10/2025 Updated created the overload to enable delH to not be specified
*)

(*A program to find a temperature schedule for parallel tempering for a graph model 
	with a given hamiltonian and input seed graph.,

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. betaLow = the lower bound of the inverse temperature,
3. betaHigh = the upper bound of the inverse temperature, the overload combines betaLow and betaHigh into a single input that is a list of temperatures.,
4. hamiltonian =  function that assigns graphs energy,
5. obs = the observable quantity in question,
6. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
7. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with the following elements:,
1. an association of the minimum energy found from the run and the states found having that energy.,
2. temperature entropy pairs,i.e., the temperature schedule,
3. an association of temperatures with minus*beta*free energy and list of energy and obs values sampled during the multihistogram process.,
4. an association of the last state of the replicas (association of temperatures with equilibrium time, correlation time, energy, magnetization, and last graph configuration)

*)
Module[{result,vCount=Length[seedGraph],hist,minusbetaFTable,chart,
       betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,
       entropyVals,tempSchedule},


hist=ECGrav`GraphMultiHistogram[seedGraph,betaLow,betaHigh,hamiltonian,obs,NN,UnlabeledVerticesYes];

(*Print[" In GraphCEITempSchedule, hist[[1]] ",hist[[1]]];*)

minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)


incrementVal=(1.0/betaLow-1.0/betaHigh)/(20*vCount);
Print["incrementVal ",incrementVal];

CvOverTtab=With[{mbf=minusbetaFTable,msrments=chart},
		ParallelTable[{i,ECGrav`CvOverT[1.0/i,mbf,msrments]},
			{i,1.0/betaHigh,1.0/betaLow,incrementVal},
			DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["CvOverT table 1 through 100 ",CvOverTtab];*)

Print["CvOverT plot "];
Print[ListLinePlot[CvOverTtab,PlotRange->All,PlotLabel->"Cv/T plot"]];


(*entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]-1}];*)

entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]}];

delS=entropyRange/vCount;

Print["entropyRange ",entropyRange, " delS ",delS ];

entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];

(*Print["entropyVals ",entropyVals];*)

tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];

Print[" tempSchedule in T (not beta) ",tempSchedule];

result={hist[[1]],tempSchedule,hist[[2]],hist[[3]]};

Remove[vCount,hist,minusbetaFTable,chart,betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule];

result
];


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Overload 2 betaTable*)


ECGrav`GraphCEITempSchedule[seedGraph_List,btTable_List,hamiltonian_,delH_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
Module[{result,vCount=Length[seedGraph],hist,minusbetaFTable,chart,betas,incrementVal,
        CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule},


hist=ECGrav`GraphMultiHistogram[seedGraph,btTable,hamiltonian,delH,obs,NN,UnlabeledVerticesYes];

(*Print[" In GraphCEITempSchedule, hist[[1]] ",hist[[1]]];*)

minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)



incrementVal=(1.0/Min[btTable]-1.0/Max[btTable])/(20*vCount);
Print["incrementVal ",incrementVal];

CvOverTtab=With[{mbf=minusbetaFTable,msrments=chart},
		ParallelTable[{i,ECGrav`CvOverT[1.0/i,mbf,msrments]},
		{i,1.0/Max[btTable],1.0/Min[btTable],incrementVal},
		DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["CvOverT table 1 through 100 ",CvOverTtab];*)

Print["CvOverT plot "];
Print[ListLinePlot[CvOverTtab,PlotRange->All,PlotLabel->"Cv/T plot"]];


(*entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]-1}];*)

entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]}];

delS=entropyRange/vCount;

Print["entropyRange ",entropyRange, " delS ",delS ];

entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];

(*Print["entropyVals ",entropyVals];*)

tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];

Print[" tempSchedule in T (not beta) as a list of {Temp,Entropy} values ",tempSchedule];

result={hist[[1]],tempSchedule,hist[[2]],hist[[3]]};

Remove[vCount,hist,minusbetaFTable,chart,betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule];

result
];


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Overload 3 betaTable and no delH*)


ECGrav`GraphCEITempSchedule[seedGraph_List,btTable_List,hamiltonian_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
Module[{result,vCount=Length[seedGraph],hist,minusbetaFTable,chart,betas,incrementVal,
        CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule},


hist=ECGrav`GraphMultiHistogram[seedGraph,btTable,hamiltonian,obs,NN,UnlabeledVerticesYes];

(*Print[" In GraphCEITempSchedule, hist[[1]] ",hist[[1]]];*)

minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)



incrementVal=(1.0/Min[btTable]-1.0/Max[btTable])/(20*vCount);
Print["incrementVal ",incrementVal];

CvOverTtab=With[{mbf=minusbetaFTable,msrments=chart},
		ParallelTable[{i,ECGrav`CvOverT[1.0/i,mbf,msrments]},
		{i,1.0/Max[btTable],1.0/Min[btTable],incrementVal},
		DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["CvOverT table 1 through 100 ",CvOverTtab];*)

Print["CvOverT plot "];
Print[ListLinePlot[CvOverTtab,PlotRange->All,PlotLabel->"Cv/T plot"]];


(*entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]-1}];*)

entropyRange=Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,Length[CvOverTtab]}];

delS=entropyRange/vCount;

Print["entropyRange ",entropyRange, " delS ",delS ];

entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];

(*Print["entropyVals ",entropyVals];*)

tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];

Print[" tempSchedule in T (not beta) as a list of {Temp,Entropy} values ",tempSchedule];

result={hist[[1]],tempSchedule,hist[[2]],hist[[3]]};

Remove[vCount,hist,minusbetaFTable,chart,betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule];

result
];


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphCEITempSchedule[args___]:=(Message[ECGrav`GraphCEITempSchedule::argerr, args];
$Failed);


(* ::Subsection:: *)
(*GraphParallelTempering*)


(* ::Item::Closed:: *)
(*GraphParallelTempering Primary*)


ECGrav`GraphParallelTempering[seedGraph_List, btTable_List,minEtoBeat_Real,
hamiltonian_,delH_,obs_,EnergyOrMag_Integer,NN_Integer,numberOfDataPoints_Integer,
UnlabeledVerticesYes_Integer]:=

(*************************************)
(***  Last updated on: 1/30/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs*)
(*
Implements Parallel tempering algorithm on graph models at several different temperatures determined by Constant Entropy Increase (CEI) temperature schedule. It equilibriates, computes correlation time, and then applies measurements. It does temperature swaps during the measurement step. It is parallelized so that equilibriation, computation of correlation time, and sweeps during measurement are all done in parallel. Temperature swaps are done on the master kernel. ,

Depends on the functions: GraphSweepReplicas, CorrelationTime, GraphEquilibriate, GraphComputeCorrelationTime., 

Inputs are:, 
1. seedGraph - adjacency matrix of the input seed graph,
2. btTable - a list of inverse temperatures for the replicas,
3. minEtoBeat -  seed minimum energy at or below which to save graphs,
4. hamiltonian =  function that assigns graphs energy,
5. delH = function that gives the change in energy when a single edge is flipped,
6. obs = the observable quantity in question,
7. EnergyOrMag = an integer (0 or 1) to specify whether to use energy (the default at 0) or magnetization (1) for computing correlation time ,
8. NN = number of independent sweeps to be carried out(so that actual number of sweeps is correlation time times NN),
9. numberOfDataPoints = number of data points of measurements of the observables to be returned.,
10. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with four objects:, 1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. an association of each replica and its beta history ,
4. the final state of each replicas ,
*)

Module[{result,vCount=Length[seedGraph],groundStates,histories,maxGStateCount=500,
replicas,replicaKeysOrderedByBeta,numRep,bt,minStates,candminE,
measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},

groundStates=<||>;
histories=<||>;
numRep=Length[btTable];

Print["Running Parallel Tempering for graph with vCount ",vCount," number of replicas ",
      numRep," betaTable ",btTable];

(*(******************************)
(**      Equilibriate          **)
(******************************)*)
Tempoutput=Association[ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,btTable[[i]],
   hamiltonian,delH,UnlabeledVerticesYes],{i,Length[btTable]},
   DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["In GraphParallelTempering after equilibriating, Tempoutput is ",Tempoutput];*)

(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

(*Print["updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states found from the equilibriation run and the temperature schedule run*)

AppendTo[groundStates,"minEnergy"->Min[Tempoutput[[All,1,"minEnergy"]]]];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)

AppendTo[groundStates,"minEstates"->Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==groundStates[[Key["minEnergy"]]]&][[All,"minEstates"]]]];

(*compare the minimum energy from the equilibriation run with that of the temperature schedule run and reset the minimum energy and states *)

If[minEtoBeat<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=minEtoBeat;
	groundStates[[Key["minEstates"]]]={}
];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)
(*Print["minStates",replicas[[Key["minEstates"]]] ];*)

(*Print["After extracting min energy and states, groundStates ",groundStates];*)

(*(******************************)
(*  Compute Correlation times *)
(******************************)*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],
		     Key["state"],Key["graph"]]],replicas[[Key[i],Key["beta"]]],hamiltonian,
		     delH,replicas[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,
		     UnlabeledVerticesYes],{i,numRep},
		     DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
	]
];


(*Print["In GraphParallelTempering after computing correlation time, Tempoutput is ",Tempoutput];*)

(*Update replicas*)
Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];

(*Print["updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
		Values[
			Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]],
	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]
		];
	];
];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)
(*Print["replicas",replicas ];
Print["groundStates ",groundStates ];*)


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{numberOfDataPoints}]|>,{i,numRep}]|>;



Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=
(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Block[{thisReplInd,nextReplInd,lowTempInd,highTempInd,delBeta,EhighTemp,ElowTemp, 
	delE,expdelBetadelE,accept,blowTemp,bhighTemp,tempHTState},

(*Print[""];
Print[""];
Print[""];
Print[""];
Print["Starting SwapReplicas with offset ",offset];

Print[" replicas ",replicas];
Print[" histories ",histories];*)

(*Print[" replicaKeysOrderedByBeta ",replicaKeysOrderedByBeta];*)

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	(*Print[" thisReplInd ",thisReplInd," nextReplInd ",nextReplInd];*)

	(*Get the low temperature and high temperature betas. 
		Note, the replica indices are sorted in increasing beta, 
		so thisRepInd has the lower beta with the exception being the last and 
		first replica indices due to cycling. Eitherway the next step will work 
		in general *)

	blowTemp=Max[replicas[[Key[thisReplInd],Key["beta"]]],
				replicas[[Key[nextReplInd],Key["beta"]]]]; 
	bhighTemp=Min[replicas[[Key[thisReplInd],Key["beta"]]],
				replicas[[Key[nextReplInd],Key["beta"]]]];
	(*Note blowTemp > bhighTemp so that Thigh > Tlow. *)

	delBeta=blowTemp-bhighTemp;

	(*Print[" blowTemp ",blowTemp," bhighTemp ",bhighTemp," delBeta ",delBeta];*)

	{lowTempInd,highTempInd}=If[replicas[[Key[thisReplInd],Key["beta"]]]>
								replicas[[Key[nextReplInd],Key["beta"]]],
								{thisReplInd,nextReplInd},{nextReplInd,thisReplInd}
								];

	(*Print[" {lowTempInd,highTempInd} = ",{lowTempInd,highTempInd}];*)

	ElowTemp=replicas[[Key[lowTempInd]]][[Key["state"]]][[Key["energy"]]];(*Note blowTemp > bhighTemp and ElowTemp is the energy of the low temperature replica*)

	EhighTemp=replicas[[Key[highTempInd]]][[Key["state"]]][[Key["energy"]]];


	delE = EhighTemp-ElowTemp;


	(*Print["lowTempInd ",lowTempInd, " blowTemp ",blowTemp, " ElowTemp ",ElowTemp];
	Print["highTempInd ",highTempInd, " bhighTemp ",bhighTemp, " EhighTemp ",EhighTemp];*)

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[lowTempInd],Key["swapTry"]]]+=1.0;
	histories[[Key[highTempInd],Key["swapTry"]]]+=1.0;


	(*Print[""];
	Print[""];
	Print["Before swap, replicas ",replicas];
	Print["Before swap, histories ",histories];*)

	accept=0;
	If[delE<=0,accept = 1,
		expdelBetadelE = Exp[-delBeta*delE];
		If[RandomReal[]<expdelBetadelE,accept =1];
	];


	(*Print[" delBeta ",delBeta, " delE ",delE, " accept ",accept];*)

	If[accept==1, (*Do swap of replicas*)

	(*Print["Accepting swap"];*)

	(*Update states: store the state of the high temp replica temporarily*)

	tempHTState=replicas[[Key[highTempInd],Key["state"]]];
	replicas[[Key[highTempInd],Key["state"]]]=replicas[[Key[lowTempInd],Key["state"]]];
	replicas[[Key[lowTempInd],Key["state"]]]=tempHTState;

	(*Print["tempHTState ",tempHTState];
	Print["replicas ",replicas];*)




(*Update histories*)
	If[Mod[numsweeps,Floor[NN/(numberOfDataPoints)]]==0,
		histories[[Key[highTempInd],Key["history"],1]]=
			histories[[Key[lowTempInd],Key["history"],-1]];
		histories[[Key[lowTempInd],Key["history"],1]]=
			histories[[Key[highTempInd],Key["history"],-1]];

		histories[[Key[highTempInd],Key["history"]]]=
			RotateLeft[histories[[Key[highTempInd],Key["history"]]]];
		histories[[Key[lowTempInd],Key["history"]]]=
			RotateLeft[histories[[Key[lowTempInd],Key["history"]]]];
	];

	(*Increase the number of swap accepts by one for both replicas*)
	histories[[Key[highTempInd],Key["swapAccept"]]]+=1.0;
	histories[[Key[lowTempInd],Key["swapAccept"]]]+=1.0 ;

	(*Print[" After swap, histories ",histories],*)

	(*, Print[" Rejecting swap "]; *)

	];


	(*Print["  After one swap, replicas ",replicas];
	Print["  histories ",histories];*)

	,{repl,1,numRep-1,2}];

	(*Print["After swap[], replicas ", replicas, " histories ",histories];*)

];




(*(******************************************
* Start taking measurements and swapping **
******************************************)*)

(*Print[" At the start of taking measurements and swapping, the state of the system is, "];*)
(*Print[" minimum energy so far ",groundStates[["minEnergy"]]];*)
(*Print[" number of states with min energy ",Length[groundStates[["minEstates"]]]];*)
(*Print[" Sample min energy states ",groundStates[["minEstates",1;;Min[10,Length[groundStates[["minEstates"]]]]]]];*)
(*Print[" replicas ",replicas];*)

numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];

(*Print[" numsweeps ",numsweeps," NN ",numsweeps];*)


measurements=Reap[
	While[numsweeps<=NN,

		If[Mod[numsweeps,printCase]==0,Print[" sweepno ",numsweeps]];

		(*Print["At the start of parallel table ",replicas];
		Print["replicas[[All,corrT]]",replicas[[All,"corrT"]]];
		Print[" those with corrT > 4 ",Select[replicas[[All,"corrT"]],#>4&]];*)

		Tempoutput=Association[
			With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
				ParallelTable[
					<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],
						Key["graph"]]],locrepl[[Key[i],Key["beta"]]],hamiltonian,delH,
						locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

					{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}
				]
			]
		];

		(*Print[" after sweep replica Tempoutput ", Tempoutput];*)


		(*Update replicas*)
		Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

		(*Update minimum energy*)
		candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

		(*Print["candminE ", candminE," minEnergy ",groundStates[[Key["minEnergy"]]]];*)

		If[candminE<groundStates[[Key["minEnergy"]]],
			groundStates[[Key["minEnergy"]]]=candminE;
			groundStates[[Key["minEstates"]]]=
			Union@@
				Values[
					Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]
					],
			If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
				groundStates[[Key["minEstates"]]]=
					Union[groundStates[[Key["minEstates"]]],
						Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]
					];
			];
		];

		(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)

		(*Print["replicas ",replicas ];*)
		(*Print["groundStates ",groundStates ];*)


		swap[numsweeps];
		(*Print["After swap, replicas",replicas ];*)

		If[Mod[numsweeps,Floor[NN/numberOfDataPoints]]==0,
			Table[
				bt=replicas[[Key[i],Key["beta"]]];
				Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],
				Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
			,{i,numRep}];
		];

		numsweeps++;

	]
][[2]];

(*Print[" measurements ",measurements ];*)

chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];

(*Print["chart",chart];*)

Do[histories[[Key[i],Key["history"]]]=
	(Select[histories[[Key[i],Key["history"]]],#>=0&]),{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)

result={ groundStates,chart,histories,replicas};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,histories,numRep,bt,minStates,candminE,measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart];

result

];


(* ::Item::Closed:: *)
(*GraphParallelTempering Overload 1 no delH*)


ECGrav`GraphParallelTempering[seedGraph_List, btTable_List,minEtoBeat_Real,
hamiltonian_,obs_,EnergyOrMag_Integer,NN_Integer,numberOfDataPoints_Integer,
UnlabeledVerticesYes_Integer]:=

(*************************************)
(***  Last updated on: 10/10/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or 
	unlabeled graphs,
*3. 10/10/2025 Updated created the overload to enable delH to not be specified,
*)
(*
Implements Parallel tempering algorithm on graph models at several different 
temperatures determined by Constant Entropy Increase (CEI) temperature schedule. 
It equilibriates, computes correlation time, and then applies measurements. 
It does temperature swaps during the measurement step. It is parallelized so that 
equilibriation, computation of correlation time, and sweeps during measurement are 
all done in parallel. Temperature swaps are done on the master kernel. ,

Depends on the functions: GraphSweepReplicas, CorrelationTime, GraphEquilibriate, GraphComputeCorrelationTime., 

Inputs are:, 
1. seedGraph - adjacency matrix of the input seed graph,
2. btTable - a list of inverse temperatures for the replicas,
3. minEtoBeat -  seed minimum energy at or below which to save graphs,
4. hamiltonian =  function that assigns graphs energy,
5. obs = the observable quantity in question,
6. EnergyOrMag = an integer (0 or 1) to specify whether to use energy (the default at 0) or magnetization (1) for computing correlation time ,
7. NN = number of independent sweeps to be carried out(so that actual number of sweeps is correlation time times NN),
8. numberOfDataPoints = number of data points of measurements of the observables to be returned.,
9. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with four objects:, 
1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. an association of each replica and its beta history ,
4. the final state of each replicas ,
*)

Module[{result,vCount=Length[seedGraph],groundStates,histories,maxGStateCount=500,
replicas,replicaKeysOrderedByBeta,numRep,bt,minStates,candminE,
measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},

groundStates=<||>;
histories=<||>;
numRep=Length[btTable];

Print["Running Parallel Tempering for graph with vCount ",vCount," number of replicas ",
      numRep," betaTable ",btTable];

(*(******************************)
(**      Equilibriate          **)
(******************************)*)
Tempoutput=Association[ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,btTable[[i]],
   hamiltonian,UnlabeledVerticesYes],{i,Length[btTable]},
   DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Print["In GraphParallelTempering after equilibriating, Tempoutput is ",Tempoutput];*)

(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

(*Print["updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states found from the equilibriation run and the temperature schedule run*)

AppendTo[groundStates,"minEnergy"->Min[Tempoutput[[All,1,"minEnergy"]]]];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)

AppendTo[groundStates,"minEstates"->Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==groundStates[[Key["minEnergy"]]]&][[All,"minEstates"]]]];

(*compare the minimum energy from the equilibriation run with that of the temperature schedule run and reset the minimum energy and states *)

If[minEtoBeat<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=minEtoBeat;
	groundStates[[Key["minEstates"]]]={}
];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)
(*Print["minStates",replicas[[Key["minEstates"]]] ];*)

(*Print["After extracting min energy and states, groundStates ",groundStates];*)

(*(******************************)
(*  Compute Correlation times *)
(******************************)*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],
		     Key["state"],Key["graph"]]],replicas[[Key[i],Key["beta"]]],hamiltonian,
		     replicas[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,
		     UnlabeledVerticesYes],{i,numRep},
		     DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
	]
];


(*Print["In GraphParallelTempering after computing correlation time, Tempoutput is ",Tempoutput];*)

(*Update replicas*)
Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];

(*Print["updated replicas ",replicas];*)

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE",candminE];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
		Values[
			Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]],
	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]
		];
	];
];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)
(*Print["replicas",replicas ];
Print["groundStates ",groundStates ];*)


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{numberOfDataPoints}]|>,{i,numRep}]|>;



Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=
(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Block[{thisReplInd,nextReplInd,lowTempInd,highTempInd,delBeta,EhighTemp,ElowTemp, 
	delE,expdelBetadelE,accept,blowTemp,bhighTemp,tempHTState},

(*Print[""];
Print[""];
Print[""];
Print[""];
Print["Starting SwapReplicas with offset ",offset];

Print[" replicas ",replicas];
Print[" histories ",histories];*)

(*Print[" replicaKeysOrderedByBeta ",replicaKeysOrderedByBeta];*)

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	(*Print[" thisReplInd ",thisReplInd," nextReplInd ",nextReplInd];*)

	(*Get the low temperature and high temperature betas. 
		Note, the replica indices are sorted in increasing beta, 
		so thisRepInd has the lower beta with the exception being the last and 
		first replica indices due to cycling. Eitherway the next step will work 
		in general *)

	blowTemp=Max[replicas[[Key[thisReplInd],Key["beta"]]],
				replicas[[Key[nextReplInd],Key["beta"]]]]; 
	bhighTemp=Min[replicas[[Key[thisReplInd],Key["beta"]]],
				replicas[[Key[nextReplInd],Key["beta"]]]];
	(*Note blowTemp > bhighTemp so that Thigh > Tlow. *)

	delBeta=blowTemp-bhighTemp;

	(*Print[" blowTemp ",blowTemp," bhighTemp ",bhighTemp," delBeta ",delBeta];*)

	{lowTempInd,highTempInd}=If[replicas[[Key[thisReplInd],Key["beta"]]]>
								replicas[[Key[nextReplInd],Key["beta"]]],
								{thisReplInd,nextReplInd},{nextReplInd,thisReplInd}
								];

	(*Print[" {lowTempInd,highTempInd} = ",{lowTempInd,highTempInd}];*)

	ElowTemp=replicas[[Key[lowTempInd]]][[Key["state"]]][[Key["energy"]]];(*Note blowTemp > bhighTemp and ElowTemp is the energy of the low temperature replica*)

	EhighTemp=replicas[[Key[highTempInd]]][[Key["state"]]][[Key["energy"]]];


	delE = EhighTemp-ElowTemp;


	(*Print["lowTempInd ",lowTempInd, " blowTemp ",blowTemp, " ElowTemp ",ElowTemp];
	Print["highTempInd ",highTempInd, " bhighTemp ",bhighTemp, " EhighTemp ",EhighTemp];*)

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[lowTempInd],Key["swapTry"]]]+=1.0;
	histories[[Key[highTempInd],Key["swapTry"]]]+=1.0;


	(*Print[""];
	Print[""];
	Print["Before swap, replicas ",replicas];
	Print["Before swap, histories ",histories];*)

	accept=0;
	If[delE<=0,accept = 1,
		expdelBetadelE = Exp[-delBeta*delE];
		If[RandomReal[]<expdelBetadelE,accept =1];
	];


	(*Print[" delBeta ",delBeta, " delE ",delE, " accept ",accept];*)

	If[accept==1, (*Do swap of replicas*)

	(*Print["Accepting swap"];*)

	(*Update states: store the state of the high temp replica temporarily*)

	tempHTState=replicas[[Key[highTempInd],Key["state"]]];
	replicas[[Key[highTempInd],Key["state"]]]=replicas[[Key[lowTempInd],Key["state"]]];
	replicas[[Key[lowTempInd],Key["state"]]]=tempHTState;

	(*Print["tempHTState ",tempHTState];
	Print["replicas ",replicas];*)




(*Update histories*)
	If[Mod[numsweeps,Floor[NN/(numberOfDataPoints)]]==0,
		histories[[Key[highTempInd],Key["history"],1]]=
			histories[[Key[lowTempInd],Key["history"],-1]];
		histories[[Key[lowTempInd],Key["history"],1]]=
			histories[[Key[highTempInd],Key["history"],-1]];

		histories[[Key[highTempInd],Key["history"]]]=
			RotateLeft[histories[[Key[highTempInd],Key["history"]]]];
		histories[[Key[lowTempInd],Key["history"]]]=
			RotateLeft[histories[[Key[lowTempInd],Key["history"]]]];
	];

	(*Increase the number of swap accepts by one for both replicas*)
	histories[[Key[highTempInd],Key["swapAccept"]]]+=1.0;
	histories[[Key[lowTempInd],Key["swapAccept"]]]+=1.0 ;

	(*Print[" After swap, histories ",histories],*)

	(*, Print[" Rejecting swap "]; *)

	];


	(*Print["  After one swap, replicas ",replicas];
	Print["  histories ",histories];*)

	,{repl,1,numRep-1,2}];

	(*Print["After swap[], replicas ", replicas, " histories ",histories];*)

];




(*(******************************************
* Start taking measurements and swapping **
******************************************)*)

(*Print[" At the start of taking measurements and swapping, the state of the system is, "];*)
(*Print[" minimum energy so far ",groundStates[["minEnergy"]]];*)
(*Print[" number of states with min energy ",Length[groundStates[["minEstates"]]]];*)
(*Print[" Sample min energy states ",groundStates[["minEstates",1;;Min[10,Length[groundStates[["minEstates"]]]]]]];*)
(*Print[" replicas ",replicas];*)

numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];

(*Print[" numsweeps ",numsweeps," NN ",numsweeps];*)


measurements=Reap[
	While[numsweeps<=NN,

		If[Mod[numsweeps,printCase]==0,Print[" sweepno ",numsweeps]];

		(*Print["At the start of parallel table ",replicas];
		Print["replicas[[All,corrT]]",replicas[[All,"corrT"]]];
		Print[" those with corrT > 4 ",Select[replicas[[All,"corrT"]],#>4&]];*)

		Tempoutput=Association[
			With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
				ParallelTable[
					<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],
						Key["graph"]]],locrepl[[Key[i],Key["beta"]]],hamiltonian,
						locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

					{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}
				]
			]
		];

		(*Print[" after sweep replica Tempoutput ", Tempoutput];*)


		(*Update replicas*)
		Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

		(*Update minimum energy*)
		candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

		(*Print["candminE ", candminE," minEnergy ",groundStates[[Key["minEnergy"]]]];*)

		If[candminE<groundStates[[Key["minEnergy"]]],
			groundStates[[Key["minEnergy"]]]=candminE;
			groundStates[[Key["minEstates"]]]=
			Union@@
				Values[
					Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]
					],
			If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
				groundStates[[Key["minEstates"]]]=
					Union[groundStates[[Key["minEstates"]]],
						Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]
					];
			];
		];

		(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)

		(*Print["replicas ",replicas ];*)
		(*Print["groundStates ",groundStates ];*)


		swap[numsweeps];
		(*Print["After swap, replicas",replicas ];*)

		If[Mod[numsweeps,Floor[NN/numberOfDataPoints]]==0,
			Table[
				bt=replicas[[Key[i],Key["beta"]]];
				Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],
				Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
			,{i,numRep}];
		];

		numsweeps++;

	]
][[2]];

(*Print[" measurements ",measurements ];*)

chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];

(*Print["chart",chart];*)

Do[histories[[Key[i],Key["history"]]]=
	(Select[histories[[Key[i],Key["history"]]],#>=0&]),{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)

result={ groundStates,chart,histories,replicas};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,histories,numRep,bt,minStates,candminE,measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart];

result

];


(* ::Item::Closed:: *)
(*GraphParallelTempering Overload 2 input replica*)


ECGrav`GraphParallelTempering[inputReplicas_Association,minEtoBeat_Real,hamiltonian_,
  delH_,obs_,EnergyOrMag_Integer,NN_Integer,numberOfDataPoints_Integer,
  UnlabeledVerticesYes_Integer]:=

(*************************************)
(***  Last updated on: 1/30/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the sweep is done with labeled or unlabeled graphs*)
(*
This overload takes as input an association that is already equilibriated and with correlation times determined and takes measurements applies Parallel tempering algorithm. It sweeps, swaps, and takes measurements. ,

Depends on the functions: GraphSweepReplicas, CorrelationTime, GraphEquilibriate, GraphComputeCorrelationTime., 

Inputs are:, 
1. inputReplicas - association of replicas with ,
2. minEtoBeat -  seed minimum energy at or below which to save graphs,
3. hamiltonian =  function that assigns graphs energy,
4. delH = function that gives the change in energy when a single edge is flipped,
5. obs = the observable quantity in question,
6. EnergyOrMag = an integer (0 or 1) to specify whether to use energy (the default at 0) for computing correlation time or magnetization (1),
7. NN = number of independent sweeps to be carried out(so that actual number of sweeps is correlation time times NN),
8. numberOfDataPoints = number of data points of measurements of the observables to be returned.,
9. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with four objects:, 1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. an association of each replica and its beta history ,
4. the final state of each replicas ,
*)

Module[{result,numRep=Length[inputReplicas],vCount=Length[inputReplicas[[1,"state","graph"]]],
  btTable=Values[inputReplicas[[All,"beta"]]],
  replicas=inputReplicas,groundStates=<|"minEnergy"->minEtoBeat,"minEstates"->{}|>,
  histories=<||>,maxGStateCount=500,replicaKeysOrderedByBeta,bt,minStates,candminE,
  measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},

Print["Running PT for graph with vCount ",vCount," numRep ",numRep," btTable ",btTable];

(*Print[" At the start of taking measurements and swapping, the state of the system is, "];
Print[" minimum energy so far ",groundStates[["minEnergy"]]];
Print[" number of states with min energy ",Length[groundStates[["minEstates"]]]];
Print[" Sample min energy states ",groundStates[["minEstates",1;;Min[10,Length[groundStates[["minEstates"]]]]]]];
Print[" replicas ",replicas];*)

replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{numberOfDataPoints}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

(*Print["histories initialized to ",histories];*)

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Block[{thisReplInd,nextReplInd,lowTempInd,highTempInd,delBeta,EhighTemp,ElowTemp, delE,expdelBetadelE,accept,blowTemp,bhighTemp,tempHTState},

(*Print[""];
Print[""];
Print[""];
Print[""];
Print["Starting SwapReplicas with offset ",offset];

Print[" replicas ",replicas];
Print[" histories ",histories];*)

(*Print[" replicaKeysOrderedByBeta ",replicaKeysOrderedByBeta];*)

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	(*Print[" thisReplInd ",thisReplInd," nextReplInd ",nextReplInd];*)

	(*Get the low temperature and high temperature betas. Note, the replica indices are sorted in increasing beta, so thisRepInd has the lower beta with the exception being the last and first replica indices due to cycling. Eitherway the next step will work in general *)

	blowTemp=Max[replicas[[Key[thisReplInd],Key["beta"]]],replicas[[Key[nextReplInd],Key["beta"]]]]; 
	bhighTemp=Min[replicas[[Key[thisReplInd],Key["beta"]]],replicas[[Key[nextReplInd],Key["beta"]]]];
	(*Note blowTemp > bhighTemp so that Thigh > Tlow. *)
	
	delBeta=blowTemp-bhighTemp;

	(*Print[" blowTemp ",blowTemp," bhighTemp ",bhighTemp," delBeta ",delBeta];*)

	{lowTempInd,highTempInd}=If[replicas[[Key[thisReplInd],Key["beta"]]]>replicas[[Key[nextReplInd],Key["beta"]]],{thisReplInd,nextReplInd},{nextReplInd,thisReplInd}];

	(*Print[" {lowTempInd,highTempInd} = ",{lowTempInd,highTempInd}];*)

	ElowTemp=replicas[[Key[lowTempInd]]][[Key["state"]]][[Key["energy"]]];(*Note blowTemp > bhighTemp and ElowTemp is the energy of the low temperature replica*)

	EhighTemp=replicas[[Key[highTempInd]]][[Key["state"]]][[Key["energy"]]];


	delE = EhighTemp-ElowTemp;


	(*Print["lowTempInd ",lowTempInd, " blowTemp ",blowTemp, " ElowTemp ",ElowTemp];
	Print["highTempInd ",highTempInd, " bhighTemp ",bhighTemp, " EhighTemp ",EhighTemp];*)

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[lowTempInd],Key["swapTry"]]]+=1.0;
	histories[[Key[highTempInd],Key["swapTry"]]]+=1.0;


	(*Print[""];
	Print[""];
	Print["Before swap, replicas ",replicas];
	Print["Before swap, histories ",histories];*)

	accept=0;
	If[delE<=0,accept = 1,
		expdelBetadelE = Exp[-delBeta*delE];
		If[RandomReal[]<expdelBetadelE,accept =1];
	];


	(*Print[" delBeta ",delBeta, " delE ",delE, " accept ",accept];*)

	If[accept==1, (*Do swap of replicas*)

	(*Print["Accepting swap"];*)

	(*Update states: store the state of the high temp replica temporarily*)

	tempHTState=replicas[[Key[highTempInd],Key["state"]]];
	replicas[[Key[highTempInd],Key["state"]]]=replicas[[Key[lowTempInd],Key["state"]]];
	replicas[[Key[lowTempInd],Key["state"]]]=tempHTState;

	(*Print["tempHTState ",tempHTState];
	Print["replicas ",replicas];*)




	(*Update histories*)
	If[Mod[numsweeps,Floor[NN/(numberOfDataPoints)]]==0,
		histories[[Key[highTempInd],Key["history"],1]]=histories[[Key[lowTempInd],Key["history"],-1]];
		histories[[Key[lowTempInd],Key["history"],1]]=histories[[Key[highTempInd],Key["history"],-1]];

		histories[[Key[highTempInd],Key["history"]]]=RotateLeft[histories[[Key[highTempInd],Key["history"]]]];
		histories[[Key[lowTempInd],Key["history"]]]=RotateLeft[histories[[Key[lowTempInd],Key["history"]]]];
	];

	(*Increase the number of swap accepts by one for both replicas*)
	histories[[Key[highTempInd],Key["swapAccept"]]]+=1.0;
	histories[[Key[lowTempInd],Key["swapAccept"]]]+=1.0 ;

	(*Print[" After swap, histories ",histories],*)

	(*, Print[" Rejecting swap "]; *)

	];


(*Print["  After one swap, replicas ",replicas];
Print["  histories ",histories];*)

,{repl,1,numRep-1,2}];

(*Print["After swap[], replicas ", replicas, " histories ",histories];*)

];




(*(******************************************
* Start taking measurements and swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
While[numsweeps<=NN,

If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

(*Print["At the start of parallel table ",replicas];
Print["replicas[[All,corrT]]",replicas[[All,"corrT"]]];
Print[" those with corrT > 4 ",Select[replicas[[All,"corrT"]],#>4&]];*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[
			<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],
				locrepl[[Key[i],Key["beta"]]],hamiltonian,delH,locrepl[[Key[i],
				Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

		{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
	]
];

(*Print[" after sweep replica Tempoutput ", Tempoutput];*)


(*Update replicas*)
Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

(*Update minimum energy*)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE ", candminE," minEnergy ",groundStates[[Key["minEnergy"]]]];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
		Values[
			Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]],
	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]
		];
	];
];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)

(*Print["replicas ",replicas ];*)
(*Print["groundStates ",groundStates ];*)


swap[numsweeps];
(*Print["After swap, replicas",replicas ];*)

If[Mod[numsweeps,Floor[NN/numberOfDataPoints]]==0,
	Table[
		bt=replicas[[Key[i],Key["beta"]]];
		Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],
			Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
	,{i,numRep}];
];

numsweeps++;

]
][[2]];

(*Print[" measurements ",measurements ];*)

chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];

(*Print["chart",chart];*)

Do[histories[[Key[i],Key["history"]]]=
	(Select[histories[[Key[i],Key["history"]]],#>=0&]),{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total 
length is no more than numberOfDataPoints*)

result={ groundStates,chart,histories,replicas};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,histories,numRep,bt,minStates,candminE,EorM,measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart];

result

];


(* ::Item:: *)
(*GraphParallelTempering Overload 3 input replica no delH*)


ECGrav`GraphParallelTempering[inputReplicas_Association,minEtoBeat_Real,hamiltonian_,
  obs_,EnergyOrMag_Integer,NN_Integer,numberOfDataPoints_Integer,
  UnlabeledVerticesYes_Integer]:=

(*************************************)
(***  Last updated on: 1/30/2025  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the sweep is done with labeled or 
	unlabeled graphs,
*3. 10/10/2025 Updated created the overload to enable delH to not be specified,
*)
(*
This overload takes as input an association that is already equilibriated and with correlation times determined and takes measurements applies Parallel tempering algorithm. It sweeps, swaps, and takes measurements. ,

Depends on the functions: GraphSweepReplicas, CorrelationTime, GraphEquilibriate, GraphComputeCorrelationTime., 

Inputs are:, 
1. inputReplicas - association of replicas with ,
2. minEtoBeat -  seed minimum energy at or below which to save graphs,
3. hamiltonian =  function that assigns graphs energy,
4. obs = the observable quantity in question,
5. EnergyOrMag = an integer (0 or 1) to specify whether to use energy (the default at 0) for computing correlation time or magnetization (1),
6. NN = number of independent sweeps to be carried out(so that actual number of sweeps is correlation time times NN),
7. numberOfDataPoints = number of data points of measurements of the observables to be returned.,
8. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with four objects:, 
1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. an association of each replica and its beta history ,
4. the final state of each replicas ,
*)

Module[{result,numRep=Length[inputReplicas],vCount=Length[inputReplicas[[1,"state","graph"]]],
  btTable=Values[inputReplicas[[All,"beta"]]],
  replicas=inputReplicas,groundStates=<|"minEnergy"->minEtoBeat,"minEstates"->{}|>,
  histories=<||>,maxGStateCount=500,replicaKeysOrderedByBeta,bt,minStates,candminE,
  measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},

Print["Running PT for graph with vCount ",vCount," numRep ",numRep," btTable ",btTable];

(*Print[" At the start of taking measurements and swapping, the state of the system is, "];
Print[" minimum energy so far ",groundStates[["minEnergy"]]];
Print[" number of states with min energy ",Length[groundStates[["minEstates"]]]];
Print[" Sample min energy states ",groundStates[["minEstates",1;;Min[10,Length[groundStates[["minEstates"]]]]]]];
Print[" replicas ",replicas];*)

replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{numberOfDataPoints}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

(*Print["histories initialized to ",histories];*)

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Block[{thisReplInd,nextReplInd,lowTempInd,highTempInd,delBeta,EhighTemp,ElowTemp, delE,expdelBetadelE,accept,blowTemp,bhighTemp,tempHTState},

(*Print[""];
Print[""];
Print[""];
Print[""];
Print["Starting SwapReplicas with offset ",offset];

Print[" replicas ",replicas];
Print[" histories ",histories];*)

(*Print[" replicaKeysOrderedByBeta ",replicaKeysOrderedByBeta];*)

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	(*Print[" thisReplInd ",thisReplInd," nextReplInd ",nextReplInd];*)

	(*Get the low temperature and high temperature betas. Note, the replica indices are sorted in increasing beta, so thisRepInd has the lower beta with the exception being the last and first replica indices due to cycling. Eitherway the next step will work in general *)

	blowTemp=Max[replicas[[Key[thisReplInd],Key["beta"]]],replicas[[Key[nextReplInd],Key["beta"]]]]; 
	bhighTemp=Min[replicas[[Key[thisReplInd],Key["beta"]]],replicas[[Key[nextReplInd],Key["beta"]]]];
	(*Note blowTemp > bhighTemp so that Thigh > Tlow. *)
	
	delBeta=blowTemp-bhighTemp;

	(*Print[" blowTemp ",blowTemp," bhighTemp ",bhighTemp," delBeta ",delBeta];*)

	{lowTempInd,highTempInd}=If[replicas[[Key[thisReplInd],Key["beta"]]]>replicas[[Key[nextReplInd],Key["beta"]]],{thisReplInd,nextReplInd},{nextReplInd,thisReplInd}];

	(*Print[" {lowTempInd,highTempInd} = ",{lowTempInd,highTempInd}];*)

	ElowTemp=replicas[[Key[lowTempInd]]][[Key["state"]]][[Key["energy"]]];(*Note blowTemp > bhighTemp and ElowTemp is the energy of the low temperature replica*)

	EhighTemp=replicas[[Key[highTempInd]]][[Key["state"]]][[Key["energy"]]];


	delE = EhighTemp-ElowTemp;


	(*Print["lowTempInd ",lowTempInd, " blowTemp ",blowTemp, " ElowTemp ",ElowTemp];
	Print["highTempInd ",highTempInd, " bhighTemp ",bhighTemp, " EhighTemp ",EhighTemp];*)

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[lowTempInd],Key["swapTry"]]]+=1.0;
	histories[[Key[highTempInd],Key["swapTry"]]]+=1.0;


	(*Print[""];
	Print[""];
	Print["Before swap, replicas ",replicas];
	Print["Before swap, histories ",histories];*)

	accept=0;
	If[delE<=0,accept = 1,
		expdelBetadelE = Exp[-delBeta*delE];
		If[RandomReal[]<expdelBetadelE,accept =1];
	];


	(*Print[" delBeta ",delBeta, " delE ",delE, " accept ",accept];*)

	If[accept==1, (*Do swap of replicas*)

	(*Print["Accepting swap"];*)

	(*Update states: store the state of the high temp replica temporarily*)

	tempHTState=replicas[[Key[highTempInd],Key["state"]]];
	replicas[[Key[highTempInd],Key["state"]]]=replicas[[Key[lowTempInd],Key["state"]]];
	replicas[[Key[lowTempInd],Key["state"]]]=tempHTState;

	(*Print["tempHTState ",tempHTState];
	Print["replicas ",replicas];*)




	(*Update histories*)
	If[Mod[numsweeps,Floor[NN/(numberOfDataPoints)]]==0,
		histories[[Key[highTempInd],Key["history"],1]]=histories[[Key[lowTempInd],Key["history"],-1]];
		histories[[Key[lowTempInd],Key["history"],1]]=histories[[Key[highTempInd],Key["history"],-1]];

		histories[[Key[highTempInd],Key["history"]]]=RotateLeft[histories[[Key[highTempInd],Key["history"]]]];
		histories[[Key[lowTempInd],Key["history"]]]=RotateLeft[histories[[Key[lowTempInd],Key["history"]]]];
	];

	(*Increase the number of swap accepts by one for both replicas*)
	histories[[Key[highTempInd],Key["swapAccept"]]]+=1.0;
	histories[[Key[lowTempInd],Key["swapAccept"]]]+=1.0 ;

	(*Print[" After swap, histories ",histories],*)

	(*, Print[" Rejecting swap "]; *)

	];


(*Print["  After one swap, replicas ",replicas];
Print["  histories ",histories];*)

,{repl,1,numRep-1,2}];

(*Print["After swap[], replicas ", replicas, " histories ",histories];*)

];




(*(******************************************
* Start taking measurements and swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
While[numsweeps<=NN,

If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

(*Print["At the start of parallel table ",replicas];
Print["replicas[[All,corrT]]",replicas[[All,"corrT"]]];
Print[" those with corrT > 4 ",Select[replicas[[All,"corrT"]],#>4&]];*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[
			<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],
				locrepl[[Key[i],Key["beta"]]],hamiltonian,locrepl[[Key[i],
				Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

		{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
	]
];

(*Print[" after sweep replica Tempoutput ", Tempoutput];*)


(*Update replicas*)
Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

(*Update minimum energy*)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

(*Print["candminE ", candminE," minEnergy ",groundStates[[Key["minEnergy"]]]];*)

If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
		Values[
			Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]],
	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]],Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]
		];
	];
];

(*Print["minEnergy",groundStates[[Key["minEnergy"]]] ];*)

(*Print["replicas ",replicas ];*)
(*Print["groundStates ",groundStates ];*)


swap[numsweeps];
(*Print["After swap, replicas",replicas ];*)

If[Mod[numsweeps,Floor[NN/numberOfDataPoints]]==0,
	Table[
		bt=replicas[[Key[i],Key["beta"]]];
		Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],
			Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
	,{i,numRep}];
];

numsweeps++;

]
][[2]];

(*Print[" measurements ",measurements ];*)

chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];

(*Print["chart",chart];*)

Do[histories[[Key[i],Key["history"]]]=
	(Select[histories[[Key[i],Key["history"]]],#>=0&]),{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total 
length is no more than numberOfDataPoints*)

result={ groundStates,chart,histories,replicas};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,histories,numRep,bt,minStates,candminE,EorM,measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart];

result

];


(* ::Item::Closed:: *)
(*GraphParallelTempering Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphParallelTempering[args___]:=(Message[ECGrav`GraphParallelTempering::argerr, args];
$Failed);


(* ::Title:: *)
(*MCSims Protect and End*)


(* End private context *)
End[]

(* Protect exported symbols *)
Protect @@ Names["ECGrav`MCSims`*"];

EndPackage[]
