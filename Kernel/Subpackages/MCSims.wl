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


(* ::Item::Closed:: *)
(*ComputeMinusBetaTimesFreeEnergy*)


(*Primary pattern*)

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
Module[{output,numTemps=Length[dat],betas=Flatten[Keys[dat]],
ntab=Table[Length[dat[[Key[i]]]],{i,Keys[dat]}],logntab,curFs,nextFs,iterNum,del,deltab},
logntab=Log[ntab*1.0];


curFs=Table[1.0,{k,1,numTemps}];

nextFs=curFs;


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
							logntab[[j]]-curFs[[j]]+(betas[[k]]-betas[[j]])*(dat[[Key[betas[[i]]],s]])
						,{j,1,numTemps}]]
					,{s,1,ntab[[i]]}]
				,{i,1,numTemps}]]
			]
	,{k,1,numTemps},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}];

	del=Sqrt[Sum[(1.0-curFs[[k]]/nextFs[[k]])^2,{k,1,numTemps}]];


	Sow[del];(* For diagnostics *)
	If[Mod[iterNum,500]==0,PrintTemporary["iterNum ",iterNum, " del ",del]];
	curFs=nextFs-(Mean[nextFs]);
	];
];


output=<|Table[Keys[dat][[i]]->curFs[[i]],{i,1,Length[betas]}]|>;

Remove[numTemps,betas,ntab,logntab,curFs,nextFs,iterNum,del,deltab];

output

];


(*Overload for external fields*)
(*Overload pattern*)
ECGrav`ComputeMinusBetaTimesFreeEnergy[dat_Association/;VectorQ[Keys[dat]],bt_Real]:=
(*Computes the value of -beta*free energy at the values of the external field values
given as the keys in the association dat at the inverse temperature bt. 

Inputs is:,
1. dat - an association of external field values (J) and the corresponding 
   list of values of the conjugate field (O) measured at these external field values.
   (i.e., the Hamiltonian is assumed to have the term +J*O)
   e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
2. bt - Real = the inverse temperature

Returns an association with the external field values (J_i) as keys and 
-(betaFixed)(F[betaFixed,J_i] as the values.
*)

Module[{output,numTotSims=Length[dat],extFieldVals=Keys[dat],
ntab=Table[Length[dat[i]],{i,Keys[dat]}],logntab,curFs,nextFs,iterNum,del,deltab},
logntab=Log[ntab*1.0];
(*curFs=Table[LogSumExp[-betas[[k]]*(dat[[Key[betas[[k]]]]])],{k,1,numTemps}];
curFs=curFs-((Max[curFs]+Min[curFs])/2);*)


curFs=Table[1.0,{k,1,numTotSims}];

nextFs=curFs;


iterNum=0;
del=200.0;

(* Collecting delta For diagnostics *)
deltab=Reap[
	While[del>10.0^(-6)&&iterNum<30000,
		iterNum++;

		nextFs=ParallelTable[
			ECGrav`LogSumExp[
				Flatten[
					Table[
						Table[
							-ECGrav`LogSumExp[
								Table[
									logntab[[j]]-curFs[[j]]+
									bt*(extFieldVals[[k]]-extFieldVals[[j]])*(dat[[Key[extFieldVals[[i]]],s]])
								,{j,1,numTotSims}]]
						,{s,1,ntab[[i]]}]
					,{i,1,numTotSims}]]
			]
		,{k,1,numTotSims},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}];


del=Sqrt[Sum[(1.0-curFs[[k]]/nextFs[[k]])^2,{k,1,numTotSims}]];

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


output=<|Table[extFieldVals[[i]]->curFs[[i]],{i,1,Length[extFieldVals]}]|>;

Remove[numTotSims,extFieldVals,
	ntab,logntab,curFs,nextFs,iterNum,del,deltab];

output

];


(*Overload pattern for multiple external fields*)
ECGrav`ComputeMinusBetaTimesFreeEnergy[dat_Association/;MatrixQ[Keys[dat]],bt_Real]:=
(*Computes the value of -beta*free energy at the values of the external field values
given as the keys in the association dat at the inverse temperature bt. 

Inputs is:,
1. dat - an association of external field values (J), assumed to be lists, and the corresponding 
   list of values of the conjugate field (O) measured at these external field values.
   (i.e., the Hamiltonian is assumed to have the term +J*O)
   e.g. <|{0.1} -> {1.1,2.3,5.2}, {2.5} -> {-2.0,-4.3,-20.1}|>
2. bt - Real = the inverse temperature

Returns an association with the external field values (J_i) as keys and 
-(betaFixed)(F[betaFixed,J_i] as the values.
*)

Module[{output,numTotSims=Length[dat],extFieldVals=Keys[dat],
ntab=Table[Length[dat[i]],{i,Keys[dat]}],logntab,curFs,nextFs,iterNum,del,deltab},
logntab=Log[ntab*1.0];
(*curFs=Table[LogSumExp[-betas[[k]]*(dat[[Key[betas[[k]]]]])],{k,1,numTemps}];
curFs=curFs-((Max[curFs]+Min[curFs])/2);*)


curFs=Table[1.0,{k,1,numTotSims}];

nextFs=curFs;


iterNum=0;
del=200.0;

(* Collecting delta For diagnostics *)
deltab=Reap[
	While[del>10.0^(-6)&&iterNum<30000,
		iterNum++;

		nextFs=ParallelTable[
			ECGrav`LogSumExp[
				Flatten[
					Table[
						Table[
							-ECGrav`LogSumExp[
								Table[
									logntab[[j]]-curFs[[j]]+
									bt*(extFieldVals[[k]]-extFieldVals[[j]]) . (dat[[Key[extFieldVals[[i]]],s]])
								,{j,1,numTotSims}]]
						,{s,1,ntab[[i]]}]
					,{i,1,numTotSims}]]
			]
		,{k,1,numTotSims},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}];


del=Sqrt[Sum[(1.0-curFs[[k]]/nextFs[[k]])^2,{k,1,numTotSims}]];


Sow[del];(* For diagnostics *)
	If[Mod[iterNum,500]==0,PrintTemporary["iterNum ",iterNum, " del ",del]];

		curFs=nextFs-(Mean[nextFs]);
	];
];


output=<|Table[extFieldVals[[i]]->curFs[[i]],{i,1,Length[extFieldVals]}]|>;

Remove[numTotSims,extFieldVals,
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

(*************************************)
(***  Last updated on: 03/13/2026  ***)
(*************************************)
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
With[{betas=Keys[minusBetaF],nValuesAssn=Length/@energyMeasurements,
	lognValuesAssn=Log[1.0*Length[#]]&/@energyMeasurements},

	ECGrav`LogSumExp[
		Flatten[
			Table[
				Table[
					-ECGrav`LogSumExp[
						Table[
							lognValuesAssn[[Key[betas[[j]]]]]-minusBetaF[[Key[betas[[j]]]]]
							+(bf-betas[[j]])*(energyMeasurements[[Key[betas[[i]]],s]])
						,{j,1,Length[betas]}
						]
					]
				,{s,1,nValuesAssn[[Key[betas[[i]]]]]}]
			,{i,1,Length[betas]}]
		]
	]
];

(* Overload Pattern *)
ECGrav`NegativeBetaTimesFreeEnergy[betaFixed_Real,targetExtField_List,
	minusBetaF_Association,conjugateExtFieldMeasurements_Association]:=
(*Computes the value of -beta*free energy at the value of the inverse temperature betaFixed and target external field targetExtField., 
Inputs are:,
1.betaFixed - the value of the inverse temperature beta at which the input associations of minusBetaF and energyMeasurements are computed,
2. targetExtField - Real value for the external field at which the extrapolated value of -beta*freeEnergy is wanted,
3. minusBetaF - an association of the external field values and the corresponding value of beta and -beta*freEnergy 
computed at those points  e.g. <|-2.1 -> -20.3, 2.5 -> 34.2|>,
4. conjugateExtFieldMeasurements - an association of the external field values and the corresponding list of energies measured at inverse temperature betaFixed and those external field values. Note, it assumes that the energy values are not multiplied by the external fiel dvalues; i.e., if the total hamiltonian is H = Sum_k J^kE^k, where J^k is the k'th external field, then the energy values in this association are the E^k's, e.g. <|-2.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the external field values which are the keys for both associations have to be equal as sets! Also, the lengths of the lists of energy measurements have to be equal for all external field values. Finally, it is assumed that the input betaFixed is the inverse temperature that both minusBetaF and energyMeasurements are computed.
*)

With[{externalFieldValues=Keys[minusBetaF],
	nValuesAssn=Length/@conjugateExtFieldMeasurements,
	lognValuesAssn=Log[1.0*Length[#]]&/@conjugateExtFieldMeasurements},

ECGrav`LogSumExp[
	Flatten[
		Table[
			Table[
				-ECGrav`LogSumExp[
					Table[
						lognValuesAssn[[Key[j]]]-minusBetaF[[Key[j]]]
						+betaFixed*(targetExtField-j) . (conjugateExtFieldMeasurements[[Key[i],s]])
					,{j,externalFieldValues}]]
			,{s,1,nValuesAssn[[Key[i]]]}]
		,{i,externalFieldValues}]]
	]
];

(* Catch-all Pattern *)
ECGrav`NegativeBetaTimesFreeEnergy[args___]:=(Message[ECGrav`NegativeBetaTimesFreeEnergy::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*ExtrapolatedExpectationValue*)


(* Primary Pattern*)
ECGrav`ExtrapolatedExpectationValue[bf_Real,minusBetaF_Association,
	energyMeasurements_Association, measuredObservableValues_Association]:=
(*************************************)
(***  Last updated on: 03/13/2026  ***)
(*************************************)
(*Computes the extrapolated expectation value of the observable at the target value of 
the inverse temperature bf. It requires four inputs:,
1. bf - Real, the value of the target inverse temperature,
2. minusBetaF - an association of inverse temperatures and the corresponding value of -beta*free energy, e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,
3. energyMeasurements - an association of inverse temperature and the corresponding list of energies measured at that beta., e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta values which are the keys for both associations have to be equal as sets! Also, the lengths of the lists of energy measurements have to be equal for all betas,
4. measuredObservableValues - an association of inverse temperature and the corresponding list of energies measured at that beta., e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta values which are the keys for both associations have to be equal as sets! Also, the lengths of the lists of energy measurements have to be equal for all betas,
*)
With[{betas=Keys[minusBetaF],nValuesAssn=Length/@energyMeasurements},
(Exp[-ECGrav`NegativeBetaTimesFreeEnergy[bf,minusBetaF,energyMeasurements]])*
	Sum[
		Sum[(measuredObservableValues[[Key[betas[[i]]],s]])/
			(Sum[
				nValuesAssn[[Key[betas[[j]]]]]*Exp[(bf-betas[[j]])*energyMeasurements[[Key[betas[[i]]]]][[s]]]*Exp[-minusBetaF[[Key[betas[[j]]]]]]
			,{j,1,Length[betas]}])
		,{s,1,nValuesAssn[[Key[betas[[i]]]]]}]
	,{i,1,Length[betas]}]
];

(* Overload Pattern*)
ECGrav`ExtrapolatedExpectationValue[betaFixed_Real,targetExtField_List,
	minusBetaF_Association,conjugateExtFieldMeasurements_Association,
	measuredObservableValues_Association]:=
(*************************************)
(***  Last updated on: 03/04/2026  ***)
(*************************************)
(*Computes the extrapolated value of the observable O at the fixed value of the 
inverse temperature betaFixed and at the value of the external field targetExtField. 
The microscopic Hamiltonian is assumed to have a term of the form J*E where J is tunable 
external field parameter and E is the conjugate observable. It requires five inputs:,
1. betaFixed - Real, the fixed value of the  inverse temperature,
2. targetExtField - Real, the value of the target external field,
3. minusBetaF - an association of the external field H and the corresponding value of -beta*free energy(H) e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,
4. conjugateExtFieldMeasurements - an association of the external field values J and the corresponding list of energies (list of E's without the J multiplying them)measured at that J. e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>,
5. measuredObservableValues - an association of the external field values J and the corresponding list of meaurements of the observable O measured at those external field values. e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the J values which are the keys for all the associations have to be equal as sets! Also, the lengths of the lists of values have to be equal for all J values,
*)
With[{extFieldVals=Keys[minusBetaF],
	nValuesAssn=Length/@conjugateExtFieldMeasurements,
	minusBetaTimesFreeEnergy=ECGrav`NegativeBetaTimesFreeEnergy[betaFixed,targetExtField,minusBetaF,conjugateExtFieldMeasurements]},

(*Print["betaFixed ",betaFixed," targetExtField ",targetExtField," minusBetaF ",minusBetaF];
Print[" extFieldVals ",extFieldVals," obsLength ",obsLength, " minusBetaTimesFreeEnergy ",minusBetaTimesFreeEnergy];
Print[" measuredObservableValues ",measuredObservableValues];*)

(Exp[-minusBetaTimesFreeEnergy])*
	Sum[
		Sum[(measuredObservableValues[[Key[i],s]])/
			(Sum[
				nValuesAssn[[Key[j]]]*Exp[-minusBetaF[[Key[j]]]]*Exp[betaFixed*(targetExtField-j) . conjugateExtFieldMeasurements[[Key[i],s]]]
			,{j,extFieldVals}])
		,{s,1,nValuesAssn[[Key[i]]]}]
	,{i,extFieldVals}]
];

(* Catch-all Pattern *)
ECGrav`ExtrapolatedExpectationValue[args___]:=(Message[ECGrav`ExtrapolatedExpectationValue::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*Constrained Probability of the Conjugate Field*)


(*Primary definition*)
ECGrav`ConstrainedProbConjugateField[betaFixed_Real,targetExtField_List/;VectorQ[targetExtField]
	,minusBetaF_Association,conjugateExtFieldMeasurements_Association]:=
(*************************************)
(***  Last updated on: 03/13/2026  ***)
(*************************************)

(*Computes the probability density function of the magnetization or conjugate field 
P(beta,M) which depends implicitely on the external field H to which M is the 
conjugate field. 

(*Notes: ,

*)

Inputs are:,
1. seedGrbetaFixed_Real - the inverse temperature,
2. targetExtField_List - the target value of the external field(s) H,
3. targetConjugateExtField_ - a list of variables e.g. M at which the value of the Landau free energy is sought,
4. minusBetaF_Association - an association of external field values as keys and -beta*free Energy at those field values returned from parallel tempering e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,
5. conjugateExtFieldMeasurements_Association - an association with external field as keys and measured values of the conjugate field as values e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the keys of minusBetaF and conjugateExtFieldMeasurements have to be the same. Also, the lengths of the lists of conjugate field measurements have to be equal for all external field values,

The Output is a function of the variable targetConjugateExtField and so can be plotted *)

Module[{externalFieldValues=Keys[minusBetaF],
	nValuesAssn=Length/@conjugateExtFieldMeasurements,
	lognValuesAssn=Log[1.0*Length[#]]&/@conjugateExtFieldMeasurements,
	LogObsLength=Log[Length[conjugateExtFieldMeasurements[[1]]]],
	logWeights,weights,weightedConjugateFiels,jointProbConjFieldDist,marginalProbConjFieldDist, 
	x,vars,minConjFieldVal,maxConjFieldVal,result},

logWeights=Flatten[
	Table[
		Table[
			-betaFixed*targetExtField . conjugateExtFieldMeasurements[[Key[i],s]]
			-ECGrav`LogSumExp[
			Table[
				lognValuesAssn[[Key[j]]]-minusBetaF[[Key[j]]]-betaFixed*(j) . (conjugateExtFieldMeasurements[[Key[i],s]])
			,{j,externalFieldValues}]]
		,{s,1,nValuesAssn[[Key[i]]]}]
	,{i,externalFieldValues}]
];

weights=Exp[logWeights-Max[logWeights]];


weightedConjugateFiels=WeightedData[Join@@Values[conjugateExtFieldMeasurements],weights];


vars=Array[x,Length[targetExtField]];

minConjFieldVal=Table[(0.8)*Min[Min/@conjugateExtFieldMeasurements[[All,All,k]]],{k,1,Length[targetExtField]}];

maxConjFieldVal=Table[(1.2)*Max[Max/@conjugateExtFieldMeasurements[[All,All,k]]],{k,1,Length[targetExtField]}];

Do[If[minConjFieldVal[[k]]<=0.0,minConjFieldVal[[k]]=-maxConjFieldVal[[k]]*(0.2)],{k,1,Length[targetExtField]}];
Do[If[maxConjFieldVal[[k]]<=0.0,maxConjFieldVal[[k]]=-minConjFieldVal[[k]]*(0.2)],{k,1,Length[targetExtField]}];


Which[
	Length[targetExtField]==1,

		marginalProbConjFieldDist=SmoothKernelDistribution[WeightedData[Flatten[Values[conjugateExtFieldMeasurements]],weights]];

		
		result={marginalProbConjFieldDist,minConjFieldVal,maxConjFieldVal},

	Length[targetExtField]>=2,

		jointProbConjFieldDist=SmoothKernelDistribution[weightedConjugateFiels];

		
		marginalProbConjFieldDist=Table[SmoothKernelDistribution[WeightedData[Join@@Values[conjugateExtFieldMeasurements[[All,All,k]]],weights]],{k,1,Length[targetExtField]}];

		
		result={Join[marginalProbConjFieldDist,{jointProbConjFieldDist}],minConjFieldVal,maxConjFieldVal};

	];
	
result

];

(* Catch-all Pattern *)
ECGrav`ConstrainedProbConjugateField[args___]:=(Message[ECGrav`ConstrainedProbConjugateField::argerr, args];
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
With[{Amsq=Am . Am},
(J/2)*(Total[Amsq,2]-Tr[Amsq])+(L/2)*Total[Am,2]
];

(* Catch-all Pattern *)
ECGrav`HIsing[args___[params___]]:=(Message[ECGrav`HIsing::argerr, args];
$Failed);


(* Primary Pattern *)
ECGrav`delHIsing[Am_List,J_,L_,a_Integer,b_Integer]:=
(* computes HIsing[Amnew] - HIsing[Am] where Amnew if found by toggling Am at 
row a and col b. HIsing = J/2 sum_{i!=j}A^2_{ij} + L/2 sum_{i!=j}A_{ij} *)
With[{n=Length[Am],togab=If[Am[[a,b]]==0,1,-1]},
togab*(J*(Total[Am[[a]]]+Total[Am[[b]]]-2Am[[a,b]])+L)
];

(* Catch-all Pattern *)
ECGrav`delHIsing[args___[params___]]:=(Message[ECGrav`delHIsing::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*HWeightedFaceCounts*)


(* Primary Pattern *)
ECGrav`HWeightedFaceCounts[Am_List,J1_,J2_,J3_,J4_, J5_]:=
(*J1*(# of vertices) + J2*(number of edges)+... + J5*(number of 5-cliques).*)
Block[{n=Length[Am],clqs=FindClique[AdjacencyGraph[Am],\[Infinity],All],maxClq,fvect},
maxClq=Max[Length/@clqs];

fvect={n,Total[Am,2]/2,Tr[MatrixPower[Am,3]]/6,0,0};


Do[fvect[[q]]=Length[Union@@(Subsets[#,{q}]&/@clqs)],{q,4,Min[5,maxClq]}];


{J1,J2,J3,J4,J5} . fvect
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
J*(Transpose[degv] . Lmat . degv)[[1,1]]];

(* Catch-all Pattern *)
ECGrav`HLaplacian[args___]:=(Message[ECGrav`HLaplacian::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*H1dCombManifold*)


(* Primary Pattern *)
ECGrav`H1dCombManifold[Amat_List,J_]:=
(*(1/N0)(-2N1+2N3+HIsing+1+)*)
With[{numV=Length[Amat],Amsq=Amat . Amat},
If[numV<=0,Return[0.0]];
(J/(numV))*(
	ECGrav`HWeightedFaceCounts[Amat,0.0,-2.0,2.0,0.0,0.0]
	+(Total[Amsq,2] - Tr[Amsq])/2.0
	+ 1.0*numV
	+ 1.0*Length[ConnectedComponents[AdjacencyGraph[Amat]]]-1.0)

];

(* Catch-all Pattern *)
ECGrav`H1dCombManifold[args___]:=(Message[ECGrav`H1dCombManifold::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*H2dCombManifold*)


(* Primary Pattern *)
ECGrav`H2dCombManifold[Amat_List,J_]:=
(*Gives 2D combinatorial manifolds as ground states*)
With[{numV=Length[Amat],vList=Range[Length[Amat]],edgeList=Select[Position[Amat,1],#[[1]]<#[[2]]&]},
		
J*(
	1.0*ECGrav`HWeightedFaceCounts[Amat,0.0,1.0,-(3.0),(4.0),0.0]

	+(1.0)*Total[Binomial[ECGrav`HyperDeg[Amat,#],2]&/@edgeList]
	
	+(1.0*numV)*Abs[Length[ConnectedComponents[AdjacencyGraph[Amat]]]-1]

	+(1.0)*Total[
		(Abs[Length[ConnectedComponents[AdjacencyGraph[#]]]-1])&/@(ECGrav`Sph[Amat,#]&/@vList)]	
	)
];

(* Catch-all Pattern *)
ECGrav`H2dCombManifold[args___]:=(Message[ECGrav`H2dCombManifold::argerr, args];
$Failed);


(* Primary Pattern *)
ECGrav`delH2dCombManifold[Amat_List,J_,i_Integer,j_Integer]:=
(*the change in H2dCombManifold when edge {i,j} is toggled*)
Module[{numV=Length[Amat],g=AdjacencyGraph[Amat],AmatNew,gNew,sphIntG,sphIntVertices},
	AmatNew=Amat;
	AmatNew[[i,j]]=AmatNew[[j,i]]=Mod[Amat[[i,j]]+1,2];
	gNew=AdjacencyGraph[AmatNew];
	sphIntG=Subgraph[g,Intersection[VertexList[ECGrav`Sph[g,i]],VertexList[ECGrav`Sph[g,j]]]];
	sphIntVertices=VertexList[sphIntG];

	J*(
		Join[{0,-(2*Amat[[i,j]]-1)},PadRight[-(2*Amat[[i,j]]-1)*ECGrav`FVector[sphIntG],3]] . {0.0,1.0,-(3.0),(4.0),0.0}
			+(1.0)*(-(2*Amat[[i,j]]-1)(Sum[ECGrav`HyperDeg[g,k],{k,Sort/@Tuples[{{i,j},sphIntVertices}]}]+Binomial[Length[sphIntVertices],2])+2*Amat[[i,j]]*Length[sphIntVertices])

		+(1.0*numV)*(Length[ConnectedComponents[AdjacencyGraph[AmatNew]]]-Length[ConnectedComponents[AdjacencyGraph[Amat]]])

		+(1.0)*Which[Length[sphIntVertices]>1&&ConnectedGraphQ[sphIntG],(*Print["case 1"];*)
					0,

				(sphIntVertices=={}&&(Total[Amat[[i]]]+Total[Amat[[j]]]>0)&&(Total[Amat[[i]]]==0||Total[Amat[[j]]]==0)),(*Print["case 2"];*)
					0,
				
				(sphIntVertices=={}&&(Total[Amat[[i]]]+Total[Amat[[j]]]==0)),(*Print["case 3"];*)
					2*(2*Amat[[i,j]]-1),

				True,(*Print["case 4"];*)
					Total[(Abs[Length[ConnectedComponents[ECGrav`Sph[gNew,#]]]-1]-Abs[Length[ConnectedComponents[ECGrav`Sph[g,#]]]-1])&/@Join[{i,j},sphIntVertices]]

			]

		)
];

(* Catch-all Pattern *)
ECGrav`delH2dCombManifold[args___]:=(Message[ECGrav`delH2dCombManifold::argerr, args];
$Failed);


(* ::Item::Closed:: *)
(*H2dPseudoCombManifold*)


(* Primary Pattern *)
ECGrav`H2dPseudoCombManifold[Amat_List,J_]:=
(*(1/N0)(-2N1+2N3+HIsing+1)*)
With[{numV=Length[Amat],vList=Range[Length[Amat]],edgeList=Select[Position[Amat,1],#[[1]]<#[[2]]&]},
J*(ECGrav`HWeightedFaceCounts[Amat,0.0,1.0,-(1.0/numV+3.0),4.0,0.0]

	+1.0*Total[Binomial[ECGrav`HyperDeg[Amat,#],2]&/@edgeList]
	+1.0*Length[ECGrav`KpathConnectedComponents[AdjacencyGraph[Amat],2]]-1.0
  )
];

(* Catch-all Pattern *)
ECGrav`H2dPseudoCombManifold[args___]:=(Message[ECGrav`H2dPseudoCombManifold::argerr, args];
$Failed);


(* ::Chapter:: *)
(*Ground State Searches*)


(* ::Section::Closed:: *)
(*Exact Searches*)


(* ::Item::Closed:: *)
(*LowEnergyStates*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`LowEnergyStates[ensemble_List,H_,level_Integer]:=
(****************************************)
(*   (* Last updated 01/15/2025. *) *)
(*   (* Note: list of states as values. *) *) 
(****************************************)
(* Calculates the 
energies of the states in the list using the hamiltonian function H. Outputs an 
ordered association with the lowest klevel energies as the keys and the states 
with those energies as the values.*)

Module[{mastertuplesParts,partialresult,minstates,numkernels},
numkernels=Max[$KernelCount,1];(* avoid division by zero when no parallel kernels are launched *)
mastertuplesParts=With[{partlength=Ceiling[Length[ensemble]/numkernels]},
Partition[ensemble,UpTo[partlength]]];
(*Print[" mastertuplesParts length ",Length/@mastertuplesParts," mastertuplesParts ",
mastertuplesParts];*)


partialresult=Join@@ParallelTable[
	TakeSmallestBy[p,H[#]&,UpTo[level]],{p,mastertuplesParts},
		DistributedContexts->{$Context, "ECGrav`PureComplexes`Private`","ECGrav`MCSims`Private`"}];
	

minstates=TakeSmallestBy[partialresult,H[#]&,UpTo[level]];

(*<|Table[GraphPlot[AdjacencyGraph[i],ImageSize->{70,70}]->H[i],{i,minstates}]|>*)

ECGrav`ChooseNonIsomorphicGraphs[#]&/@(GroupBy[minstates,H[#]&])

];

(* Catch-all Pattern *)
ECGrav`LowEnergyStates[args___]:=(Message[ECGrav`LowEnergyStates::argerr, args];
$Failed);


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
,{i,edgeList}];


negatives=Select[deltaEtable,#<=0&];

If[Length[negatives]>0,
flipSpin=Keys[TakeSmallest[negatives,1]][[1]];


Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2],
Return[0](* unable to find a flip that lowers the energy*)
];

1 (*success*)
];

printCase=Floor[(cutoff*1.0)/5.0];

numsteps=0;
stepagain=1;

While[stepagain==1&&numsteps<cutoff,
If[Mod[numsteps,printCase]==0,PrintTemporary["step number ",numsteps]];

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
,{i,edgeList}];

computeWeights[beta];
flipSpin=RandomChoice[Values[weightsTable]->Keys[weightsTable]];


Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2];


curE=curE+deltaEtable[flipSpin];
If[curE<minE,
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
If[Mod[numsteps,printCase]==0,PrintTemporary["step number ",numsteps]];


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
,{i,edgeList}];

computeWeights[beta];
flipSpin=RandomChoice[Values[weightsTable]->Keys[weightsTable]];


Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2];


curE=curE+deltaEtable[flipSpin];
If[curE<minE,
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
If[Mod[numsteps,printCase]==0,PrintTemporary["step number ",numsteps]];


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
ECGrav`SimulatedAnnealing[seedAmat_List, hamiltonian_[hparams___],delH_[delHparams___],betai_Real,betaf_Real,roundLength_Integer,NN_Integer]:=                     
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

curE = SetPrecision[hamiltonian[Amcur,hparams],precision];
minE = curE;
excitedEnergy=N[minE+SetPrecision[1000.0,precision],precision];
minStates={Amcur};
excitedStates={Amcur};

beta = betai;
rate = SetPrecision[Exp[Log[betai/betaf]*roundLength/NN],precision];

 
(*Print["state ", state];
Print["finalstate ", finalState];*)

step[bta_]:=
Block[{delE,expdelE,flipSpin,accept},

	flipSpin=RandomChoice[edgeList];
	delE = N[delH[Amcur,delHparams,flipSpin[[1]],flipSpin[[2]]],precision];
		
	accept = 0;
	If[delE<=0,accept = 1,
		expdelE = N[Exp[-SetPrecision[bta,precision]*delE],precision];
		If[RandomReal[]<=expdelE,accept =1];
	];

If[accept ==0, Return[0]];

If[accept==1,


Amcur[[flipSpin[[1]],flipSpin[[2]]]]=Amcur[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amcur[[flipSpin[[1]],flipSpin[[2]]]]+1,2];

(*
(*For diagnostic*)
	If[Abs[curE+delE-hamiltonian[Amcur,hparams]]>0.000001,
		Print["updated E ",curE+delE, " h[Amnext] ",hamiltonian[Amcur,hparams]," Amnext ", Amcur, " flipped spin ",flipSpin];
	];
(*End diagnostic*)
*)


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

If[Mod[numsteps,printCase]==0,PrintTemporary["sweep number ",numsteps]];

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
ECGrav`SimulatedAnnealing[seedAmat_List, hamiltonian_[hparams___],betai_Real,betaf_Real,roundLength_Integer,NN_Integer]:=                     

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


Module[{ nn=Length[seedAmat],result,maxNumOfSavedStates=400,
	(*precision=100,*)precision=MachinePrecision,Amcur,Amnext,edgeList,minE,
		excitedEnergy,curE,minStates,excitedStates,beta,rate,step,numsteps,printCase},


Amcur=seedAmat;
Amnext=seedAmat;
edgeList=Subsets[Range[nn],{2}];

curE = SetPrecision[hamiltonian[Amcur,hparams],precision];
minE = curE;
excitedEnergy=N[minE+SetPrecision[1000.0,precision],precision];
minStates={Amcur};
excitedStates={Amcur};

beta = betai;
rate = SetPrecision[Exp[Log[betai/betaf]*roundLength/NN],precision];

 
(*Print["state ", state];
Print["finalstate ", finalState];*)

step[bta_Real]:=
Block[{delE,expdelE,flipSpin,accept},

flipSpin=RandomChoice[edgeList];
Amnext[[flipSpin[[1]],flipSpin[[2]]]]=Amnext[[flipSpin[[2]],flipSpin[[1]]]]=Mod[Amnext[[flipSpin[[1]],flipSpin[[2]]]]+1,2];


delE = N[hamiltonian[Amnext,hparams],precision]-curE;
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

If[Mod[numsteps,printCase]==0,PrintTemporary["sweep number ",numsteps]];

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


(* ::Chapter:: *)
(*Exact Expectation Value Calculations*)


(* ::Section::Closed:: *)
(*Expectation Value*)


(* ::Item::Closed:: *)
(*ExactExpectationValue*)


(* :Code Section: *)

(* Primary Pattern *)
ECGrav`ExactExpectationValue[ensemble_List,degenValues_List/;NumericQ[degenValues[[1,1]]],function_,Hamiltonian_,beta_]:=
(*(*****************************)
(* Last Updated: 01/31/2025  *)
(*****************************)*)
(*Notes: ,
1. 01/30/2025 update added the ability to include a degeneracy function that adds 
	additional weight to each graph.,
2. 01/31/2025 update changed ParallelTable to parallel map which has a 
	significantly better memory management for large ensembles*)
(*This program computes the expectation value of the observable function(s). 
It does the same thing as the program ExpectationValue but in parallel. Energies in 
Boltzman weights are computed using the Hamiltonian. degeneracy is vector of functions 
that give the different weights for a given graph. Beta is inverse temperature. 
Outputs a list of lists with as many items as the length of the degeneracy list. 
Each item has in itself the partition function, the free energy, and the expectation 
values of the functions inputed all of which are computed with the corresponding 
degeneracy function. *)

Block[{result,n,nprint,boltzmanWeights,partitionZ,freeF, funcvalues},
n=Length[ensemble];
nprint=Floor[n/5];
result = 0.0;

boltzmanWeights=ParallelMap[
Exp[-beta*Hamiltonian[#]]&,
ensemble,DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`","ECGrav`MCSims`Private`"}];


funcvalues=Transpose[
	ParallelMap[
		Through[function[#]]&,ensemble,
		DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`","ECGrav`MCSims`Private`"}]
	];


partitionZ=Total/@Table[degenValues[[i]]*boltzmanWeights,{i,1,Length[degenValues]}];

freeF=(-1.0/beta)*Log[#]&/@partitionZ;

result=Table[(1.0/partitionZ[[i]])*Table[Total[(degenValues[[i]]*boltzmanWeights*funcvalues[[j]])],{j,1,Length[funcvalues]}],{i,1,Length[degenValues]}];


(*Print["boltzmanWeights ", boltzmanWeights];
Print["funcvalues ", funcvalues];
Print["partitionZ ", partitionZ];
Print["freeF ", freeF];
Print["expvalue ", result];*)

Table[{partitionZ[[i]],freeF[[i]],result[[i]]},{i,1,Length[degenValues]}]
];

ECGrav`ExactExpectationValue[ensemble_List,degenValues_List/;NumericQ[degenValues[[1,1]]],obsValues_List/;NumericQ[obsValues[[1,1]]],Hamiltonian_,beta_]:=
(*(*****************************)
(* Last Updated: 02/07/2025  *)
(*****************************)*)
(*Notes: ,
1. 01/30/2025 update added the ability to include a degeneracy function that adds 
	additional weight to each graph.,
2. 01/31/2025 update changed ParallelTable to parallel map which has a significantly 
	better memory management for large ensembles*)
(*This program computes the expectation value of the observable function(s). 
 Energies in Boltzman weights are computed using the Hamiltonian. ,
Inputs are:,
1. ensemble,
2. degenValues = a list of list of graph multiplicities or degeneracies corresponding 
	to different cases (e.g. labeled vs unlabeled). The list should have as many lists 
	as the number of cases and each element in the list should have as many elements
	as the length of the ensemble,
3. obsValues = a list of list of observable values. Should have the same length as the 
	number of observables and each element of the list should be as long as the 
	ensemble with nth element of the ith list corresponding to the value of the ith 
	observable on the nth graph in the ensemble,
4. Hamiltonian = hamiltonian that maps graphs to reals,
5. beta = inverse temperature 
Outputs a list of lists with as many items as the length of the degeneracy list. Each 
item has in itself the partition function, the free energy, and the expectation values 
of the functions inputed all of which are computed with the corresponding degeneracy 
function. *)

Block[{result,n,nprint,boltzmanWeights,partitionZ,freeF},
n=Length[ensemble];
nprint=Floor[n/5];
result = 0.0;

boltzmanWeights=ParallelMap[
Exp[-beta*Hamiltonian[#]]&,
ensemble,DistributedContexts->{$Context,"ECGrav`PureComplexes`Private`","ECGrav`MCSims`Private`"}];


partitionZ=Total/@Table[degenValues[[i]]*boltzmanWeights,{i,1,Length[degenValues]}];

freeF=(-1.0/beta)*Log[#]&/@partitionZ;

result=Table[(1.0/partitionZ[[i]])*Table[Total[(degenValues[[i]]*boltzmanWeights*obsValues[[j]])],{j,1,Length[obsValues]}],{i,1,Length[degenValues]}];


(*Print["boltzmanWeights ", boltzmanWeights];
Print["obsValues ", obsValues];
Print["partitionZ ", partitionZ];
Print["freeF ", freeF];
Print["expvalue ", result];*)

Table[{partitionZ[[i]],freeF[[i]],result[[i]]},{i,1,Length[degenValues]}]
];

(* Catch-all Pattern *)
ECGrav`ExactExpectationValue[args___]:=(Message[ECGrav`ExactExpectationValue::argerr, args];
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
Module[{vCount=Length[seedAmat],maxEdgeCount,expDEVals,MCStep,MCSweep,numsweeps,Amat,AmatNew,printVal,result},
maxEdgeCount=vCount (vCount-1)/2;
printVal=Floor[maxSweep/5];

Amat = seedAmat;
AmatNew = Amat;


MCStep[]:=Module[{l,lp,m,i,j,deltaE,expBetaDE,accept},(*one MC flip attempt*)
Do[
AmatNew = Amat;
l=RandomInteger[{1,maxEdgeCount}];

lp=maxEdgeCount+1-l;
m=0;
While[lp>0,
m++;
lp-=vCount-m;
];
j=vCount+1-m;
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
If[Mod[numsweeps,printVal]==0,PrintTemporary["Now at sweep number ",numsweeps]];


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
(*GraphSweepReplica Primary*)


ECGrav`GraphSweepReplica[seedGraph_List,beta_Real,hamiltonian_[hparams___],
	delH_[delHparams___],NN_Integer,minEToBeat_Real,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 11/16/2025  ***)
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

data=<|"graph"->seedGraph,"energy"->hamiltonian[seedGraph,hparams], "mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>;

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

delE=delH[data[[Key["graph"]]],delHparams,row,col];


newAmat=data[[Key["graph"]]];
newAmat[[row,col]]=newAmat[[col,row]]=Mod[newAmat[[row,col]]+1,2];

Which[UnlabeledVerticesYes==0,selectionProb=1,UnlabeledVerticesYes==1,
selectionProb=(*Automorphism group order of new graph over auromorphism group order of the old*)
GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[newAmat]]]/GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[data[[Key["graph"]]]]]],True,Abort[];
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


data[[Key["graph"]]]=newAmat;
data[[Key["energy"]]]+= delE;
data[[Key["mag"]]]+= (2.0/(vCount(vCount-1)))*(2*data[[Key["graph"]]][[row,col]]-1);
(*Note the spin has already been flipped.*)

(*Note, below since flip is accepted, the energy of the new state is curE + delE. 
We cannot reassign curE because it is passed by value i.e., is is immutable, hence the use of curE + delE *)

(*(*For diagnostic*)
	If[Abs[data[[Key["energy"]]]-hamiltonian[data[[Key["graph"]]],hparams]]>0.000001,
		Print["For step i ",i," stored Energy ",data[[Key["energy"]]], " computed from the hamiltonian ",
			hamiltonian[data[[Key["graph"]]],hparams]," graph ", data[[Key["graph"]]], 
			" flipped spin ",{row,col}];
	];
(*End diagnostic*)*)


If[curE+delE<minE,minE=curE+delE;minStates={data[[Key["graph"]]]},
If[curE+delE==minE&&Length[minStates]<maxGStateCount,minStates=DeleteDuplicates[Union[minStates,{data[[Key["graph"]]]}]]]];

];
];

Do[step[];

,{i,NN*vCount (vCount-1)/2}];


result={<|"minEnergy"->minE,"minEstates"->minStates|>,data};

Remove[vCount,minE, minStates,maxGStateCount,expDelETable,data,step];


result

];


(*Overload without delH *)
ECGrav`GraphSweepReplica[seedGraph_List,beta_Real,hamiltonian_[hparams___],NN_Integer,minEToBeat_Real,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 11/16/2025  ***)
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

data=<|"graph"->seedGraph,"energy"->hamiltonian[seedGraph,hparams], "mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>;

step[]:=(*Performs one spin flip step*)
Module[{
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

delE=hamiltonian[newAmat,hparams]-curE;


Which[UnlabeledVerticesYes==0,selectionProb=1,UnlabeledVerticesYes==1,
selectionProb=(*Automorphism group order of new graph over auromorphism group order of the old*)
GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[newAmat]]]/GroupOrder[GraphAutomorphismGroup[AdjacencyGraph[data[[Key["graph"]]]]]],True,Abort[];
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


data[[Key["graph"]]]=newAmat;
data[[Key["energy"]]]+= delE;
data[[Key["mag"]]]+= (2.0/(vCount(vCount-1)))*(2*data[[Key["graph"]]][[row,col]]-1);
(*Note the spin has already been flipped.*)

(*(*For diagnostic*)
	If[Abs[data[[Key["energy"]]]-hamiltonian[data[[Key["graph"]]],hparams]]>0.000001,
		Print["stored Energy ",data[[Key["energy"]]], " computed from the hamiltonian ",
			hamiltonian[data[[Key["graph"]]],hparams]," graph ", data[[Key["graph"]]], 
			" flipped spin ",{row,col}];
	];
(*End diagnostic*)
*)
(*(*For diagnostic*)
	If[Abs[curE+delE-hamiltonian[newAmat,hparams]]>0.000001,
		Print["updated E ",curE+delE, " hamiltonian[newAmat] ",
			hamiltonian[newAmat,hparams]," newAmat ", newAmat, 
			" flipped spin ",{row,col}];
	];
(*End diagnostic*)
*)
(*Note, below since flip is accepted, the energy of the new state is curE + delE. We cannot reassign curE because it is passed by value i.e., is is immutable, hence the use of curE + delE *)

If[curE+delE<minE,minE=curE+delE;minStates={data[[Key["graph"]]]},
If[curE+delE==minE&&Length[minStates]<maxGStateCount,minStates=DeleteDuplicates[Union[minStates,{data[[Key["graph"]]]}]]]];

];
];

Do[step[];

,{i,NN*vCount (vCount-1)/2}];


result={<|"minEnergy"->minE,"minEstates"->minStates|>,data};

Remove[vCount,minE, minStates,maxGStateCount,expDelETable,data,step];


result

];


(* ::Item::Closed:: *)
(*GraphSweepReplica Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphSweepReplica[args___[argparams___]]:=(Message[ECGrav`GraphSweepReplica::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*GraphEquilibriate*)


(* ::Item::Closed:: *)
(*GraphEquilibriate *)


ECGrav`GraphEquilibriate[seedGraph_List, beta_Real, hamiltonian_[hparams___],delH_[delHparams___],UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 03/04/2026  ***)
(*************************************)
(*Notes: ,
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
data=<|"minEnergy" ->hamiltonian[seedGraph,hparams],"minEstates"->{seedGraph},
"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph,hparams],"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>,"empty"-><|"graph"->Table[0,{i,vCount},{j,vCount}],"energy" ->0.0,"mag"->0.0|>,
"random"-><|"graph"->With[{rm=RandomInteger[1,{vCount,vCount}]},Mod[rm+Transpose[rm],2]],"energy" ->0.0,"mag"->0.0|>|>;

data[[Key["empty"],Key["energy"]]]=hamiltonian[data[[Key["empty"],Key["graph"]]],hparams];
data[[Key["random"],Key["energy"]]]=hamiltonian[data[[Key["random"],Key["graph"]]],hparams];
data[[Key["random"],Key["mag"]]]=Total[Flatten[data[[Key["random"],Key["graph"]]]]]*1.0/(vCount(vCount-1));


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

Sow[{data[[Key["state"],Key["energy"]]],data[[Key["empty"],Key["energy"]]],data[[Key["random"],Key["energy"]]]}];

Do[


sweepOutput=ECGrav`GraphSweepReplica[data[[Key[replicaName]]][[Key["graph"]]],beta,
				hamiltonian[hparams],delH[delHparams],1,data[[Key["minEnergy"]]],
				UnlabeledVerticesYes];


data[[Key[replicaName]]]=sweepOutput[[2]];

(*Update minimum energy states after the sweep*)
If[sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
	data[[Key["minEnergy"]]]=sweepOutput[[1,Key["minEnergy"]]];
	data[[Key["minEstates"]]]=sweepOutput[[1,Key["minEstates"]]],
	If[sweepOutput[[1,Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<=maxGStateCount,
		data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1,Key["minEstates"]]]]
	];
];
,{replicaName,{"state","empty","random"}}];


Entable[[All,1]]={data[[Key["state"]]][[Key["energy"]]],data[[Key["empty"]]][[Key["energy"]]],data[[Key["random"]]][[Key["energy"]]]};


If[numsweeps>outWinLength,

sqMeanEMat=Table[(Mean[Entable[[i]][[1;;inWinLength]]]-Mean[Entable[[j]][[(outWinLength-inWinLength);;outWinLength]]])^2,{i,3},{j,3}];

(*sqMeanEMat is a three by three matrix of the squared means of difference in energy 
within and across the tracks at the beginning and end of outWinLength.
At equilibrium, these 9 mumbers should be randomly distributed with mean of 0. 
Their fluctuation from 0 should be within the variance of the newer (late time) 
variance in energy for the equilibriation to exit. sqMeanPairwiseDiff is the mean 
of these 9 numbers*)

sqMeanPairwiseDiff=Mean[Flatten[sqMeanEMat]];

meanLateVar=Mean[Table[Variance[Entable[[i]][[1;;inWinLength]]],{i,3}]];(*Mean variance of the newer energies*)


If[Abs[sqMeanPairwiseDiff]<meanLateVar,
	eqlTime=numsweeps-outWinLength+inWinLength;
	Break[]
];(*Exit the while loop and go to the Do loop*)

If[Abs[Tr[sqMeanEMat]]<0.00001,
	(*Print["exiting because stuck in a metastable state, 
		sqMeanEMat is ",MatrixForm[sqMeanEMat]];*)
	eqlTime=numsweeps-outWinLength+inWinLength;
	Break[]
](*Exit the while loop and go to the Do loop*)

];

Entable=RotateRight[#,1]&/@Entable;
(*Shifts every entry to the right by one cyclically. 
The new data will be written on the first slot, so newer data is at the beginning.*)

];
][[2,1]];


(********************
*  For diagnostics *
********************)

(*Print["AllEnMat ",Transpose[AllEnMat]];*)
(*Print["numsweeps ",numsweeps];*)
(*Print["In IsingEquilibriate, eqlTime ",eqlTime, " Length[AllEnMat] ",Length[AllEnMat]];*)

Print[ListLinePlot[Transpose[AllEnMat[[1;;Min[Length[AllEnMat],2*eqlTime]]]],
	PlotRange->All,PlotLabel->"t from 1 to 2 times eqlT for beta "<>ToString[beta]<>
	" hparams { "<>StringJoin[Riffle[ToString/@{hparams},", "]]<>" }"]];

(*Print["data",data];*)

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,
		<|"beta"->beta, "externalField"->{hparams}, "eqlT"->eqlTime,
		"state"->data[[Key["state"]]]|>};

Remove[{vCount,data,maxGStateCount,sweepOutput,Entable,outWinLength,inWinLength,
	AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime,maxNumSweeps}];

result

];


(* Overload without delH *)
ECGrav`GraphEquilibriate[seedGraph_List, beta_Real, hamiltonian_[hparams___],UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 3/04/2026  ***)
(*************************************)
(*Notes: ,
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
data=<|"minEnergy" ->hamiltonian[seedGraph,hparams],"minEstates"->{seedGraph},
"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph,hparams],"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>,"empty"-><|"graph"->Table[0,{i,vCount},{j,vCount}],"energy" ->0.0,"mag"->0.0|>,
"random"-><|"graph"->With[{rm=RandomInteger[1,{vCount,vCount}]},Mod[rm+Transpose[rm],2]],"energy" ->0.0,"mag"->0.0|>|>;

data[[Key["empty"],Key["energy"]]]=hamiltonian[data[[Key["empty"],Key["graph"]]],hparams];
data[[Key["random"],Key["energy"]]]=hamiltonian[data[[Key["random"],Key["graph"]]],hparams];
data[[Key["random"],Key["mag"]]]=Total[Flatten[data[[Key["random"],Key["graph"]]]]]*1.0/(vCount(vCount-1));


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

	Sow[{data[[Key["state"],Key["energy"]]],data[[Key["empty"],Key["energy"]]],
		data[[Key["random"],Key["energy"]]]}];

	Do[


		sweepOutput=ECGrav`GraphSweepReplica[data[[Key[replicaName],Key["graph"]]],
						beta,hamiltonian[hparams],1,data[[Key["minEnergy"]]],
						UnlabeledVerticesYes];


	data[[Key[replicaName]]]=sweepOutput[[2]];

	(*Update minimum energy states after the sweep*)
	If[sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
		data[[Key["minEnergy"]]]=sweepOutput[[1,Key["minEnergy"]]];
		data[[Key["minEstates"]]]=sweepOutput[[1,Key["minEstates"]]],
		If[sweepOutput[[1,Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<=maxGStateCount,
			data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1,Key["minEstates"]]]]
		];
	];
,{replicaName,{"state","empty","random"}}];


Entable[[All,1]]={data[[Key["state"],Key["energy"]]],data[[Key["empty"],Key["energy"]]],
	data[[Key["random"],Key["energy"]]]};


If[numsweeps>outWinLength,
	sqMeanEMat=Table[(Mean[Entable[[i]][[1;;inWinLength]]]-Mean[Entable[[j]][[(outWinLength-inWinLength);;outWinLength]]])^2,{i,3},{j,3}];

	(*sqMeanEMat is a three by three matrix of the squared means of difference in 
	energy within and across the tracks at the beginning and end of outWinLength.
	At equilibrium, these 9 mumbers should be randomly distributed with mean of 0. 
	Their fluctuation from 0 should be within the variance of the newer (late time) 
	variance in energy for the equilibriation to exit. sqMeanPairwiseDiff is the mean 
	of these 9 numbers*)

	sqMeanPairwiseDiff=Mean[Flatten[sqMeanEMat]];

	meanLateVar=Mean[Table[Variance[Entable[[i]][[1;;inWinLength]]],{i,3}]];(*Mean variance of the newer energies*)


	If[Abs[sqMeanPairwiseDiff]<meanLateVar,
	eqlTime=numsweeps-outWinLength+inWinLength;
	Break[]
	];(*Exit the while loop and go to the Do loop*)

	If[Abs[Tr[sqMeanEMat]]<0.00001,
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

Print[ListLinePlot[Transpose[AllEnMat[[1;;Min[Length[AllEnMat],2*eqlTime]]]],
	PlotRange->All,PlotLabel->"t from 1 to 2 times eqlT for beta "<>ToString[beta]<>
	" hparams { "<>StringJoin[Riffle[ToString/@{hparams},", "]]<>" }"]];

(*Print["data",data];*)

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,
		<|"beta"->beta,"externalField"->{hparams},"eqlT"->eqlTime,
			"state"->data[[Key["state"]]]|>};

Remove[{vCount,data,maxGStateCount,sweepOutput,Entable,outWinLength,inWinLength,
	AllEnMat,sqMeanEMat,sqMeanPairwiseDiff,meanLateVar,numsweeps,eqlTime,maxNumSweeps}];

result

];


(* ::Item::Closed:: *)
(*GraphEquilibriate Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphEquilibriate[args___]:=(Message[ECGrav`GraphEquilibriate::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*GraphComputeCorrelationTime*)


(* ::Item::Closed:: *)
(*GraphComputeCorrelationTime Primary*)


ECGrav`GraphComputeCorrelationTime[seedGraph_List,beta_Real,hamiltonian_[hparams___],
	delH_[delHparams___],eqlT_Integer, minEToBeat_Real,EnergyOrMag_Integer,UnlabeledVerticesYes_Integer]:=
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


data=<|"minEnergy" ->minEToBeat,
		"minEstates"->{},
		"beta"->beta, "externalField"->{hparams}, "eqlT"->eqlT, "corrT"->2,
		"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph,hparams],
		"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>|>;


numsweeps=30*eqlT;

PrintTemporary["computing correlation time at beta ",beta, " hparams ",{hparams}, 
	" using Energy or Magnetization ",EorM, " numsweeps ",numsweeps];


empair=Reap[EorMTable=
	Table[
		If[Mod[i,Ceiling[numsweeps/5.0]]==0,PrintTemporary[" sweepno ",i]];
			sweepOutput=ECGrav`GraphSweepReplica[data[[Key["state"]]][[Key["graph"]]],beta,
									hamiltonian[hparams],delH[delHparams],1,
									data[[Key["minEnergy"]]],UnlabeledVerticesYes];


	data[[Key["state"]]]=sweepOutput[[2]];

	(*Update minimum energy states after the sweep*)
	If[
		sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
		data[[Key["minEnergy"]]]=sweepOutput[[1,Key["minEnergy"]]];
		data[[Key["minEstates"]]]=sweepOutput[[1,Key["minEstates"]]],
		If[sweepOutput[[1,Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<maxGStateCount,
			data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1,Key["minEstates"]]]]
		];
	];

	data[[Key["state"],Key[EorM]]]
	,{i,numsweeps}
	];
];

(*For diagnostics*)

(*Table of energy/magnetizations*)


If[Length[DeleteDuplicates[EorMTable]]>1,

(*If the magnetizations have all equal value, then correlation time can not be 
computed, so this loop will be exited with the default corrT left at 2 *)

	norm=ECGrav`CorrelationTime[0,EorMTable];


	If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)


	corrTable=Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,
		PrintTemporary["      computing corrT at t = ",t," hparams ",{hparams}]];
		ECGrav`CorrelationTime[t,EorMTable],{t,0,numsweeps-10}
	]/norm;


tmax=(FirstPosition[corrTable,_?(#<=0&),numsweeps-10])[[1]]; (* A place to stop the integration for calculation of correlation time is when the autocorrelation value first becomes negative*)


Print[ListLinePlot[corrTable[[1;;Min[Length[corrTable],4*tmax]]],PlotRange->Full,
	PlotLabel->"t vs auto correlation for beta "<>ToString[beta]<>" hparams { "<>
	StringJoin[Riffle[ToString/@{hparams},", "]]<>" }"]];

data[[Key["corrT"]]]=
Max[Ceiling[Sum[corrTable[[t]],{t,tmax}]],2]

];

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,
	data[[3;;All]]};

Remove[vCount,data,maxGStateCount,sweepOutput,EorMTable,corrTable,tmax,norm,numsweeps,
	EorM ,empair,maxNumSweeps];

result

];


(* Overload without delH *)
ECGrav`GraphComputeCorrelationTime[seedGraph_List,beta_Real,hamiltonian_[hparams___],eqlT_Integer, 
	minEToBeat_Real,EnergyOrMag_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 3/04/2026  ***)
(*************************************)

(*Notes: ,
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

Module[{result,vCount=Length[seedGraph],data,maxGStateCount,sweepOutput,EorMTable,
	corrTable,tmax,norm,numsweeps,EorM = "energy",empair,maxNumSweeps=5000},

	maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

	Which[EnergyOrMag==0,EorM = "energy",EnergyOrMag==1,EorM = "mag"];


data=<|"minEnergy" ->minEToBeat,
		"minEstates"->{},"beta"->beta,"externalField"->{hparams},
		"eqlT"->eqlT,"corrT"->2,
		"state"-><|"graph"->seedGraph,"energy" ->hamiltonian[seedGraph,hparams],
		"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>|>;


numsweeps=30*eqlT;

PrintTemporary["computing correlation time at beta ",beta, " hparams ",{hparams}, 
	" using Energy or Magnetization ",EorM, " numsweeps ",numsweeps];


empair=Reap[EorMTable=
Table[
If[Mod[i,Ceiling[numsweeps/5.0]]==0,PrintTemporary[" sweepno ",i]];
sweepOutput=ECGrav`GraphSweepReplica[data[[Key["state"]]][[Key["graph"]]],beta,
									hamiltonian[hparams],1,data[[Key["minEnergy"]]],
									UnlabeledVerticesYes];


data[[Key["state"]]]=sweepOutput[[2]];

(*Update minimum energy states after the sweep*)
If[
	sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
	data[[Key["minEnergy"]]]=sweepOutput[[1,Key["minEnergy"]]];
	data[[Key["minEstates"]]]=sweepOutput[[1,Key["minEstates"]]],
	If[sweepOutput[[1,Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<maxGStateCount,
		data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1,Key["minEstates"]]]]
	];
];

(*Sow[{hamiltonian[data[[Key["state"]]][[Key["graph"]]],hparams],
Total[Flatten[data[[Key["state"]]][[Key["graph"]]]]]*1.0/(vCount(vCount-1))}];*)

data[[Key["state"]]][[Key[EorM]]]
,{i,numsweeps}
];
];

(*For diagnostics*)

(*Table of energy/magnetizations*)


If[Length[DeleteDuplicates[EorMTable]]>1,

(*If the magnetizations have all equal value, then correlation time can not be computed, so this loop will be exited with the default corrT left at 2 *)

norm=ECGrav`CorrelationTime[0,EorMTable];


If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)


corrTable=Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,
	PrintTemporary["      computing corrT at t = ",t," hparams ",{hparams}]];
	ECGrav`CorrelationTime[t,EorMTable],{t,0,numsweeps-10}
]/norm;


tmax=(FirstPosition[corrTable,_?(#<=0&),numsweeps-10])[[1]]; (* A place to stop the integration for calculation of correlation time is when the autocorrelation value first becomes negative*)

Print[ListLinePlot[corrTable[[1;;Min[Length[corrTable],4*tmax]]],PlotRange->Full,
	PlotLabel->"t vs auto correlation for beta "<>ToString[beta]<>" hparams { "<>
	StringJoin[Riffle[ToString/@{hparams},", "]]<>" }"]];

data[[Key["corrT"]]]=
Max[Ceiling[Sum[corrTable[[t]],{t,tmax}]],2]

];

result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,data[[3;;All]]};

Remove[vCount,data,maxGStateCount,sweepOutput,EorMTable,corrTable,tmax,norm,numsweeps,EorM ,empair,maxNumSweeps];

result

];


(* ::Item::Closed:: *)
(*GraphComputeCorrelationTime with operators*)


ECGrav`GraphComputeCorrelationTime[seedGraph_List,beta_Real,hamiltonian_[hparams___],
	delH_[delHparams___],eqlT_Integer, minEToBeat_Real,operators_/;MatchQ[operators,{__Function}],
	UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 3/31/2026  ***)
(*************************************)

(*Notes: *)
(*
Takes an equilibriated graph model with equilibriation time, graph configurations, beta, and minEToBeat and computes the correlation time., 

Inputs are:,
1. seedGraph = List, adjacency matrix of the seed graph,
2. beta = Real, inverse temperature,
3. hamiltonian = formula, the hamiltonian, 
4. delH = a formula for delta E (when one edge is flipped to expedite computation),  
5. eqlT = equilibriation time.,6. minEToBeat = energy value such that if a state (or states) are found with energy lower than this then the lowest such energy and configurations with that energy will be saved.,
6. minEToBeat = energy value such that if a state (or states) are found with energy lower than this then the lowest such energy and configurations with that energy will be saved.,
7. operators = functions whose correlation time is to be computed,
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

Module[{result,vCount=Length[seedGraph],data,maxGStateCount,sweepOutput,observablesTable,
	fluctuatingObservableIndices,norm,corrTable,tmaxVals,corrTValues,numsweeps},

	maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

	data=<|"minEnergy" ->minEToBeat,"minEstates"->{},"beta"->beta,"externalField"->{hparams},
			"eqlT"->eqlT,"corrT"->2,"corrTValues"->Table[2,{Length[operators]+2}],"state"-><|"graph"->seedGraph,
			"energy" ->hamiltonian[seedGraph,hparams],
			"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>|>;

	numsweeps=30*eqlT;

	PrintTemporary["computing correlation time at beta ",beta," hparams ",{hparams}, " numsweeps ",numsweeps];

	
	observablesTable=
	Transpose[
		Table[
			If[Mod[i,Ceiling[numsweeps/5.0]]==0,PrintTemporary[" sweepno ",i]];
				sweepOutput=ECGrav`GraphSweepReplica[data[[Key["state"],Key["graph"]]],beta,
								hamiltonian[hparams],delH[delHparams],1,data[[Key["minEnergy"]]],
								UnlabeledVerticesYes];

			data[[Key["state"]]]=sweepOutput[[2]];

			(*Update minimum energy states after the sweep*)
			If[
				sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
				data[[Key["minEnergy"]]]=sweepOutput[[1]][[Key["minEnergy"]]];
				data[[Key["minEstates"]]]=sweepOutput[[1]][[Key["minEstates"]]],
				If[sweepOutput[[1]][[Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<maxGStateCount,
					data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1]][[Key["minEstates"]]]]
				];
			];


			Flatten[{data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["mag"]]],
					Through[operators[data[[Key["state"],Key["graph"]]]]]}]
		,{i,1,numsweeps}]
	];


	(*If the measured observable values have all equal value, then correlation time can not be computed, 
		so fluctuating selects the measurements whose correlation time can be computed. If 
		none of the measurements fluctuate and are all stuck at a single value, the program
		will return the default values corrT -> 2 and corrTValues->{} *)
	fluctuatingObservableIndices=
	Table[
		If[Length[DeleteDuplicates[observablesTable[[i]]]]>1,i,
			Which[
				i==1,
				Message[ECGrav`GraphComputeCorrelationTime::stuck, "the energy", observablesTable[[1,1]]],
				i==2,
				Message[ECGrav`GraphComputeCorrelationTime::stuck, "the magnetization", observablesTable[[2,1]]],
				i>2,
				Message[ECGrav`GraphComputeCorrelationTime::stuck, operators[[i-2]], observablesTable[[i-2,1]]]];
		Nothing],
	{i,1,Length[observablesTable]}];


	If[fluctuatingObservableIndices=={},
		Message[ECGrav`GraphComputeCorrelationTime::alldefault,
			DeleteDuplicates[#]&/@observablesTable,data[[Key["corrTValues"]]]],
		corrTable=Table[

			norm=ECGrav`CorrelationTime[0,observablesTable[[i]]];
			If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)

			Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,
				PrintTemporary["      computing corrT at t = ",t]];
				ECGrav`CorrelationTime[t,observablesTable[[i]]],
			{t,0,numsweeps-10}]/norm


		,{i,fluctuatingObservableIndices}];
	
	
		tmaxVals=Table[(FirstPosition[cT,_?(#<=0&),numsweeps-10])[[1]],{cT,corrTable}]; 
		(* A place to stop the integration for calculation of correlation time is when the 
		autocorrelation value first becomes negative*)
	
	
		Print[
			ListLinePlot[
				Table[corrTable[[i,1;;Min[Length[corrTable[[i]]],4*tmaxVals[[i]]]]]
				,{i,1,Length[corrTable]}],
				PlotLegends->Flatten[{"energy","mag",ToString/@operators}],
				PlotRange->Full,
				PlotLabel->"t vs auto correlation for beta "<>ToString[beta]<>" hparams "<>ToString[{hparams}]
			]
		];

		corrTValues=Table[Max[Ceiling[Sum[corrTable[[i,t]],{t,1,tmaxVals[[i]]}]],2],{i,1,Length[corrTable]}];

		
		data[[Key["corrT"]]]=Max[corrTValues];
		Do[data[[Key["corrTValues"],j]]=corrTValues[[First@Flatten[Position[fluctuatingObservableIndices,j]]]],{j,fluctuatingObservableIndices}];
		(*data[[Key["corrTValues"]]]=Table[If[MemberQ[i,]data[[Key["corrTValues"],i]],corrTValues,{i,1,Length[data[[Key["corrTValues"]]]]}]];*)

	];

	result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,data[[3;;All]]};

	Remove[vCount,data,maxGStateCount,sweepOutput,observablesTable,fluctuatingObservableIndices,norm,corrTable,tmaxVals,corrTValues
		,numsweeps];

	result

];


(* Overload without delH *)

ECGrav`GraphComputeCorrelationTime[seedGraph_List,beta_Real,hamiltonian_[hparams___],
	eqlT_Integer, minEToBeat_Real,operators_/;MatchQ[operators,{__Function}],
	UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 3/31/2026  ***)
(*************************************)

(*Notes: *)
(*
Takes an equilibriated graph model with equilibriation time, graph configurations, beta, and minEToBeat and computes the correlation time., 

Inputs are:,
1. seedGraph = List, adjacency matrix of the seed graph,
2. beta = Real, inverse temperature,
3. hamiltonian = formula, the hamiltonian,  
4. eqlT = equilibriation time.,6. minEToBeat = energy value such that if a state (or states) are found with energy lower than this then the lowest such energy and configurations with that energy will be saved.,
5. minEToBeat = energy value such that if a state (or states) are found with energy lower than this then the lowest such energy and configurations with that energy will be saved.,
6. operators = functions whose correlation time is to be computed,
7. UnlabeledVerticesYes = Integer:  0 means no selection probability to make the graphs unlabeled so graphs are labeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with two associations,
1. the minimum energy visited throughout the equilibriation and the states with that energy. If multiple states have degenerate minimum energy, they will all be included; but if multiple identical adjacency graphs are found, only one unique adjacency graph is kept. ,
2. The second association has the inverse temperature, equilibriation time, 
    correlation timeand the final state visited,   which itself is an association which includes the adjacency matrix, 
   magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.,

Depends on the functions: GraphSweepReplicas, CorrelationTime., 

 *)

Module[{result,vCount=Length[seedGraph],data,maxGStateCount,sweepOutput,observablesTable,
	norm,corrTable,tmaxVals,corrTValues,numsweeps},

	maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

	data=<|"minEnergy" ->minEToBeat,"minEstates"->{},"beta"->beta,"externalField"->{hparams},
			"eqlT"->eqlT,"corrT"->2,"corrTValues"->{},"state"-><|"graph"->seedGraph,
			"energy" ->hamiltonian[seedGraph,hparams],
			"mag"->Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>|>;

	numsweeps=30*eqlT;

	PrintTemporary["computing correlation time at beta ",beta," hparams ",{hparams}, " numsweeps ",numsweeps];

	observablesTable=
	Transpose[
		Table[
			If[Mod[i,Ceiling[numsweeps/5.0]]==0,PrintTemporary[" sweepno ",i]];
				sweepOutput=ECGrav`GraphSweepReplica[data[[Key["state"],Key["graph"]]],beta,
								hamiltonian[hparams],1,data[[Key["minEnergy"]]],
								UnlabeledVerticesYes];

			data[[Key["state"]]]=sweepOutput[[2]];

			(*Update minimum energy states after the sweep*)
			If[
				sweepOutput[[1,Key["minEnergy"]]]<data[[Key["minEnergy"]]],
				data[[Key["minEnergy"]]]=sweepOutput[[1]][[Key["minEnergy"]]];
				data[[Key["minEstates"]]]=sweepOutput[[1]][[Key["minEstates"]]],
				If[sweepOutput[[1]][[Key["minEnergy"]]]==data[[Key["minEnergy"]]]&&Length[data[[Key["minEstates"]]]]<maxGStateCount,
					data[[Key["minEstates"]]]=Union[data[[Key["minEstates"]]],sweepOutput[[1]][[Key["minEstates"]]]]
				];
			];


			Flatten[{data[[Key["state"],Key["energy"]]],data[[Key["state"],Key["mag"]]],
					Through[operators[data[[Key["state"],Key["graph"]]]]]}]
		,{i,1,numsweeps}]
	];


	corrTable=Table[
		If[Length[DeleteDuplicates[oTbl]]>1,

		(*If the magnetizations have all equal value, then correlation time can not be computed, 
		so this loop will be exited with the default corrT left at 2 *)

		norm=ECGrav`CorrelationTime[0,oTbl];
		If[norm==0.0,norm=1.0]; (*the time 0 correlation can be 0 sometimes*)

		Table[If[Mod[t,Ceiling[numsweeps/5.0]]==0,
			PrintTemporary["      computing corrT at t = ",t]];
			ECGrav`CorrelationTime[t,oTbl],
		{t,0,numsweeps-10}]/norm

		]

	,{oTbl,observablesTable}];

	tmaxVals=Table[(FirstPosition[cT,_?(#<=0&),numsweeps-10])[[1]],{cT,corrTable}]; 
	(* A place to stop the integration for calculation of correlation time is when the 
	autocorrelation value first becomes negative*)

	
	corrTValues=Table[Max[Ceiling[Sum[corrTable[[i,t]],{t,1,tmaxVals[[i]]}]],2],{i,1,Length[corrTable]}];

	data[[Key["corrT"]]]=Max[corrTValues];

	data[[Key["corrTValues"]]]=corrTValues;

	result={<|"minEnergy"->data[[Key["minEnergy"]]],"minEstates"->data[[Key["minEstates"]]]|>,data[[3;;All]]};

	Remove[vCount,data,maxGStateCount,sweepOutput,observablesTable,norm,corrTable,tmaxVals,corrTValues
		,numsweeps];

	result

];


(* ::Item::Closed:: *)
(*GraphComputeCorrelationTime Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphComputeCorrelationTime[args___]:=(Message[ECGrav`GraphComputeCorrelationTime::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*GraphMultiHistogram*)


(* ::Item::Closed:: *)
(*GraphMultiHistogram Primary betaMin to betaHigh*)


ECGrav`GraphMultiHistogram[seedGraph_List,betaLow_Real,betaHigh_Real,
	hamiltonian_[hparams___],delH_[delHparams___],obs_/;MatchQ[obs,{__Function}],NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 10/10/2025  ***)
(*************************************)

(*Implements the Multiple Histogram Method for the graph models to get a smooth plot 
of the quantity obs as a function of inverse temperature ranging from betaLow to 
betaHigh. This does not have replica swaps.,

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
groundStates=<|"minEnergy"->hamiltonian[seedGraph,hparams],"minEstates"->{seedGraph}|>;


maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

btTable={{betaLow,betaHigh}};

defaultRatio={Exp[Log[betaHigh/betaLow]/(vCount*1.0)],Exp[Log[betaLow/betaHigh]/(vCount*1.0)]};
(*A default value for increment of beta if the specific heat happens to be 0. In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective, hence the factor of vCount.*)


stopnum=1;

While[stopnum<500,

stopnum++;


(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[
		i->ECGrav`GraphEquilibriate[seedGraph,i,hamiltonian[hparams],
									delH[delHparams],UnlabeledVerticesYes]
		,{i,btTable[[-1]]},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Union[replicas,Tempoutput[[All,2]]];

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[
			i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]],i,
								hamiltonian[hparams],delH[delHparams],locrepl[[Key[i],Key["eqlT"]]],
								locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes]
			,{i,btTable[[-1]]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,btTable[[-1]]}];


(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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
If[Mod[numsweeps,Ceiling[NN/5.0]]==0,PrintTemporary[" sweepno ",numsweeps]];


Tempoutput=Association[With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
	ParallelTable[
	(*repNumSweeps=replicas[[Key[i]]][[Key["corrT"]]];*)
	(* Each replica will be swept corrT times *)


	<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],i,
			hamiltonian[hparams],delH[delHparams],locrepl[[Key[i],Key["corrT"]]],
			locMinEtoBeat,UnlabeledVerticesYes]|>,

	{i,btTable[[-1]]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*Update replicas*)
Do[replicas[[Key[i]]][[Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,btTable[[-1]]}];

(* Extract the minimum energy and corresponding states *)

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

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


If[Mod[numsweeps,1]==0,
Table[

Sow[Flatten[{numsweeps,i,replicas[[Key[i]]][[Key["state"]]][[Key["energy"]]],Through[obs[replicas[[Key[i]]][[Key["state"]]][[Key["graph"]]]]]}],i]
,{i,btTable[[-1]]}];
];

]][[2]];


chart=Union[chart,AssociationThread[btTable[[-1]],measurements[[1;;2]]]];


If[(btTable[[-1,1]]>btTable[[-1,2]]),Break[]];

(* sqrt of the specific heat at current temp*)
curRootSpecificHeat=Table[Sqrt[ECGrav`SpecificHeat[chart[[Key[i]]][[All,3]],1,i]],{i,btTable[[-1]]}];
(*Total specific heat, not specific heat per site!*)


curRatios=
Table[If[curRootSpecificHeat[[k]]!=0,(1.0+1.0/curRootSpecificHeat[[k]])^k,defaultRatio[[k]]],{k,{1,-1}}];
(*when k=1 it multiplies, and divides when k=2. Ensures that the first factor is greater than 1, and the second less than one.*)

newBetas={btTable[[-1,1]]*curRatios[[1]],btTable[[-1,2]]*curRatios[[2]]};


AppendTo[btTable,newBetas];


];


btTable=Sort[Flatten[btTable]];


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


(*mBetaF=computeBFs[energyAsn];*)

mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[energyAsn];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas};

Remove[vCount,edgeCount,groundStates,maxGStateCount,replicas,Tempoutput,btTable,defaultRatio,measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,curRootSpecificHeat,curRatios,newBetas,computeBFs,energyAsn,mBetaF];

result

];


(*Overload betaLow to betaHigh without delH*)
ECGrav`GraphMultiHistogram[seedGraph_List,betaLow_Real,betaHigh_Real,
	hamiltonian_[hparams___],
	obs_/;MatchQ[obs,{__Function}],NN_Integer,UnlabeledVerticesYes_Integer]:=
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

Module[{result,vCount=Length[seedGraph],edgeCount,groundStates,maxGStateCount,replicas,
	Tempoutput,btTable,defaultRatio,
	EnergyOrMag=1 (*comuting correlation time using energy (0) or magnetization (1)*),
	measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,curRootSpecificHeat,
	curRatios,newBetas,computeBFs,mBetaF},

edgeCount=vCount (vCount-1)/2;

replicas=<||>;
chart=<||>;
groundStates=<|"minEnergy"->hamiltonian[seedGraph,hparams],"minEstates"->{seedGraph}|>;


maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

btTable={{betaLow,betaHigh}};

defaultRatio={Exp[Log[betaHigh/betaLow]/(vCount*1.0)],Exp[Log[betaLow/betaHigh]/(vCount*1.0)]};
(*A default value for increment of beta if the specific heat happens to be 0. In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective, hence the factor of vCount.*)


stopnum=1;

While[stopnum<500,

stopnum++;


(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[
		i->ECGrav`GraphEquilibriate[seedGraph,i,hamiltonian[hparams],
				UnlabeledVerticesYes]
		,{i,btTable[[-1]]},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Union[replicas,Tempoutput[[All,2]]];

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[
			i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]],i,
					hamiltonian[hparams],locrepl[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes]
			,{i,btTable[[-1]]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,btTable[[-1]]}];


(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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
If[Mod[numsweeps,Ceiling[NN/5.0]]==0,PrintTemporary[" sweepno ",numsweeps]];


Tempoutput=Association[With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
	ParallelTable[
	(*repNumSweeps=replicas[[Key[i]]][[Key["corrT"]]];*)
	(* Each replica will be swept corrT times *)


	<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],i,
		hamiltonian[hparams],locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

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


If[Mod[numsweeps,1]==0,
Table[

Sow[Flatten[{numsweeps,i,replicas[[Key[i]]][[Key["state"]]][[Key["energy"]]],Through[obs[replicas[[Key[i]]][[Key["state"]]][[Key["graph"]]]]]}],i]
,{i,btTable[[-1]]}];
];

]][[2]];


(*Print[" btTable ",btTable ];
Print[" measurements ",measurements ];*)

chart=Union[chart,AssociationThread[btTable[[-1]],measurements[[1;;2]]]];


If[(btTable[[-1,1]]>btTable[[-1,2]]),Break[]];

(* sqrt of the specific heat at current temp*)
curRootSpecificHeat=Table[Sqrt[ECGrav`SpecificHeat[chart[[Key[i]]][[All,3]],1,i]],{i,btTable[[-1]]}];
(*Total specific heat, not specific heat per site!*)


curRatios=
Table[If[curRootSpecificHeat[[k]]!=0,(1.0+1.0/curRootSpecificHeat[[k]])^k,defaultRatio[[k]]],{k,{1,-1}}];
(*when k=1 it multiplies, and divides when k=2. Ensures that the first factor is greater than 1, and the second less than one.*)


newBetas={btTable[[-1,1]]*curRatios[[1]],btTable[[-1,2]]*curRatios[[2]]};


AppendTo[btTable,newBetas];


];


btTable=Sort[Flatten[btTable]];


(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each temp, i.e., (-beta*F). It is the
* same as computing the partition function, *   
* but better due to overflow issues   *
*        *)
(**********************************************)
*)


mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[chart[[All,All,3]]];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas};

Remove[vCount,edgeCount,groundStates,maxGStateCount,replicas,Tempoutput,btTable,defaultRatio,measurements,numsweeps,stopnum, candminE,repNumSweeps,chart,curRootSpecificHeat,curRatios,newBetas,computeBFs,mBetaF];

result

];


(* ::Item::Closed:: *)
(*GraphMultiHistogram Overload 1 betaTable*)


ECGrav`GraphMultiHistogram[seedGraph_List,btTable_List/;VectorQ[btTable],hamiltonian_[hparams___],
	delH_[delHparams___],obs_/;MatchQ[obs,{__Function}],NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 12/02/2026  ***)
(*************************************)

(*An overload which implements the Multiple Histogram Method for the graph models 
to get a smooth plot of the quantity obs as a function of inverse temperatures 
given by btTable. It applies replica swaps during the last measurement stage.,

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
3. the replicas at the last step,
4. an association of the histories of each replica *)

Module[{result,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,numRep,histories,bt,Tempoutput,swap,EnergyOrMag=1 (*comuting correlation time using energy (0) or magnetization (1)*),measurements,numsweeps,stopnum, candminE,printCase,repNumSweeps,chart,computeBFs,energyAsn,mBetaF},


groundStates=<||>;
histories=<||>;
numRep=Length[btTable];

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)


(*In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective.*)


(*
(********************)
(*   Equilibriate   *)
(********************)
*)

Tempoutput=Association[ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,btTable[[i]],hamiltonian[hparams],delH[delHparams],UnlabeledVerticesYes],{i,Length[btTable]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

(* Extract the minimum energy and corresponding states *)

AppendTo[groundStates,"minEnergy"->Min[Tempoutput[[All,1,"minEnergy"]]]];

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


AppendTo[groundStates,"minEstates"->Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==groundStates[[Key["minEnergy"]]]&][[All,"minEstates"]]]];

(*compare the minimum energy from the equilibriation run with that of the temperature schedule run and reset the minimum energy and states *)


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]],replicas[[Key[i],Key["beta"]]],hamiltonian[hparams],delH[delHparams],replicas[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes],{i,numRep},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
]
];


(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];


(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;


Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Module[{thisReplInd,nextReplInd,deltaBetaDeltaE,accept,tempThisReplicaState},

(*Print[""];
Print[""];
Print[""];
Print[""];
Print["Starting SwapReplicas with offset ",offset];

Print[" replicas ",replicas];
Print[" histories ",histories];*)


Do[
thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];


deltaBetaDeltaE= (replicas[[Key[thisReplInd],Key["beta"]]]-replicas[[Key[nextReplInd],Key["beta"]]])*(replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]-replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]);


(*Increase the number of swap tries by one for both replicas*)

histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;
histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;


(*Print[""];
Print[""];
Print["Before swap, replicas ",replicas];
Print["Before swap, histories ",histories];*)


accept=0;
If[deltaBetaDeltaE>0,accept = 1,
If[RandomReal[]<Exp[deltaBetaDeltaE],accept =1];
];


If[accept==1, (*Do swap of replicas*)


(*Update states: store the state of the high temp replica temporarily*)

tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

(*Print["tempThisReplicaState ",tempThisReplicaState];
Print["replicas ",replicas];*)


(*Update histories*)

histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

(*Increase the number of swap accepts by one for both replicas*)
histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;


];


(*Print["  After one swap, replicas ",replicas];
Print["  histories ",histories];*)

,{repl,1,numRep-1,2}];


];


(*(******************************************
* Take NN measurements with swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
While[numsweeps<NN,


If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];


Tempoutput=Association[With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
ParallelTable[

<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],locrepl[[Key[i],Key["beta"]]],hamiltonian[hparams],delH[delHparams],locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

{i,numRep},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
]
];

(*Print["Before breaking : "];
Print["replicas ",replicas];
Print["gstates ",groundStates];*)

(*Update replicas*)
Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

(* Extract the minimum energy and corresponding states *)

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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


swap[numsweeps];


If[Mod[numsweeps,1]==0,
Table[
bt=replicas[[Key[i],Key["beta"]]];Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
,{i,numRep}];
];

numsweeps++;

]][[2]];


(*Print[" btTable ",btTable ];
Print[" measurements ",measurements ];*)

(*chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];*)


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[btTable]]]}]|>;


Do[histories[[Key[i],Key["history"]]]=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]]),{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)


(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each temp, i.e., (-beta*F). It is the
* same as computing the partition function, *   
* but better due to overfloq issues   *
*        *)
(**********************************************)
*)


mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[chart[[All,All,3]]];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"beta"->i,"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas,histories};


Remove[groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,numRep,histories,bt,Tempoutput,swap,EnergyOrMag,measurements,numsweeps,stopnum, candminE,printCase,repNumSweeps,chart,computeBFs,energyAsn,mBetaF];

result

];


(*Overload without delH *)
ECGrav`GraphMultiHistogram[seedGraph_List,btTable_List/;VectorQ[btTable],hamiltonian_[hparams___],
obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 12/02/2026  ***)
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
3. the replicas at the last step,
4. an association of the histories of each replica *)

Module[{result,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,numRep,
	histories,bt,Tempoutput,swap,EnergyOrMag=1 (*comuting correlation time using energy (0) or magnetization (1)*),measurements,numsweeps,stopnum, candminE,printCase,repNumSweeps,chart,computeBFs,energyAsn,mBetaF},


groundStates=<||>;
histories=<||>;
numRep=Length[btTable];

maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)


(*In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective.*)


(*
(********************)
(*   Equilibriate   *)
(********************)
*)

Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,btTable[[i]],hamiltonian[hparams],
		UnlabeledVerticesYes],{i,Length[btTable]},
		DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

(* Extract the minimum energy and corresponding states *)

AppendTo[groundStates,"minEnergy"->Min[Tempoutput[[All,1,"minEnergy"]]]];

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


AppendTo[groundStates,"minEstates"->Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==groundStates[[Key["minEnergy"]]]&][[All,"minEstates"]]]];

(*compare the minimum energy from the equilibriation run with that of the temperature schedule run and reset the minimum energy and states *)


(*
(***************************************)
(*  Compute energy correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],
			Key["graph"]]],replicas[[Key[i],Key["beta"]]],hamiltonian[hparams],
			replicas[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,
			UnlabeledVerticesYes],{i,numRep},
			DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];


(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];


(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;


Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Module[{thisReplInd,nextReplInd,deltaBetaDeltaE,accept,tempThisReplicaState},

(*Print[""];
Print[""];
Print[""];
Print[""];
Print["Starting SwapReplicas with offset ",offset];

Print[" replicas ",replicas];
Print[" histories ",histories];*)


Do[
thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];


deltaBetaDeltaE= (replicas[[Key[thisReplInd],Key["beta"]]]-replicas[[Key[nextReplInd],Key["beta"]]])*(replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]-replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]);


(*Increase the number of swap tries by one for both replicas*)

histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;
histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;


(*Print[""];
Print[""];
Print["Before swap, replicas ",replicas];
Print["Before swap, histories ",histories];*)


accept=0;
If[deltaBetaDeltaE>0,accept = 1,
If[RandomReal[]<Exp[deltaBetaDeltaE],accept =1];
];


If[accept==1, (*Do swap of replicas*)


(*Update states: store the state of the high temp replica temporarily*)

tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

(*Print["tempThisReplicaState ",tempThisReplicaState];
Print["replicas ",replicas];*)


(*Update histories*)

histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

(*Increase the number of swap accepts by one for both replicas*)
histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;


];


(*Print["  After one swap, replicas ",replicas];
Print["  histories ",histories];*)

,{repl,1,numRep-1,2}];


];


(*(******************************************
* Take NN measurements with swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
While[numsweeps<NN,


If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];


Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
	ParallelTable[
	<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]],
		locrepl[[Key[i],Key["beta"]]],hamiltonian[hparams],locrepl[[Key[i],
		Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

	{i,numRep},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
	]
];

(*Print["Before breaking : "];
Print["replicas ",replicas];
Print["gstates ",groundStates];*)

(*Update replicas*)
Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

(* Extract the minimum energy and corresponding states *)

candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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


swap[numsweeps];


If[Mod[numsweeps,1]==0,
	Table[
	bt=replicas[[Key[i],Key["beta"]]];
	Sow[
		Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],
		Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
	,{i,numRep}];
];

numsweeps++;

]][[2]];


(*Print[" btTable ",btTable ];
Print[" measurements ",measurements ];*)

(*chart=AssociationThread[btTable,measurements[[1;;Length[btTable]]]];*)


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[btTable]]]}]|>;


Do[histories[[Key[i],Key["history"]]]=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]]),{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)


(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each temp, i.e., (-beta*F). It is the
* same as computing the partition function, *   
* but better due to overfloq issues   *
*        *)
(**********************************************)
*)


mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[chart[[All,All,3]]];

result={groundStates,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"beta"->i,"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>,replicas,histories};


Remove[groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,numRep,histories,
	bt,Tempoutput,swap,EnergyOrMag,measurements,numsweeps,stopnum, candminE,printCase,
	repNumSweeps,chart,computeBFs,energyAsn,mBetaF];

result

];


(* ::Item::Closed:: *)
(*GraphMultiHistogram Overload 2 one fixed beta value, multiple external fields given as a matrix*)


(*Overload with multiple external fields*)
ECGrav`GraphMultiHistogram[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	delH_[delHparams___],externalFieldTable_List/;MatrixQ[externalFieldTable],
	conjugateObs_/;MatchQ[conjugateObs,{__Function}],obs_/;MatchQ[obs,{__Function}],NN_Integer,
	UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 4/03/2026  ***)
(*************************************)

(*An overload which implements the Multiple Histogram Method for the graph models to 
get a smooth plot of the quantities obs as a function of multiple external parameters at a 
fixed inverse temperature. It applies replica swaps during the measurement stage.,

(*Notes: ,

*)

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. bt - Real = a fixed value of the inverse temperature,
3. hamiltonian[hparams] =  the hamiltonian function with its parameters hparams that 
	assigns graphs energy,
4. delH[delHparams] = a delta-E function that gives the change in energy when a 
	single edge is flipped with its parameters delHparams,
5. externalFieldTable - List = a list of list of external field values {{J^1_1,J^2_1},{J^1_2,J^2_2}...},
6. conjugateObs_ = a list of  functions O which are conjugate to the external fields 
	J^i i.e., there is a term in the microscopic hamiltonian of the form Sum_i(J^i*O^i(s)),
7. obs = a list of various observables,
8. NN = number of independent sweeps (so that actual number of sweeps is correlation 
	time times NN),
9. UnlabeledVerticesYes = 0 means no selection probability to make the graphs 
	unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with three entries:,
1. an association of the minimum energy found from the run and the states found 
	having that energy.,
2. an association with external field table values as keys and values of 
	negative*beta*free energy and values for the energies and observables at each 
	external field value,
3. the replicas at the last step,
4. an association of the histories of each replica *)

Module[{result,groundStates,maxGStateCount,replicas,replicaKeysOrderedExternalField,numRep,
	histories,swap,Tempoutput,measurements,numsweeps,stopnum, candminE,printCase,repNumSweeps,
	chart,mBetaF},


groundStates=<||>;
histories=<||>;
numRep=Length[externalFieldTable];


maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

(*In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective.*)

PrintTemporary[" Starting multihistogram with replica swap and delH and at beta "
	,bt," external field parameters ",externalFieldTable];


(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[
		i->ECGrav`GraphEquilibriate[seedGraph,bt
					,Apply[hamiltonian,externalFieldTable[[i]]]
					,Apply[delH,externalFieldTable[[i]]],UnlabeledVerticesYes]
	,{i,Length[externalFieldTable]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

(* Extract the minimum energy and corresponding states *)

groundStates=Tempoutput[[All,1]];

(*
(***************************************)
(*  Compute correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	ParallelTable[
		i->ECGrav`GraphComputeCorrelationTime[
			replicas[[Key[i],Key["state"],Key["graph"]]],bt
			,Apply[hamiltonian,externalFieldTable[[i]]]
			,Apply[delH,externalFieldTable[[i]]],replicas[[Key[i],Key["eqlT"]]]
			,groundStates[[Key[i],Key["minEnergy"]]],conjugateObs,UnlabeledVerticesYes]
		,{i,numRep},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}
	]
];


(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];


(*Update ground states*)

Do[

	(* Extract the minimum energy and corresponding states *)
	candminE=Tempoutput[[Key[i],1,"minEnergy"]];


	If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
		groundStates[[Key[i],Key["minEnergy"]]]=candminE;
		groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
		If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
			groundStates[[Key[i],Key["minEstates"]]]=
				Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
		];
	];

,{i,numRep}];

replicaKeysOrderedExternalField=Keys[Sort[replicas[[All,"externalField"]]]];

histories=<|Table[i-><|"externalField"->replicas[[Key[i],Key["externalField"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["externalField"]]],{i,numRep}];

(*(*********************)
(*  The swap function  *)
(*********************)*)

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-beta(delta )(delta E)). *)
Module[{thisReplInd,nextReplInd,betaDelHDelConj,accept,tempThisReplicaState},

Do[
	thisReplInd=replicaKeysOrderedExternalField[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedExternalField[[Mod[offset+repl,numRep]+1]];


	betaDelHDelConj=
		bt*((replicas[[Key[thisReplInd],Key["externalField"]]]
				-replicas[[Key[nextReplInd],Key["externalField"]]])
			. (Through[conjugateObs[replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]]]
				-Through[conjugateObs[replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]]])
			);


	histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;
	histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;

	accept=0;
	If[betaDelHDelConj>0,accept = 1,
		If[RandomReal[]<Exp[betaDelHDelConj],accept =1];
	];


	If[accept==1, (*Do swap of replicas*)

		(*Update states: store the state of the high temp replica temporarily*)

		tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
		replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
		replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

		(*Update the energies to reflect the new external field values*)
	
		replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]
			=Apply[hamiltonian
				,Join[{replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]}
					,replicas[[Key[thisReplInd],Key["externalField"]]]]];
		replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]
			=Apply[hamiltonian
				,Join[{replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]}
					,replicas[[Key[nextReplInd],Key["externalField"]]]]];


		(*Update histories*)

		histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
		histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

		histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
		histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

		(*Increase the number of swap accepts by one for both replicas*)
		histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
		histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;


	];


,{repl,1,numRep-1,2}];

];


(*
(****************************************)
(* Take NN measurements at the           *
* equilibriated temperatures to compute *
* specific heats *)
(****************************************)
*)

numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];

measurements=Reap[
While[numsweeps<=NN,

	If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

	Tempoutput=Association[
		ParallelTable[
			(* Each replica will be swept corrT times *)
			<|i->ECGrav`GraphSweepReplica[replicas[[Key[i],Key["state"],Key["graph"]]]
				,bt,Apply[hamiltonian,externalFieldTable[[i]]],Apply[delH,externalFieldTable[[i]]]
				,replicas[[Key[i],Key["corrT"]]],groundStates[[Key[i],Key["minEnergy"]]]
				,UnlabeledVerticesYes]|>,

			{i,numRep},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}
			]
	];


	(*Update replicas*)
	Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i],2]],{i,numRep}];

	(*Update ground states*)

	Do[

		(* Extract the minimum energy and corresponding states *)
		candminE=Tempoutput[[Key[i],1,"minEnergy"]];


		If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
			groundStates[[Key[i],Key["minEnergy"]]]=candminE;
			groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
			If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
				groundStates[[Key[i],Key["minEstates"]]]=
					Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
			];
		];

	,{i,numRep}];

	(*(*For diagnostic*)
	Print[""];
	Print[""];
	Print[""];
	Print[""];Print[""];Print[""];
	Print[" Energy check before swap "];
	Table[
		If[Abs[hamiltonian[replicas[[Key[i],Key["state"],Key["graph"]]],replicas[[Key[i],Key["externalField"]]]]
		-replicas[[Key[i],Key["state"],Key["energy"]]]]>=0.00001
		,Print[" Found mismatch for replica index i = ",i]]
	,{i,numRep}];
	(*End diagnostic*)
	*)

	swap[numsweeps];

	(*(*For diagnostic*)
	Print[" Energy check after swap "];
	Table[
		If[Abs[hamiltonian[replicas[[Key[i],Key["state"],Key["graph"]]],replicas[[Key[i],Key["externalField"]]]]
		-replicas[[Key[i],Key["state"],Key["energy"]]]]>=0.00001
		,Print[" Found mismatch for replica index i = ",i," state ",replicas[[Key[i]]]]]
	,{i,numRep}];
	(*End diagnostic*)
	*)


	If[Mod[numsweeps,1]==0,
	Table[
		Sow[{numsweeps,externalFieldTable[[i]],replicas[[Key[i],Key["state"],Key["energy"]]]
			,Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]]
			,replicas[[Key[i],Key["state"],Key["energy"]]]
				-(externalFieldTable[[i]] .
					Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]])
			,Through[obs[replicas[[Key[i],Key["state"],Key["graph"]]]]]},i]
	,{i,numRep}];
	];

	numsweeps++;

]][[2]];


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[externalFieldTable]]]}]|>;


Do[If[histories[[Key[i],Key["swapAccept"]]]<Length[histories[[Key[i]
	,Key["history"]]]],histories[[Key[i],Key["history"]]]
		=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]])]
,{i,numRep}];
(*remove the -1's in the initiation *)

(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each value of the externalFieldTable, i.e., (-beta*F). *)
(**********************************************)
*)

mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[chart[[All,All,4]],bt];

result={<|Table[replicas[[Key[i],Key["externalField"]]]->groundStates[[Key[i]]],{i,Keys[replicas]}]|>
	,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"beta"->bt,"externalField"->i,"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>
	,replicas,histories};

Remove[groundStates,maxGStateCount,replicas,replicaKeysOrderedExternalField,numRep,histories,swap,Tempoutput,EnergyOrMag ,measurements,numsweeps,stopnum, candminE,printCase,repNumSweeps,chart,mBetaF];

result

];


(*Overload with multiple external fields no delH*)
ECGrav`GraphMultiHistogram[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	externalFieldTable_List/;MatrixQ[externalFieldTable],
	conjugateObs_/;MatchQ[conjugateObs,{__Function}],obs_/;MatchQ[obs,{__Function}],NN_Integer,
	UnlabeledVerticesYes_Integer]:=
(*************************************)
(***  Last updated on: 3/04/2026  ***)
(*************************************)

(*An overload which implements the Multiple Histogram Method for the graph models to 
get a smooth plot of the quantities obs as a function of multiple external parameters at a 
fixed inverse temperature. It applies replica swaps during the measurement stage.,

(*Notes: ,

*)

Inputs are:,
1. seedGraph - adjacency matrix of the input seed graph,
2. bt - Real = a fixed value of the inverse temperature,
3. hamiltonian[hparams] =  the hamiltonian function with its parameters hparams that 
	assigns graphs energy,
4. externalFieldTable - List = a list of list of external field values {{J^1_1,J^2_1},{J^1_2,J^2_2}...},
5. conjugateObs_ = a list of  functions O which are conjugate to the external fields 
	J^i i.e., there is a term in the microscopic hamiltonian of the form Sum_i(J^i*O^i(s)),
6. obs = a list of various observables,
7. NN = number of independent sweeps (so that actual number of sweeps is correlation 
	time times NN),
8. UnlabeledVerticesYes = 0 means no selection probability to make the graphs 
	unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with three entries:,
1. an association of the minimum energy found from the run and the states found 
	having that energy.,
2. an association with external field table values as keys and values of 
	negative*beta*free energy and values for the energies and observables at each 
	external field value,
3. the replicas at the last step,
4. an association of the histories of each replica *)

Module[{result,groundStates,maxGStateCount,replicas,replicaKeysOrderedExternalField,numRep,histories,
	swap,Tempoutput,measurements,numsweeps,stopnum, candminE,printCase,repNumSweeps,chart,mBetaF},


groundStates=<||>;
histories=<||>;
numRep=Length[externalFieldTable];


maxGStateCount=500;(*Maximum count for lowest energy states to be saved.*)

(*In general, sqrt[number of sites] replicas are needed for parallel tempering to be effective.*)

PrintTemporary[" Starting multihistogram with replical swap with delH and with beta "
	,bt," external field parameters ",externalFieldTable];


(*
(********************)
(*   Equilibriate   *)
(********************)
*)


Tempoutput=Association[
	ParallelTable[
		i->ECGrav`GraphEquilibriate[seedGraph,bt
					,Apply[hamiltonian,externalFieldTable[[i]]]
					,UnlabeledVerticesYes]
	,{i,Length[externalFieldTable]},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

(* Extract the minimum energy and corresponding states *)

groundStates=Tempoutput[[All,1]];


(*
(***************************************)
(*  Compute correlation times   *)
(***************************************)
*)

Tempoutput=Association[
	ParallelTable[
		i->ECGrav`GraphComputeCorrelationTime[
			replicas[[Key[i],Key["state"],Key["graph"]]],bt
			,Apply[hamiltonian,externalFieldTable[[i]]]
			,replicas[[Key[i],Key["eqlT"]]]
			,groundStates[[Key[i],Key["minEnergy"]]],conjugateObs,UnlabeledVerticesYes]
		,{i,numRep},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}
	]
];


(*
(*Update replicas*)
*)

Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];


(*Update ground states*)

Do[

	(* Extract the minimum energy and corresponding states *)
	candminE=Tempoutput[[Key[i],1,"minEnergy"]];


	If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
		groundStates[[Key[i],Key["minEnergy"]]]=candminE;
		groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
		If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
			groundStates[[Key[i],Key["minEstates"]]]=
				Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
		];
	];

,{i,numRep}];

replicaKeysOrderedExternalField=Keys[Sort[replicas[[All,"externalField"]]]];

histories=<|Table[i-><|"externalField"->replicas[[Key[i],Key["externalField"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["externalField"]]],{i,numRep}];

(*(*********************)
(*  The swap function  *)
(*********************)*)

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-beta(delta )(delta E)). *)
Module[{thisReplInd,nextReplInd,betaDelHDelConj,accept,tempThisReplicaState},

Do[
	thisReplInd=replicaKeysOrderedExternalField[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedExternalField[[Mod[offset+repl,numRep]+1]];


	betaDelHDelConj=
		bt*((replicas[[Key[thisReplInd],Key["externalField"]]]
				-replicas[[Key[nextReplInd],Key["externalField"]]])
			. (Through[conjugateObs[replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]]]
				-Through[conjugateObs[replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]]])
			);


	histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;
	histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;

	accept=0;
	If[betaDelHDelConj>0,accept = 1,
		If[RandomReal[]<Exp[betaDelHDelConj],accept =1];
	];


	If[accept==1, (*Do swap of replicas*)

		(*Update states: store the state of the high temp replica temporarily*)

		tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
		replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
		replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

		(*Update the energies to reflect the new external field values*)
	
		replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]
			=Apply[hamiltonian
				,Join[{replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]}
					,replicas[[Key[thisReplInd],Key["externalField"]]]]];
		replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]
			=Apply[hamiltonian
				,Join[{replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]}
					,replicas[[Key[nextReplInd],Key["externalField"]]]]];

		(*Update histories*)

		histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
		histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

		histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
		histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

		(*Increase the number of swap accepts by one for both replicas*)
		histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
		histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;


	];


,{repl,1,numRep-1,2}];

];


(*
(****************************************)
(* Take NN measurements at the           *
* equilibriated temperatures to compute *
* specific heats *)
(****************************************)
*)

numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];

measurements=Reap[
While[numsweeps<=NN,

	If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

	Tempoutput=Association[
		ParallelTable[
			(* Each replica will be swept corrT times *)
			<|i->ECGrav`GraphSweepReplica[replicas[[Key[i],Key["state"],Key["graph"]]]
				,bt,Apply[hamiltonian,externalFieldTable[[i]]]
				,replicas[[Key[i],Key["corrT"]]],groundStates[[Key[i],Key["minEnergy"]]]
				,UnlabeledVerticesYes]|>,

			{i,numRep},DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}
			]
	];


	(*Update replicas*)
	Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i],2]],{i,numRep}];

	(*Update ground states*)

	Do[

		(* Extract the minimum energy and corresponding states *)
		candminE=Tempoutput[[Key[i],1,"minEnergy"]];


		If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
			groundStates[[Key[i],Key["minEnergy"]]]=candminE;
			groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
			If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
				groundStates[[Key[i],Key["minEstates"]]]=
					Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
			];
		];

	,{i,numRep}];

	(*(*For diagnostic*)
	Print[""];
	Print[""];
	Print[""];
	Print[""];Print[""];Print[""];
	Print[" Energy check before swap "];
	Table[
		If[Abs[hamiltonian[replicas[[Key[i],Key["state"],Key["graph"]]],replicas[[Key[i],Key["externalField"]]]]
		-replicas[[Key[i],Key["state"],Key["energy"]]]]>=0.00001
		,Print[" Found mismatch for replica index i = ",i]]
	,{i,numRep}];
	(*End diagnostic*)
	*)

	swap[numsweeps];

	(*(*For diagnostic*)
	Print[" Energy check after swap "];
	Table[
		If[Abs[hamiltonian[replicas[[Key[i],Key["state"],Key["graph"]]],replicas[[Key[i],Key["externalField"]]]]
		-replicas[[Key[i],Key["state"],Key["energy"]]]]>=0.00001
		,Print[" Found mismatch for replica index i = ",i," state ",replicas[[Key[i]]]]]
	,{i,numRep}];
	(*End diagnostic*)
	*)


	If[Mod[numsweeps,1]==0,
	Table[
		Sow[{numsweeps,externalFieldTable[[i]],replicas[[Key[i],Key["state"],Key["energy"]]]
			,Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]]
			,replicas[[Key[i],Key["state"],Key["energy"]]]
				-(externalFieldTable[[i]] .
					Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]])
			,Through[obs[replicas[[Key[i],Key["state"],Key["graph"]]]]]},i]
	,{i,numRep}];
	];

	numsweeps++;

]][[2]];


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[externalFieldTable]]]}]|>;


Do[If[histories[[Key[i],Key["swapAccept"]]]<Length[histories[[Key[i]
	,Key["history"]]]],histories[[Key[i],Key["history"]]]
		=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]])]
,{i,numRep}];
(*remove the -1's in the initiation *)
(*
(*********************************************)
(*  Compute -beta times the free energies at   
* each value of the externalFieldTable, i.e., (-beta*F). *)
(**********************************************)
*)


mBetaF=ECGrav`ComputeMinusBetaTimesFreeEnergy[chart[[All,All,4]],bt];

result={<|Table[replicas[[Key[i],Key["externalField"]]]->groundStates[[Key[i]]],{i,Keys[replicas]}]|>
	,<|Table[i-><|"minusBetaF"->mBetaF[[Key[i]]],"beta"->bt,"externalField"->i,"data"->chart[[Key[i]]]|>,{i,Keys[chart]}]|>
	,replicas,histories};

Remove[groundStates,maxGStateCount,replicas,replicaKeysOrderedExternalField,numRep,histories,swap,Tempoutput,EnergyOrMag ,measurements,numsweeps,stopnum, candminE,printCase,repNumSweeps,chart,mBetaF];

result

];


(* ::Item::Closed:: *)
(*GraphMultiHistogram Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphMultiHistogram[args___]:=(Message[ECGrav`GraphMultiHistogram::argerr, args];
$Failed);


(* ::Subsection::Closed:: *)
(*GraphCEITempSchedule*)


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Primary*)


ECGrav`GraphCEITempSchedule[seedGraph_List,betaLow_Real,betaHigh_Real,
	hamiltonian_[hparams___],
	delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
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


hist=ECGrav`GraphMultiHistogram[seedGraph,betaLow,betaHigh,hamiltonian[hparams],
		delH[delHparams],obs,NN,UnlabeledVerticesYes];


minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)


incrementVal=(1.0/betaLow-1.0/betaHigh)/(20*vCount);


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


entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];


tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];


result={hist[[1]],tempSchedule,hist[[2]],hist[[3]]};

Remove[vCount,hist,minusbetaFTable,chart,betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule];

result
];


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Overload 1 no delH*)


ECGrav`GraphCEITempSchedule[seedGraph_List,betaLow_Real,betaHigh_Real,hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
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


hist=ECGrav`GraphMultiHistogram[seedGraph,betaLow,betaHigh,hamiltonian[hparams],obs,NN,UnlabeledVerticesYes];


minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)


incrementVal=(1.0/betaLow-1.0/betaHigh)/(20*vCount);


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


entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];


tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];


result={hist[[1]],tempSchedule,hist[[2]],hist[[3]]};

Remove[vCount,hist,minusbetaFTable,chart,betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule];

result
];


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Overload 2 betaTable*)


ECGrav`GraphCEITempSchedule[seedGraph_List,btTable_List,hamiltonian_[hparams___],
	delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
Module[{result,vCount=Length[seedGraph],hist,minusbetaFTable,chart,betas,incrementVal,
        CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule},


hist=ECGrav`GraphMultiHistogram[seedGraph,btTable,hamiltonian[hparams],
		delH[delHparams],obs,NN,UnlabeledVerticesYes];


minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)


incrementVal=(1.0/Min[btTable]-1.0/Max[btTable])/(20*vCount);


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


entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];


tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];


result={hist[[1]],tempSchedule,hist[[2]],hist[[3]]};

Remove[vCount,hist,minusbetaFTable,chart,betas,incrementVal,CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule];

result
];


(* ::Item::Closed:: *)
(*GraphCEITempSchedule Overload 3 betaTable and no delH*)


ECGrav`GraphCEITempSchedule[seedGraph_List,btTable_List,hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer]:=
Module[{result,vCount=Length[seedGraph],hist,minusbetaFTable,chart,betas,incrementVal,
        CvOverTtab,entropyRange,delS,entropyTable,entropyVals,tempSchedule},


hist=ECGrav`GraphMultiHistogram[seedGraph,btTable,hamiltonian[hparams],obs,NN,UnlabeledVerticesYes];


minusbetaFTable=hist[[2,All,"minusBetaF"]];
chart=hist[[2,All,"data",All,3]];

betas=Keys[minusbetaFTable];

(*Print[" After multihistogram, betas ", betas];
Print[" After multihistogram, minusbetaFTable ",minusbetaFTable];
Print[" After multihistogram, measurements ", chart];*)


incrementVal=(1.0/Min[btTable]-1.0/Max[btTable])/(20*vCount);


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


entropyTable=Table[{CvOverTtab[[k,1]],Sum[CvOverTtab[[i,2]]*incrementVal,{i,1,k}]},{k,1,Length[CvOverTtab]}];(*integral of Cv/T dT gives entropy*)

(*Print["entropyTable ",entropyTable];*)
Print["entropyTable plot "];
Print[ListLinePlot[entropyTable,PlotRange->All,PlotLabel->"entropy vs temp plot"]];

(*entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,1,vCount}]]];*)

entropyVals=Flatten[Nearest[entropyTable[[All,2]],Table[i*delS,{i,0,vCount}]]];


tempSchedule=DeleteDuplicates[Table[Flatten[Select[entropyTable,#[[2]]==i&]],{i,entropyVals}]];


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
(*GraphCTLSchedule*)


(* ::Item::Closed:: *)
(*GraphCTLSchedule for beta with input beta list*)


ECGrav`GraphCTLSchedule[seedGraph_List,btTable_List,hamiltonian_[hparams___],
	delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer,
	numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 02/16/2026  ***)
(*************************************)
(*Notes: , *)
(*A program to find an inverse temperature schedule for parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach.,

Inputs are:,
1. seedGraph_List - adjacency matrix of the input seed graph,
2. btTable_List = a list of beta values for an initial run,
3. hamiltonian_[hparams___] =  the hamiltonian with its parameters that assigns graphs energy,
4. delH_[delHparams___] = function that gives the change in energy when a single edge is flipped,
5. obs = the observable quantity(ies) to be measured,
6. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
7. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,
8. numReplicas_Integer = the number of replicas that will be used for parallel tempering. 

Outputs a list with the following elements:,
1. the inverse temperature schedule,
2. an association of the output from the MultiHistogram run)

*)
Module[{result,vCount=Length[seedGraph],hist,minusbetaFAssn,betaMin,betaMax,
	thermodynamicSpeed,chartEnergy,incrementValue,thermodynamicSpeedTable,
	interpolatedThermodynamicSpeed,totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,betaSchedule,b},


hist=ECGrav`GraphMultiHistogram[seedGraph,btTable,hamiltonian[hparams],
	delH[delHparams],obs,NN,UnlabeledVerticesYes];


minusbetaFAssn=hist[[2,All,"minusBetaF"]];
chartEnergy=hist[[2,All,"data",All,3]];

betaMin=Min[Keys[minusbetaFAssn]];
betaMax=Max[Keys[minusbetaFAssn]];


thermodynamicSpeed[b_]:=Sqrt[ECGrav`ExtrapolatedExpectationValue[b,minusbetaFAssn,chartEnergy,chartEnergy^2]
		-(ECGrav`ExtrapolatedExpectationValue[b,minusbetaFAssn,chartEnergy,chartEnergy])^2];


incrementValue=(betaMax-betaMin)/(10.0*vCount);


thermodynamicSpeedTable=ParallelTable[{b,thermodynamicSpeed[b]},
	{b,betaMin,betaMax*(1.1),incrementValue},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}];

interpolatedThermodynamicSpeed=Interpolation[thermodynamicSpeedTable];

Print[Plot[interpolatedThermodynamicSpeed[b],{b,betaMin,betaMax},
	PlotRange->All,PlotLabel->"interpolated thermodynamic speed"]];

thermodynamicDistanceTable=Table[{b,NIntegrate[interpolatedThermodynamicSpeed[x],
		{x,betaMin,b}]},{b,betaMin,betaMax*(1.1),incrementValue}];


interpolatedThermodynamicLength=Interpolation[thermodynamicDistanceTable];

totalThermodynamicLength = NIntegrate[interpolatedThermodynamicSpeed[x],
	{x,betaMin,betaMax}];

Print[Plot[interpolatedThermodynamicLength[b],{b,betaMin,betaMax},
	PlotRange->All,PlotLabel->" interpolated thermodynamic distance "]];

(*Print[" total thermodynamic length ",totalThermodynamicLength];*)

targetSteps=Table[i*(totalThermodynamicLength/(numReplicas-1)),{i,0,numReplicas-1}];

betaSchedule=Table[b/. FindRoot[interpolatedThermodynamicLength[b]==targetSteps[[i]]
	,{b,betaMin+(betaMax-betaMin)*(i-1)/(numReplicas-1),betaMin,betaMax*(1.1)}],{i,1,numReplicas}];


result={betaSchedule,hist};

Remove[vCount,hist,minusbetaFAssn,betaMin,betaMax,thermodynamicSpeed,chartEnergy,
	incrementValue,thermodynamicSpeedTable,interpolatedThermodynamicSpeed,
	totalThermodynamicLength,thermodynamicDistanceTable,interpolatedThermodynamicLength,
	targetSteps,betaSchedule,b];

result
];

ECGrav`GraphCTLSchedule[seedGraph_List,btTable_List,hamiltonian_[hparams___],
	obs_,NN_Integer,UnlabeledVerticesYes_Integer,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 02/16/2026  ***)
(*************************************)
(*Notes: , *)
(*A program to find an inverse temperature schedule for parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach.,

Inputs are:,
1. seedGraph_List - adjacency matrix of the input seed graph,
2. btTable_List = a list of beta values for an initial run,
3. hamiltonian_[hparams___] =  the hamiltonian with its parameters that assigns graphs energy,
4. obs = the observable quantity(ies) to be measured,
5. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
6. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,
7. numReplicas_Integer = the number of replicas that will be used for parallel tempering.,

Outputs a list with the following elements:,
1. the inverse temperature schedule,
2. an association of the output from the MultiHistogram run)
*)

Module[{result,vCount=Length[seedGraph],hist,minusbetaFAssn,betaMin,betaMax,
	thermodynamicSpeed,chartEnergy,incrementValue,thermodynamicSpeedTable,
	interpolatedThermodynamicSpeed,totalThermodynamicLength,
	thermodynamicDistanceTable,interpolatedThermodynamicLength,targetSteps,
	betaSchedule,b},


hist=ECGrav`GraphMultiHistogram[seedGraph,btTable,hamiltonian[hparams],obs,NN,
	UnlabeledVerticesYes];

minusbetaFAssn=hist[[2,All,"minusBetaF"]];
chartEnergy=hist[[2,All,"data",All,3]];

betaMin=Min[Keys[minusbetaFAssn]];
betaMax=Max[Keys[minusbetaFAssn]];


thermodynamicSpeed[b_]:=Sqrt[ECGrav`ExtrapolatedExpectationValue[b,minusbetaFAssn,chartEnergy,chartEnergy^2]
	-(ECGrav`ExtrapolatedExpectationValue[b,minusbetaFAssn,chartEnergy,chartEnergy])^2];


incrementValue=(betaMax-betaMin)/(10.0*vCount);

thermodynamicSpeedTable=ParallelTable[{b,thermodynamicSpeed[b]}
	,{b,betaMin,betaMax*(1.1),incrementValue},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}];

interpolatedThermodynamicSpeed=Interpolation[thermodynamicSpeedTable];

Print[Plot[interpolatedThermodynamicSpeed[b],{b,betaMin,betaMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic speed"]];

thermodynamicDistanceTable=Table[{b,NIntegrate[interpolatedThermodynamicSpeed[x],
	{x,betaMin,b}]},{b,betaMin,betaMax*(1.1),incrementValue}];


interpolatedThermodynamicLength=Interpolation[thermodynamicDistanceTable];

totalThermodynamicLength = NIntegrate[interpolatedThermodynamicSpeed[x],
	{x,betaMin,betaMax}];

Print[Plot[interpolatedThermodynamicLength[b],{b,betaMin,betaMax},
	PlotRange->All,PlotLabel->" interpolated thermodynamic distance "]];


targetSteps=Table[i*(totalThermodynamicLength/(numReplicas-1)),{i,0,numReplicas-1}];


betaSchedule=Table[b/. FindRoot[interpolatedThermodynamicLength[b]==targetSteps[[i]],
	{b,betaMin+(betaMax-betaMin)*(i-1)/(numReplicas-1),betaMin,betaMax*(1.1)}],
	{i,1,numReplicas}];


result={betaSchedule,hist};


Remove[vCount,hist,minusbetaFAssn,betaMin,betaMax,thermodynamicSpeed,chartEnergy,
	incrementValue,thermodynamicSpeedTable,interpolatedThermodynamicSpeed,
	totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,betaSchedule,b];

result
];


(* ::Item::Closed:: *)
(*GraphCTLSchedule for  beta with input from a previous run*)


ECGrav`GraphCTLSchedule[betaMin_Real,betaMax_Real,minusbetaF_Association,
	energyValues_Association,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 02/16/2026  ***)
(*************************************)
(*Notes: , *)
(*An overload  to find an inverse temperature schedule for parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach. It takes as an input the return from a previous run of GraphMultihistogram ,

Inputs are:,
1. betaMin_Real = the minimum beta value for the schedule to be built,
2. betaMax_Real = the maximum beta value for the schedule to be built,
3. minusbetaF_Association = an association of the inverse temperature and corresponding values of -beta*freeEnergy from a previous run of GraphMultihistograph,
4. energyValues_Association = an association of the inverse temperature and corresponding table of energy measurements from a previous run of GraphMultihistograph,
5. numReplicas_Integer = the number of replicas that will be used for parallel tempering. 

Outputs a single list of the beta schedule.*)

Module[{result,thermodynamicSpeed,incrementValue,thermodynamicSpeedTable,
	interpolatedThermodynamicSpeed,totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,betaSchedule,b},


thermodynamicSpeed[b_]:=Sqrt[ECGrav`ExtrapolatedExpectationValue[b,minusbetaF,energyValues,energyValues^2]
	-(ECGrav`ExtrapolatedExpectationValue[b,minusbetaF,energyValues,energyValues])^2];

incrementValue=(betaMax-betaMin)/(100.0);

thermodynamicSpeedTable=ParallelTable[{b,thermodynamicSpeed[b]}
	,{b,betaMin,betaMax*(1.1),incrementValue},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}];

interpolatedThermodynamicSpeed=Interpolation[thermodynamicSpeedTable];

(*Print["interpolatedThermodynamicSpeed plot "];*)
Print[Plot[interpolatedThermodynamicSpeed[b],{b,betaMin,betaMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic speed"]];

thermodynamicDistanceTable=Table[{b,NIntegrate[interpolatedThermodynamicSpeed[x],
	{x,betaMin,b}]},{b,betaMin,betaMax*(1.1),incrementValue}];

interpolatedThermodynamicLength=Interpolation[thermodynamicDistanceTable];

totalThermodynamicLength = NIntegrate[interpolatedThermodynamicSpeed[x],
	{x,betaMin,betaMax}];

Print[Plot[interpolatedThermodynamicLength[b],{b,betaMin,betaMax},PlotRange->All,
	PlotLabel->" interpolated thermodynamic distance "]];

targetSteps=Table[i*(totalThermodynamicLength/(numReplicas-1)),{i,0,numReplicas-1}];


betaSchedule=Table[b/. FindRoot[interpolatedThermodynamicLength[b]==targetSteps[[i]],{b,betaMin+(betaMax-betaMin)*(i-1)/(numReplicas-1),betaMin,betaMax*(1.1)}],{i,1,numReplicas}];


result=betaSchedule;


Remove[incrementValue,thermodynamicSpeedTable,interpolatedThermodynamicSpeed,
	totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,betaSchedule,b];

result
];


(* ::Item::Closed:: *)
(*GraphCTLSchedule for a single external field at a fixed beta with input external field  list*)


ECGrav`GraphCTLSchedule[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	delH_[delHparams___],
	externalFieldTable_List/;(MatrixQ[externalFieldTable]&&Length[externalFieldTable[[1]]]==1),
	conjugateObs_/;MatchQ[conjugateObs,{_Function}],
	obs_/;MatchQ[obs,{__Function}],NN_Integer,
	UnlabeledVerticesYes_Integer,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 03/05/2026  ***)
(*************************************)
(*Notes: , *)
(*A program to find a schedule of external field values for use in parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach.,

Inputs are:,
1. seedGraph_List - adjacency matrix of the input seed graph,
2. bt_Real = a fixed value of the inverse temperature,
3. hamiltonian_[hparams___] =  the hamiltonian with its parameters that assigns graphs energy,
4. delH_[delHparams___] = function that gives the change in energy when a single edge is flipped,
5. externalFieldTable_List = a list of real number values of the external field for use in the initial run (e.g. a table of H values for the Ising model),
6. conjugateObs_ = a function, the conjugate observable to the external field (e.g. the magnetization formula in the Ising model),
7. obs = the observable quantity(ies) to be measured,
8. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
9. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,
10. numReplicas_Integer = the number of replicas that will be used for parallel tempering. 

Outputs a list with the following elements:,
1. the external field schedule (a list containing the external field values, 
	edge list for swaps, and an association matching external field values to edge labels),
2. an association of the output from the MultiHistogram run)
*)

Module[{result,vCount=Length[seedGraph],hist,minusbetaFAssn,hMin,hMax,
	thermodynamicSpeed,chartConjugateField,incrementValue,thermodynamicSpeedTable,
	interpolatedThermodynamicSpeed,totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,externalFieldSchedule,h,
	neighbors,edges,replicaLabels},

PrintTemporary[" Starting GraphCTLSchedule for external field at a fixed beta ",bt
	," with input external field schedule, ",externalFieldTable," numReplicas "
	,numReplicas];
	
hist=ECGrav`GraphMultiHistogram[seedGraph,bt,hamiltonian[hparams],delH[delHparams],
	externalFieldTable,conjugateObs,obs,NN,UnlabeledVerticesYes];

minusbetaFAssn=hist[[2,All,"minusBetaF"]];
chartConjugateField=hist[[2,All,"data",All,4]];


hMin=Min[Keys[minusbetaFAssn]];
hMax=Max[Keys[minusbetaFAssn]];


thermodynamicSpeed[h_]:=
	Sqrt[ECGrav`ExtrapolatedExpectationValue[bt,{h},minusbetaFAssn,
							chartConjugateField,chartConjugateField[[All,All,1]]^2]
		-(ECGrav`ExtrapolatedExpectationValue[bt,{h},minusbetaFAssn,chartConjugateField,
					chartConjugateField[[All,All,1]]])^2];


incrementValue=(hMax-hMin)/(10.0*vCount);

thermodynamicSpeedTable=ParallelTable[{h,thermodynamicSpeed[h]}
	,{h,hMin-0.2*Abs[hMax-hMin],hMax+0.2*Abs[hMax-hMin],incrementValue},
	DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}];

interpolatedThermodynamicSpeed=Interpolation[thermodynamicSpeedTable];

Print[Plot[interpolatedThermodynamicSpeed[h],{h,hMin,hMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic speed"]];

thermodynamicDistanceTable=Table[{h,NIntegrate[interpolatedThermodynamicSpeed[x],
	{x,hMin,h}]},{h,hMin,hMax,incrementValue}];


totalThermodynamicLength = NIntegrate[interpolatedThermodynamicSpeed[x],{x,hMin,hMax}];

interpolatedThermodynamicLength=Interpolation[thermodynamicDistanceTable];

Print[Plot[interpolatedThermodynamicLength[h],{h,hMin,hMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic distance"]];


targetSteps=Table[i*(totalThermodynamicLength/(numReplicas-1)),{i,0,numReplicas-1}];

externalFieldSchedule=Table[h/. FindRoot[interpolatedThermodynamicLength[h]==targetSteps[[i]]
	,{h,hMin+(hMax-hMin)*(i-1)/(numReplicas-1),hMin,hMax}],{i,1,numReplicas}];

externalFieldSchedule={#}&/@externalFieldSchedule;

neighbors=Join[{{2,3}},Table[{i-1,i+1},{i,2,numReplicas-1}],{{numReplicas-2,numReplicas-1}}];

edges=DeleteDuplicates@(
		Sort/@Flatten[Table[UndirectedEdge[i,j],{i,numReplicas},{j,neighbors[[i]]}],1]);


replicaLabels=<|Table[i->externalFieldSchedule[[i]],{i,numReplicas}]|>;

result={{externalFieldSchedule,edges,replicaLabels},hist};

Remove[vCount,hist,minusbetaFAssn,hMin,hMax,thermodynamicSpeed,chartConjugateField,
	totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,externalFieldSchedule,h,
	neighbors,edges,replicaLabels];

result
];


ECGrav`GraphCTLSchedule[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	externalFieldTable_List/;(MatrixQ[externalFieldTable]&&Length[externalFieldTable[[1]]]==1),
	conjugateObs_/;MatchQ[conjugateObs,{_Function}],
	obs_/;MatchQ[obs,{__Function}],NN_Integer,
	UnlabeledVerticesYes_Integer,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 03/05/2026  ***)
(*************************************)
(*Notes: , *)
(*A program to find a schedule of external field values for use in parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach.,

Inputs are:,
1. seedGraph_List - adjacency matrix of the input seed graph,
2. bt_Real = a fixed value of the inverse temperature,
3. hamiltonian_[hparams___] =  the hamiltonian with its parameters that assigns graphs energy,
4. externalFieldTable_List = a list of real number values of the external field for use in the initial run (e.g. a table of H values for the Ising model),
5. conjugateObs_ = a function, the conjugate observable to the external field (e.g. the magnetization formula in the Ising model),
6. obs = the observable quantity(ies) to be measured,
7. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
8. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,
9. numReplicas_Integer = the number of replicas that will be used for parallel tempering. 

Outputs a list with the following elements:,
1. the external field schedule (a list containing the external field values, 
	edge list for swaps, and an association matching external field values to edge labels),
2. an association of the output from the MultiHistogram run)
*)

Module[{result,vCount=Length[seedGraph],hist,minusbetaFAssn,hMin,hMax,
	thermodynamicSpeed,chartConjugateField,incrementValue,thermodynamicSpeedTable,
	interpolatedThermodynamicSpeed,totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,externalFieldSchedule,h,
	neighbors,edges,replicaLabels},


hist=ECGrav`GraphMultiHistogram[seedGraph,bt,hamiltonian[hparams],
	externalFieldTable,conjugateObs,obs,NN,UnlabeledVerticesYes];

minusbetaFAssn=hist[[2,All,"minusBetaF"]];
chartConjugateField=hist[[2,All,"data",All,4]];


hMin=Min[Keys[minusbetaFAssn]];
hMax=Max[Keys[minusbetaFAssn]];


thermodynamicSpeed[h_]:=Sqrt[ECGrav`ExtrapolatedExpectationValue[bt,{h},minusbetaFAssn,chartConjugateField,chartConjugateField[[All,All,1]]^2]
	-(ECGrav`ExtrapolatedExpectationValue[bt,{h},minusbetaFAssn,chartConjugateField,chartConjugateField[[All,All,1]]])^2];


incrementValue=(hMax-hMin)/(10.0*vCount);

thermodynamicSpeedTable=ParallelTable[{h,thermodynamicSpeed[h]}
	,{h,hMin-0.2*Abs[hMax-hMin],hMax+0.2*Abs[hMax-hMin],incrementValue},
	DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}];

interpolatedThermodynamicSpeed=Interpolation[thermodynamicSpeedTable];

Print[Plot[interpolatedThermodynamicSpeed[h],{h,hMin,hMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic speed"]];

thermodynamicDistanceTable=Table[{h,NIntegrate[interpolatedThermodynamicSpeed[x],
	{x,hMin,h}]},{h,hMin,hMax,incrementValue}];


totalThermodynamicLength = NIntegrate[interpolatedThermodynamicSpeed[x],{x,hMin,hMax}];

interpolatedThermodynamicLength=Interpolation[thermodynamicDistanceTable];

Print[Plot[interpolatedThermodynamicLength[h],{h,hMin,hMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic distance"]];


targetSteps=Table[i*(totalThermodynamicLength/(numReplicas-1)),{i,0,numReplicas-1}];

externalFieldSchedule=Table[h/. FindRoot[interpolatedThermodynamicLength[h]==targetSteps[[i]]
	,{h,hMin+(hMax-hMin)*(i-1)/(numReplicas-1),hMin,hMax}],{i,1,numReplicas}];

externalFieldSchedule={#}&/@externalFieldSchedule;

neighbors=Join[{{2,3}},Table[{i-1,i+1},{i,2,numReplicas-1}],{{numReplicas-2,numReplicas-1}}];

edges=DeleteDuplicates@(
		Sort/@Flatten[Table[UndirectedEdge[i,j],{i,numReplicas},{j,neighbors[[i]]}],1]);


replicaLabels=<|Table[i->externalFieldSchedule[[i]],{i,numReplicas}]|>;

result={{externalFieldSchedule,edges,replicaLabels},hist};

Remove[vCount,hist,minusbetaFAssn,hMin,hMax,thermodynamicSpeed,chartConjugateField,
	totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,externalFieldSchedule,h,
	neighbors,edges,replicaLabels];

result
];


(* ::Item::Closed:: *)
(*GraphCTLSchedule for a single external field at a fixed beta with input from a previous run*)


ECGrav`GraphCTLSchedule[bt_Real,hMin_Real,hMax_Real,minusbetaF_Association,
	conjugateFieldMeasurements_Association,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 03/05/2026  ***)
(*************************************)
(*Notes: , *)
(*An overload  to find a schedule of external field values for use in parallel 
tempering that keeps acceptance ratios constant based on the constant thermodynamic 
length approach. It takes as an input the return from a previous run of 
GraphMultihistogram,

Inputs are:,
1. bt_Real = a fixed value of the inverse temperature,
2. hMin_Real = the minimum value of the external field for the schedule to be built,
3. hMax_Real = the maximum value of the external field for the schedule to be built,
4. minusbetaF_Association = an association of external field values and corresponding values of -beta*freeEnergy from a previous run of GraphMultihistograph,
5. conjugateFieldMeasurements_Association = an association of external field values and the corresponding table of the conjugate observable measurements from a previous run of GraphMultihistograph,
6. numReplicas_Integer = the number of replicas that will be used for parallel tempering.,

Outputs a single list of the external field schedule.*)

Module[{result,thermodynamicSpeed,incrementValue,thermodynamicSpeedTable,
interpolatedThermodynamicSpeed,totalThermodynamicLength,thermodynamicDistanceTable,
interpolatedThermodynamicLength,targetSteps,externalFieldSchedule,h,
neighbors,edges,replicaLabels},

PrintTemporary[" Starting GraphCTLSchedule for external field at a fixed beta ",bt
	," with input from a previous run of GraphMultiHistogram, numReplicas "
	,numReplicas];	

thermodynamicSpeed[h_]:=Sqrt[ECGrav`ExtrapolatedExpectationValue[bt,{h},minusbetaF,
	conjugateFieldMeasurements,conjugateFieldMeasurements[[All,All,1]]^2]
	-(ECGrav`ExtrapolatedExpectationValue[bt,{h},minusbetaF,conjugateFieldMeasurements,
	conjugateFieldMeasurements[[All,All,1]]])^2];

incrementValue=(hMax-hMin)/(100.0);

thermodynamicSpeedTable=ParallelTable[{h,thermodynamicSpeed[h]}
	,{h,hMin-0.2*Abs[hMax-hMin],hMax+0.2*Abs[hMax-hMin],incrementValue}
	,DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}];

interpolatedThermodynamicSpeed=Interpolation[thermodynamicSpeedTable];

Print[Plot[interpolatedThermodynamicSpeed[h],{h,hMin,hMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic speed"]];

thermodynamicDistanceTable=Table[{h,NIntegrate[interpolatedThermodynamicSpeed[x]
	,{x,hMin,h}]},{h,hMin,hMax,incrementValue}];

interpolatedThermodynamicLength=Interpolation[thermodynamicDistanceTable];

totalThermodynamicLength = NIntegrate[interpolatedThermodynamicSpeed[x],{x,hMin,hMax}];

Print[Plot[interpolatedThermodynamicLength[h],{h,hMin,hMax},PlotRange->All,
	PlotLabel->"interpolated thermodynamic distance"]];


targetSteps=Table[i*(totalThermodynamicLength/(numReplicas-1)),{i,0,numReplicas-1}];

externalFieldSchedule=Table[h/. FindRoot[interpolatedThermodynamicLength[h]==targetSteps[[i]]
	,{h,hMin+(hMax-hMin)*(i-1)/(numReplicas-1),hMin,hMax}],{i,1,numReplicas}];

externalFieldSchedule={#}&/@externalFieldSchedule;

neighbors=Join[{{2,3}},Table[{i-1,i+1},{i,2,numReplicas-1}],{{numReplicas-2,numReplicas-1}}];

edges=DeleteDuplicates@(
	Sort/@Flatten[Table[UndirectedEdge[i,j],{i,numReplicas},{j,neighbors[[i]]}],1]);


replicaLabels=<|Table[i->externalFieldSchedule[[i]],{i,numReplicas}]|>;

result={externalFieldSchedule,edges,replicaLabels};

Remove[thermodynamicSpeed,totalThermodynamicLength,thermodynamicDistanceTable,
	interpolatedThermodynamicLength,targetSteps,externalFieldSchedule,h,
	neighbors,edges,replicaLabels];

result
];


(* ::Item::Closed:: *)
(*GraphCTLSchedule for multiple external fields at a fixed beta with input external field  matrix*)


(*Primary pattern*)
ECGrav`GraphCTLSchedule[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	delH_[delHparams___],
	externalFieldTable_List/;(MatrixQ[externalFieldTable]&&Length[externalFieldTable[[1]]]>1),
	conjugateObs_/;MatchQ[conjugateObs,{__Function}],
	obs_/;MatchQ[obs,{__Function}],NN_Integer,
	UnlabeledVerticesYes_Integer,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 03/06/2026  ***)
(*************************************)
(*Notes: , *)
(*A program to find a schedule of external field values for use in parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach.,

Inputs are:,
1. seedGraph_List - adjacency matrix of the input seed graph,
2. bt_Real = a fixed value of the inverse temperature,
3. hamiltonian_[hparams___] =  the hamiltonian with its parameters that assigns graphs energy,
4. delH_[delHparams___] = function that gives the change in energy when a single edge is flipped,
5. externalFieldTable_List = a list of real number values of the external field for use in the initial run (e.g. a table of H values for the Ising model),
6. conjugateObs_ = a function, the conjugate observable to the external field (e.g. the magnetization formula in the Ising model),
7. obs = the observable quantity(ies) to be measured,
8. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
9. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,
10. numReplicas_Integer = the number of replicas that will be used for parallel tempering. default is set to the number of kernels,

Outputs a list with the following elements:,
1. the external field schedule (a list containing the external field values, 
	edge list for swaps, and an association matching external field values to edge labels),
2. an association of the output from the MultiHistogram run)

*)

Module[
{result,hist,minusbetaFAssn,numVars=Length[conjugateObs],hMins,hMaxs,
	chartConjugateField,numInterpolatingSamples=20,pSoft=0.98,h,hVars,
	sigmaInterpolationAsn,sigmaMat,metricDistance,rho,densityPoints,dist,
	samples,centers,edges,neighbors,numNeighbors,replicasDistMat,
	replicaLabels},

PrintTemporary[" Running GraphCTLScheduler for multiple external fields with input external 
	fields",externalFieldTable," at beta ",bt," numReplicas ",numReplicas];

hist=ECGrav`GraphMultiHistogram[seedGraph,bt,hamiltonian[hparams],delH[delHparams],
		externalFieldTable,conjugateObs,obs,NN,UnlabeledVerticesYes];


minusbetaFAssn=hist[[2,All,"minusBetaF"]];
chartConjugateField=hist[[2,All,"data",All,4]];

hMins=Min/@Transpose[Keys[minusbetaFAssn]];
hMaxs=Max/@Transpose[Keys[minusbetaFAssn]];


(* ===============1) BUILD the thermodynamic metric \[CapitalSigma](H1,H2,...)================*)

hVars=Array[h,numVars];

(*Create interpolating functions from samples of MBAR data for spead*)
PrintTemporary[" Collecting ",numInterpolatingSamples^numVars," sample points using MBAR for 
	interpolation of thermodynamic metric"];


sigmaInterpolationAsn=<|
	Table[
		Table[
			{i,j}->(Thread[{hVars,(hMins-0.2*Abs[hMaxs-hMins]),(hMaxs+0.2*Abs[hMaxs-hMins]),
						(hMaxs-hMins)/numInterpolatingSamples}]
				/.{y__}:>Join@@ParallelTable[
					{hVars,Max[0.0,ECGrav`ExtrapolatedExpectationValue[bt,hVars,minusbetaFAssn,
						chartConjugateField,chartConjugateField[[All,All,i]]
						*chartConjugateField[[All,All,j]]]]},y,
						DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}])
			,{i,1,j}]
	,{j,1,numVars}]|>;


sigmaInterpolationAsn=Interpolation/@sigmaInterpolationAsn;


(* Build Sigma Matrix, S_{a,b}=b^2Cov(M^a,M^b)*)
sigmaMat[x__?NumericQ]:=With[
	{mat=Table[
		Table[sigmaInterpolationAsn[[Key[Sort[{i,j}]]]][x],{i,1,numVars}]
		,{j,1,numVars}]},

	bt^2*mat
];


(* Density function*)
(*soften extreme spikes so one region doesn't eat all points.Set p=1 for pure rho;
	p<1 flattens.*)


rho[x__?NumericQ]:=(Sqrt[Max[Det[sigmaMat[x]],0.0]])^(pSoft);

(* Metric distance function 
	(*Metric distance squared between two field points using midpoint metric:d^2=\[CapitalDelta]H^T (\[Beta]^2 \[CapitalSigma](mid)) \[CapitalDelta]H*)*)

metricDistance[x_List,y_List]:=With[{del=x-y,mid=(x+y)/2.0},

	Sqrt[Max[0.0,del . (sigmaMat@@mid) . del]]

];

(* =================2) Sample Thermodynamic surface ===================*)
PrintTemporary["Sampling thermodynamic density..."];

(*Create a fine grid of points and their weights*)

densityPoints=
	Thread[{hVars,hMins,hMaxs,(hMaxs-hMins)/20}]
		/.{y__}:>Flatten[(ParallelTable[{hVars,rho@@hVars},y,
							DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]),1];

Print[" rho = Sqrt[det[g]] plot ",ListPlot3D[Flatten/@densityPoints,Mesh->All]];
(*Print[" densityPoints ListPlot ",ListPointPlot3D[Flatten/@densityPoints]];*)

(*Create the continuous distribution using weights*)
dist=SmoothKernelDistribution[WeightedData[densityPoints[[All,1]],densityPoints[[All,2]]]];

Print["Smoothened rho contours ",Thread[{hVars,hMins,hMaxs}]
	/.{y__}:>ContourPlot[PDF[dist,hVars],y,PlotLegends->Automatic]];

(*2. Sample from the continuous distribution*)
PrintTemporary["Sampling 20,000*numConjFields points from distribution..."];
samples=RandomVariate[dist,20000*numVars];


(*Clip samples to ensure they stay within physical bounds*)
samples=Select[samples,(And@@Thread[hMins<=#<=hMaxs])&];

(*3. Use FindClusters to find Centroidal Voronoi Centers*)
PrintTemporary["Clustering into ",numReplicas," replicas..."];

(*KMeans clustering finds the centers that minimize the variance*)
centers=Mean/@FindClusters[samples,numReplicas,
		Method->{"KMeans","MaxIterations"->200,"Standardize"->False},
		DistanceFunction->metricDistance];

centers=Select[centers,(And@@Thread[hMins<=#<=hMaxs])&];

centers=LexicographicSort[Sort/@centers,NumericalOrder];


(* ==========4) DEFINE SWAP NEIGHBORS USING THE LOCAL METRIC=============*)
(*Build a neighbor list connect each center to its nn nearest neighbors by metric 
	distance*)
numNeighbors=numVars*2; (*typical twice the dimension*)
replicasDistMat=Table[If[i==j,Infinity,-1.0],{i,numReplicas},{j,numReplicas}];

Do[
	Do[replicasDistMat[[i,j]]=replicasDistMat[[j,i]]=metricDistance[centers[[i]],centers[[j]]]
	,{i,1,j-1 }]
,{j,1,numReplicas}];


neighbors=Table[Ordering[replicasDistMat[[i]],numNeighbors],{i,numReplicas}];


(*Edge list for swaps (undirected,unique) Note that the swap graph is not guaranteed to be regular*)
edges=DeleteDuplicates@(
	Sort/@Flatten[Table[UndirectedEdge[i,j],{i,numReplicas},{j,neighbors[[i]]}],1]);


replicaLabels=<|Table[i->centers[[i]],{i,numReplicas}]|>;

result={{centers,edges,replicaLabels},hist};

Remove[hist,minusbetaFAssn,hMins,hMaxs,chartConjugateField,numVars,
	numInterpolatingSamples,pSoft,h,hVars,sigmaInterpolationAsn,sigmaMat,metricDistance,
	rho,densityPoints,dist,samples,centers,edges,neighbors,numNeighbors,
	replicasDistMat,replicaLabels];

result
];


(*Overload without delH *)
ECGrav`GraphCTLSchedule[seedGraph_List,bt_Real,hamiltonian_[hparams___],externalFieldTable_List/;(MatrixQ[externalFieldTable]&&Length[externalFieldTable[[1]]]>1),conjugateObs_/;MatchQ[conjugateObs,{__Function}],obs_/;MatchQ[obs,{__Function}],NN_Integer,UnlabeledVerticesYes_Integer,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 02/25/2026  ***)
(*************************************)
(*Notes: , *)
(*A program to find a schedule of external field values for use in parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach.,

Inputs are:,
1. seedGraph_List - adjacency matrix of the input seed graph,
2. bt_Real = a fixed value of the inverse temperature,
3. hamiltonian_[hparams___] =  the hamiltonian with its parameters that assigns graphs energy,
4. externalFieldTable_List = a list of real number values of the external field for use in the initial run (e.g. a table of H values for the Ising model),
5. conjugateObs_ = a function, the conjugate observable to the external field (e.g. the magnetization formula in the Ising model),
6. obs = the observable quantity(ies) to be measured,
7. NN = number of independent sweeps (so that actual number of sweeps is correlation time times NN),
8. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,
9. numReplicas_Integer = the number of replicas that will be used for parallel tempering. default is set to the number of kernels,

Outputs a list with the following elements:,
1. the external field schedule (a list containing the external field values, 
	edge list for swaps, and an association matching external field values to edge labels),
2. an association of the output from the MultiHistogram run)

*)

Module[
{result,hist,minusbetaFAssn,numVars=Length[conjugateObs],hMins,hMaxs,
	chartConjugateField,numInterpolatingSamples=20,pSoft=0.98,h,hVars,
	sigmaInterpolationAsn,sigmaMat,metricDistance,rho,densityPoints,dist,
	samples,centers,edges,neighbors,numNeighbors,replicasDistMat,
	replicaLabels},

PrintTemporary[" Running GraphCTLScheduler for multiple external fields with input external 
	fields",externalFieldTable," at beta ",bt," numReplicas ",numReplicas];

hist=ECGrav`GraphMultiHistogram[seedGraph,bt,hamiltonian[hparams],
		externalFieldTable,conjugateObs,obs,NN,UnlabeledVerticesYes];


minusbetaFAssn=hist[[2,All,"minusBetaF"]];
chartConjugateField=hist[[2,All,"data",All,4]];

hMins=Min/@Transpose[Keys[minusbetaFAssn]];
hMaxs=Max/@Transpose[Keys[minusbetaFAssn]];


(* ===============1) BUILD the thermodynamic metric \[CapitalSigma](H1,H2,...)================*)

hVars=Array[h,numVars];

(*Create interpolating functions from samples of MBAR data for spead*)
PrintTemporary[" Collecting ",numInterpolatingSamples^numVars," sample points using MBAR for 
	interpolation of thermodynamic metric"];


sigmaInterpolationAsn=<|
	Table[
		Table[
			{i,j}->(Thread[{hVars,(hMins-0.2*Abs[hMaxs-hMins]),(hMaxs+0.2*Abs[hMaxs-hMins]),
						(hMaxs-hMins)/numInterpolatingSamples}]
				/.{y__}:>Join@@ParallelTable[
					{hVars,Max[0.0,ECGrav`ExtrapolatedExpectationValue[bt,hVars,minusbetaFAssn,
						chartConjugateField,chartConjugateField[[All,All,i]]
						*chartConjugateField[[All,All,j]]]]},y,
						DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}])
			,{i,1,j}]
	,{j,1,numVars}]|>;


sigmaInterpolationAsn=Interpolation/@sigmaInterpolationAsn;


(* Build Sigma Matrix, S_{a,b}=b^2Cov(M^a,M^b)*)
sigmaMat[x__?NumericQ]:=With[
	{mat=Table[
		Table[sigmaInterpolationAsn[[Key[Sort[{i,j}]]]][x],{i,1,numVars}]
		,{j,1,numVars}]},

	bt^2*mat
];


(* Density function*)
(*soften extreme spikes so one region doesn't eat all points.Set p=1 for pure rho;
	p<1 flattens.*)


rho[x__?NumericQ]:=(Sqrt[Max[Det[sigmaMat[x]],0.0]])^(pSoft);

(* Metric distance function 
	(*Metric distance squared between two field points using midpoint metric:d^2=\[CapitalDelta]H^T (\[Beta]^2 \[CapitalSigma](mid)) \[CapitalDelta]H*)*)

metricDistance[x_List,y_List]:=With[{del=x-y,mid=(x+y)/2.0},

	Sqrt[Max[0.0,del . (sigmaMat@@mid) . del]]

];

(* =================2) Sample Thermodynamic surface ===================*)
PrintTemporary["Sampling thermodynamic density..."];

(*Create a fine grid of points and their weights*)

densityPoints=
	Thread[{hVars,hMins,hMaxs,(hMaxs-hMins)/20}]
		/.{y__}:>Flatten[(ParallelTable[{hVars,rho@@hVars},y,
							DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]),1];

Print[" rho = Sqrt[det[g]] plot ",ListPlot3D[Flatten/@densityPoints,Mesh->All]];
(*Print[" densityPoints ListPlot ",ListPointPlot3D[Flatten/@densityPoints]];*)

(*Create the continuous distribution using weights*)
dist=SmoothKernelDistribution[WeightedData[densityPoints[[All,1]],densityPoints[[All,2]]]];

Print["Smoothened rho contours ",Thread[{hVars,hMins,hMaxs}]
	/.{y__}:>ContourPlot[PDF[dist,hVars],y,PlotLegends->Automatic]];

(*2. Sample from the continuous distribution*)
PrintTemporary["Sampling 20,000*numConjFields points from distribution..."];
samples=RandomVariate[dist,20000*numVars];


(*Clip samples to ensure they stay within physical bounds*)
samples=Select[samples,(And@@Thread[hMins<=#<=hMaxs])&];

(*3. Use FindClusters to find Centroidal Voronoi Centers*)
PrintTemporary["Clustering into ",numReplicas," replicas..."];

(*KMeans clustering finds the centers that minimize the variance*)
centers=Mean/@FindClusters[samples,numReplicas,
		Method->{"KMeans","MaxIterations"->200,"Standardize"->False},
		DistanceFunction->metricDistance];

centers=LexicographicSort[Sort/@centers,NumericalOrder];


(* ==========4) DEFINE SWAP NEIGHBORS USING THE LOCAL METRIC=============*)
(*Build a neighbor list connect each center to its nn nearest neighbors by metric 
	distance*)
numNeighbors=numVars*2; (*typical twice the dimension*)
replicasDistMat=Table[If[i==j,Infinity,-1.0],{i,numReplicas},{j,numReplicas}];

Do[
	Do[replicasDistMat[[i,j]]=replicasDistMat[[j,i]]=metricDistance[centers[[i]],centers[[j]]]
	,{i,1,j-1 }]
,{j,1,numReplicas}];


neighbors=Table[Ordering[replicasDistMat[[i]],numNeighbors],{i,numReplicas}];


(*Edge list for swaps (undirected,unique) Note that the swap graph is not guaranteed to be regular*)
edges=DeleteDuplicates@(
	Sort/@Flatten[Table[UndirectedEdge[i,j],{i,numReplicas},{j,neighbors[[i]]}],1]);


replicaLabels=<|Table[i->centers[[i]],{i,numReplicas}]|>;

result={{centers,edges,replicaLabels},hist};

Remove[hist,minusbetaFAssn,hMins,hMaxs,chartConjugateField,numVars,
	numInterpolatingSamples,pSoft,h,hVars,sigmaInterpolationAsn,sigmaMat,metricDistance,
	rho,densityPoints,dist,samples,centers,edges,neighbors,numNeighbors,
	replicasDistMat,replicaLabels];

result
];


(* ::Item::Closed:: *)
(*GraphCTLSchedule for multiple external fields at a fixed beta with input from a previous run of Multihistogram*)


(*Primary pattern*)
ECGrav`GraphCTLSchedule[bt_Real,hMins_List/;VectorQ[hMins],hMaxs_List/;VectorQ[hMaxs],minusbetaF_Association,
	conjugateFieldMeasurements_Association,numReplicas_Integer]:=
(*************************************)
(***  Last updated on: 03/05/2026  ***)
(*************************************)
(*Notes: , *)
(*An overload  to find a schedule of external field values for use in parallel tempering that keeps acceptance ratios constant based on the constant thermodynamic length approach. It takes as an input the return from a previous run of GraphMultihistogram ,

Inputs are:,
1. bt_Real = a fixed value of the inverse temperature,
2. hMin_Real = the minimum value of the external field for the schedule to be built,
2. hMax_Real = the maximum value of the external field for the schedule to be built,
3. minusbetaF_Association = an association of external field values and corresponding values of -beta*freeEnergy from a previous run of GraphMultihistograph,
4. conjugateFieldMeasurements_Association = an association of external field values and the corresponding table of the conjugate observable measurements from a previous run of GraphMultihistograph,
5. numReplicas_Integer = the number of replicas that will be used for parallel tempering. default is set to the number of kernels,

Outputs the external field schedule (a list containing the external field values, 
	edge list for swaps, and an association matching external field values to edge labels).*)

Module[{result,numVars=Length[hMins],numInterpolatingSamples=20,pSoft=0.96,h,hVars,
	sigmaInterpolationAsn,sigmaMat,metricDistance,rho,densityPoints,dist,samples,
	centers,edges,swapGraph,neighbors,numNeighbors,replicasDistMat,replicaLabels},


PrintTemporary[" Running GraphCTLScheduler for multiple external fields with input from a 
	previous run of GraphMultihistogram at beta ",bt," with hMins ",hMins," hMaxs "
	,hMaxs];


(* =========================1) BUILD \[CapitalSigma](H1,H2,...)=========================*)

hVars=Array[h,numVars];


(*Create interpolating functions from samples of MBAR data for spead*)
PrintTemporary[" Collecting ",numInterpolatingSamples^numVars," sample points using MBAR for 
	interpolation of Fisher matrix"];


sigmaInterpolationAsn=<|
Table[
	Table[
		{i,j}->(Thread[{hVars,(hMins-0.2*Abs[hMaxs-hMins]),(hMaxs+0.2*Abs[hMaxs-hMins]),(hMaxs-hMins)/numInterpolatingSamples}]
			/.{y__}:>Join@@ParallelTable[
					{hVars, Max[0.0,ECGrav`ExtrapolatedExpectationValue[bt,hVars,
							minusbetaF,conjugateFieldMeasurements,conjugateFieldMeasurements[[All,All,i]]
							*conjugateFieldMeasurements[[All,All,j]]]]},y,
							DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}])
		,{i,1,j}]
,{j,1,numVars}]|>;


sigmaInterpolationAsn=Interpolation/@sigmaInterpolationAsn;


(* Build Sigma Matrix, S_{a,b}=b^2Cov(M^a,M^b)*)
sigmaMat[x__?NumericQ]:=
	With[{mat=Table[
		Table[sigmaInterpolationAsn[[Key[Sort[{i,j}]]]][x],{i,1,numVars}]
	,{j,1,numVars}]},

	bt^2*mat
];


(* Density function*)
(*soften extreme spikes so one region doesn't eat all points.Set p=1 for pure rho;p<1 flattens.*)


rho[x__?NumericQ]:=(Sqrt[Max[Det[sigmaMat[x]],0.0]])^(pSoft);

(* Metric distance function 
(*Metric distance squared between two field points using midpoint metric:d^2=\[CapitalDelta]H^T (\[Beta]^2 \[CapitalSigma](mid)) \[CapitalDelta]H*)*)

metricDistance[x_List,y_List]:=With[{del=x-y,mid=(x+y)/2.0},
	Sqrt[Max[0.0,del . (sigmaMat@@mid) . del]]
];

(* =================2) Sample Thermodynamic surface ===================*)
PrintTemporary["Sampling thermodynamic density..."];

(*Create a fine grid of points and their weights*)

densityPoints=Thread[
	{hVars,hMins,hMaxs,(hMaxs-hMins)/20}]/.{y__}:>Flatten[
		(ParallelTable[{hVars,rho@@hVars},y,DistributedContexts->{$Context,"ECGrav`MCSims`Private`"}]),1];

Print[" rho = Sqrt[det[g]] Plot ",ListPlot3D[Flatten/@densityPoints,Mesh->All]];
(*Print[" densityPoints ListPlot ",ListPointPlot3D[Flatten/@densityPoints]];*)

(*Create the continuous distribution using weights*)
dist=SmoothKernelDistribution[WeightedData[densityPoints[[All,1]],densityPoints[[All,2]]]];

Print["Smooth pdf of rho contours ",
	Thread[{hVars,hMins,hMaxs}]
		/.{y__}:>ContourPlot[PDF[dist,hVars],y,PlotLegends->Automatic]
];

(*2. Sample from the continuous distribution*)
PrintTemporary["Sampling 20,000*numConjFields points from distribution..."];

samples=RandomVariate[dist,20000*numVars];


(*Clip samples to ensure they stay within physical bounds*)
samples=Select[samples,(And@@Thread[hMins<=#<=hMaxs])&];


(*3. Use FindClusters to find Centroidal Voronoi Centers*)
PrintTemporary["Clustering into ",numReplicas," replicas..."];

(*KMeans clustering finds the centers that minimize the variance*)
centers=Mean/@FindClusters[samples,numReplicas,
	Method->{"KMeans","MaxIterations"->200,"Standardize"->False},
	DistanceFunction->metricDistance];

centers=Select[centers,(And@@Thread[hMins<=#<=hMaxs])&];

centers=LexicographicSort[Sort/@centers,NumericalOrder];


(* ==========4) DEFINE SWAP NEIGHBORS USING THE LOCAL METRIC=============*)
(*Build a neighbor list connect each center to its nn nearest neighbors by metric distance*)

numNeighbors=numVars*2; (*typical twice the dimension*)
replicasDistMat=Table[If[i==j,Infinity,-1.0],{i,numReplicas},{j,numReplicas}];

Do[
	Do[replicasDistMat[[i,j]]=replicasDistMat[[j,i]]=metricDistance[centers[[i]],centers[[j]]]
	,{i,1,j-1 }]
,{j,1,numReplicas}];

neighbors=Table[Ordering[replicasDistMat[[i]],numNeighbors],{i,numReplicas}];


(*Edge list for swaps (undirected,unique) Note that the swap graph is not guaranteed to 
	be regular*)
edges=DeleteDuplicates@(
	Sort/@Flatten[
		Table[UndirectedEdge[i,j],{i,numReplicas},{j,neighbors[[i]]}],1]);

swapGraph=Graph[Range[numReplicas],edges,VertexCoordinates->Thread[Range[numReplicas]->centers]];


replicaLabels=<|Table[i->centers[[i]],{i,numReplicas}]|>;

result={centers,edges,replicaLabels};

Remove[numVars,numInterpolatingSamples,pSoft,h,hVars,sigmaInterpolationAsn,sigmaMat,
	metricDistance,rho,densityPoints,dist,samples,centers,edges,swapGraph,neighbors,
	numNeighbors,replicasDistMat,replicaLabels];

result
];


(* ::Item::Closed:: *)
(*GraphCTLSchedule Catch-all*)


(* Catch-all Pattern *)
ECGrav`GraphCTLSchedule[args___]:=(Message[ECGrav`GraphCTLSchedule::argerr, args];
$Failed);


(* ::Subsection:: *)
(*GraphParallelTempering*)


(* ::Item::Closed:: *)
(*GraphParallelTempering Primary*)


ECGrav`GraphParallelTempering[seedGraph_List, btTable_List,
	hamiltonian_[hparams___],delH_[delHparams___],obs_,EnergyOrMag_Integer,
	NN_Integer,UnlabeledVerticesYes_Integer,minEtoBeat_Real]:=

(*************************************)
(***  Last updated on: 02/16/2026  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs,
*3. 02/16/2026 Update - made code more efficient by exporting the computeMBF and cleaning up the swap functions*)
(*
Implements Parallel tempering algorithm on graph models at several different temperatures. It does not determine the temperature schedule here, it assumes the input beta list comes from a scheduler such as Constant Entropy Increase (CEI) or Constant Thermodynamic Length (CTL) temperature schedule. It equilibriates, computes correlation time, and then applies measurements. It does temperature swaps during the measurement step. It is parallelized so that equilibriation, computation of correlation time, and sweeps during measurement are all done in parallel. Temperature swaps are done on the master kernel. ,

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
9. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with four objects:, 1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. the final state of each replicas ,
4. an association of each replica and its beta history ,
*)

Module[{result,vCount=Length[seedGraph],groundStates,histories,maxGStateCount=500,
	replicas,replicaKeysOrderedByBeta,numRep,bt,minStates,candminE,measurements,
	numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},

groundStates=<||>;
histories=<||>;
numRep=Length[btTable];


(*(******************************)
(**      Equilibriate          **)
(******************************)*)
Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,btTable[[i]],
		hamiltonian[hparams],delH[delHparams],UnlabeledVerticesYes],{i,Length[btTable]}
	,DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Tempoutput[[All,2]];


(* Extract the minimum energy and corresponding states found from the equilibriation 
	run and the temperature schedule run*)
AppendTo[groundStates,"minEnergy"->Min[Tempoutput[[All,1,"minEnergy"]]]];

AppendTo[groundStates,"minEstates"->Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==groundStates[[Key["minEnergy"]]]&][[All,"minEstates"]]]];

(*compare the minimum energy from the equilibriation run with that of the temperature 
schedule run and reset the minimum energy and states *)

If[minEtoBeat<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=minEtoBeat;
	groundStates[[Key["minEstates"]]]={}
];

(*(******************************)
(*  Compute Correlation times *)
(******************************)*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
	ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]]
		,replicas[[Key[i],Key["beta"]]],hamiltonian[hparams],delH[delHparams]
		,replicas[[Key[i],Key["eqlT"]]],locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes]
		,{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}
		]
	]
];


(*Update replicas*)
Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];


(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


If[candminE<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=candminE;
	groundStates[[Key["minEstates"]]]=
	Union@@
		Values[
			Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]
			],
	If[candminE==groundStates[[Key["minEnergy"]]]&&Length[groundStates[[Key["minEstates"]]]]<=maxGStateCount,
		groundStates[[Key["minEstates"]]]=Union[groundStates[[Key["minEstates"]]]
		,Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==candminE&][[All,"minEstates"]]]
		];
	];
];


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, 
	"history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;


Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Module[{thisReplInd,nextReplInd,deltaBetaDeltaE,accept,tempThisReplicaState},

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	deltaBetaDeltaE= (replicas[[Key[thisReplInd],Key["beta"]]]-replicas[[Key[nextReplInd],Key["beta"]]])*(replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]-replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]);

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;
	histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;


accept=0;
If[deltaBetaDeltaE>0,accept = 1,
	If[RandomReal[]<Exp[deltaBetaDeltaE],accept =1];
];


If[accept==1, (*Do swap of replicas*)

	(*Update states: store the state of the high temp replica temporarily*)

	tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
	replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
	replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

	histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
	histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

	histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
	histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

	(*Increase the number of swap accepts by one for both replicas*)
	histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
	histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

];

,{repl,1,numRep-1,2}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)

numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];

measurements=Reap[
While[numsweeps<=NN,

If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
	ParallelTable[
		<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]]
			,locrepl[[Key[i],Key["beta"]]],hamiltonian[hparams],delH[delHparams]
			,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,

	{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
	]
];


(*Update replicas*)
Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

(*Update minimum energy*)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

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

swap[numsweeps];

If[Mod[numsweeps,1]==0,
	Table[
		bt=replicas[[Key[i],Key["beta"]]];
		Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]]
		,Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
	,{i,numRep}];
];

numsweeps++;

]][[2]];


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[btTable]]]}]|>;

Do[histories[[Key[i],Key["history"]]]=(histories[[Key[i],Key["history"],
	-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]]),{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)

result={ groundStates,chart,replicas,histories};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,
	histories,numRep,bt,minStates,candminE,measurements,numsweeps,Tempoutput,swap,
	printCase,repNumSweeps,chart];

result

];


ECGrav`GraphParallelTempering[seedGraph_List, btTable_List,
	hamiltonian_[hparams___],obs_,EnergyOrMag_Integer,NN_Integer,UnlabeledVerticesYes_Integer,minEtoBeat_Real]:=

(*************************************)
(***  Last updated on: 02/16/2026  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs,
*3. 02/16/2026 Update - made code more efficient by exporting the computeMBF and cleaning up the swap functions*)
(*
Implements Parallel tempering algorithm on graph models at several different temperatures. It does not determine the temperature schedule here, it assumes the input beta list comes from a scheduler such as Constant Entropy Increase (CEI) or Constant Thermodynamic Length (CTL) temperature schedule. It equilibriates, computes correlation time, and then applies measurements. It does temperature swaps during the measurement step. It is parallelized so that equilibriation, computation of correlation time, and sweeps during measurement are all done in parallel. Temperature swaps are done on the master kernel. ,

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

Outputs a list with four objects:, 1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. the final state of each replicas ,
4. an association of each replica and its beta history ,
*)

Module[{result,vCount=Length[seedGraph],groundStates,histories,maxGStateCount=500,
	replicas,replicaKeysOrderedByBeta,numRep,bt,minStates,candminE,measurements,
	numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},

groundStates=<||>;
histories=<||>;
numRep=Length[btTable];


(*(******************************)
(**      Equilibriate          **)
(******************************)*)
Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,btTable[[i]],hamiltonian[hparams]
		,UnlabeledVerticesYes],{i,Length[btTable]},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];


(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

(* Extract the minimum energy and corresponding states found from the equilibriation 
run and the temperature schedule run*)
AppendTo[groundStates,"minEnergy"->Min[Tempoutput[[All,1,"minEnergy"]]]];

AppendTo[groundStates,"minEstates"->Union@@Values[Select[Tempoutput[[All,1]],#[[Key["minEnergy"]]]==groundStates[[Key["minEnergy"]]]&][[All,"minEstates"]]]];

(*compare the minimum energy from the equilibriation run with that of the temperature 
schedule run and reset the minimum energy and states *)

If[minEtoBeat<groundStates[[Key["minEnergy"]]],
	groundStates[[Key["minEnergy"]]]=minEtoBeat;
	groundStates[[Key["minEstates"]]]={}
];


(*(******************************)
(*  Compute Correlation times *)
(******************************)*)

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[i->ECGrav`GraphComputeCorrelationTime[locrepl[[Key[i],Key["state"],Key["graph"]]]
		,replicas[[Key[i],Key["beta"]]],hamiltonian[hparams],replicas[[Key[i],Key["eqlT"]]]
		,locMinEtoBeat,EnergyOrMag,UnlabeledVerticesYes],{i,numRep}
		,DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
	]
];


(*Update replicas*)
Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];

(* Extract the minimum energy and corresponding states *)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0
	, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;


Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Module[{thisReplInd,nextReplInd,deltaBetaDeltaE,accept,tempThisReplicaState},

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	deltaBetaDeltaE= (replicas[[Key[thisReplInd],Key["beta"]]]-replicas[[Key[nextReplInd],Key["beta"]]])*(replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]-replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]);

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;
	histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;

	accept=0;
	If[deltaBetaDeltaE>0,accept = 1,
		If[RandomReal[]<Exp[deltaBetaDeltaE],accept =1];
	];


	If[accept==1, (*Do swap of replicas*)

	(*Update states: store the state of the high temp replica temporarily*)

		tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
		replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
		replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

	(*Update histories*)

		histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
		histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

		histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
		histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

		(*Increase the number of swap accepts by one for both replicas*)
		histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
		histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

	];

	,{repl,1,numRep-1,2}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];

measurements=Reap[
	While[numsweeps<=NN,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		Tempoutput=Association[
			With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},

		ParallelTable[
			<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]]
			,locrepl[[Key[i],Key["beta"]]],hamiltonian[hparams]
			,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>
			,{i,numRep}
			,DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
			]
		];

	(*Update replicas*)
	Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

	(*Update minimum energy*)
	candminE=Min[Tempoutput[[All,1,"minEnergy"]]];


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

	swap[numsweeps];

	If[Mod[numsweeps,1]==0,
		Table[
			bt=replicas[[Key[i],Key["beta"]]];Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
		,{i,numRep}];
	];

	numsweeps++;

]][[2]];


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[btTable]]]}]|>;


Do[histories[[Key[i],Key["history"]]]=(histories[[Key[i],Key["history"]
		,-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]]),{i,numRep}];
		(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)

result={ groundStates,chart,replicas,histories};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,
	histories,numRep,bt,minStates,candminE,measurements,numsweeps,Tempoutput,swap,
	printCase,repNumSweeps,chart];

result

];


(* ::Item:: *)
(*Overload GraphParallelTempering with previous beta replica as input*)


ECGrav`GraphParallelTempering[inputReplicas_Association/;(Length[DeleteDuplicates[inputReplicas[[All,"beta"]]]]>1&&Length[DeleteDuplicates[inputReplicas[[All,"externalField"]]]]==1)
	,hamiltonian_[hparams___],delH_[delHparams___],obs_
	,NN_Integer,UnlabeledVerticesYes_Integer,minEtoBeat_Real]:=

(*************************************)
(***  Last updated on: 02/16/2026  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the sweep is done with labeled or unlabeled graphs,
*3. 02/16/2026 Update - made code more efficient by exporting the computeMBF and cleaning up the swap functions*)
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

Outputs a list with four objects:, 
1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. the final state of each replicas ,
4. an association of each replica and its beta history ,

*)

Module[{result,numRep=Length[inputReplicas],vCount=Length[inputReplicas[[1,"state","graph"]]]
	,btTable=Values[inputReplicas[[All,"beta"]]],replicas=inputReplicas
	,groundStates=<|"minEnergy"->minEtoBeat,"minEstates"->{}|>,histories=<||>
	,maxGStateCount=500,replicaKeysOrderedByBeta,bt,minStates,candminE,measurements
	,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Module[{thisReplInd,nextReplInd,deltaBetaDeltaE,accept,tempThisReplicaState},

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	deltaBetaDeltaE= (replicas[[Key[thisReplInd],Key["beta"]]]-replicas[[Key[nextReplInd],Key["beta"]]])*(replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]-replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]);

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;
	histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;

	accept=0;
	If[deltaBetaDeltaE>0,accept = 1,
		If[RandomReal[]<Exp[deltaBetaDeltaE],accept =1];
	];


	If[accept==1, (*Do swap of replicas*)

		(*Update states: store the state of the high temp replica temporarily*)

		tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
		replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
		replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

		(*Update histories*)

		histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
		histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

		histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
		histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

		(*Increase the number of swap accepts by one for both replicas*)
		histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
		histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

	];


,{repl,1,numRep-1,2}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
While[numsweeps<=NN,

	If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

	Tempoutput=Association[
		With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
			ParallelTable[
				<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]]
				,locrepl[[Key[i],Key["beta"]]],hamiltonian[hparams],delH[delHparams]
				,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,
			{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
		]
	];

	(*Update replicas*)
	Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

	(*Update minimum energy*)
	candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

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

	swap[numsweeps];

	If[Mod[numsweeps,1]==0,
		Table[
			bt=replicas[[Key[i],Key["beta"]]];Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]],Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
		,{i,numRep}];
	];

	numsweeps++;

]][[2]];

chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[btTable]]]}]|>;

Do[histories[[Key[i],Key["history"]]]=(histories[[Key[i],Key["history"]
	,-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]]),{i,numRep}];
	(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)

result={ groundStates,chart,replicas,histories};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,histories
	,numRep,bt,minStates,candminE,measurements,numsweeps,Tempoutput,swap,printCase
	,repNumSweeps,chart];

result

];


ECGrav`GraphParallelTempering[inputReplicas_Association/;(Length[DeleteDuplicates[inputReplicas[[All,"beta"]]]]>1&&Length[DeleteDuplicates[inputReplicas[[All,"externalField"]]]]==1)
	,minEtoBeat_Real,hamiltonian_[hparams___],obs_,NN_Integer
	,UnlabeledVerticesYes_Integer]:=

(*************************************)
(***  Last updated on: 02/16/2026  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the sweep is done with labeled or unlabeled graphs,
*3. 02/16/2026 Update - made code more efficient by exporting the computeMBF and cleaning up the swap functions*)
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

Outputs a list with four objects:, 1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. the final state of each replicas ,
4. an association of each replica and its beta history ,

*)

Module[{result,numRep=Length[inputReplicas],vCount=Length[inputReplicas[[1,"state","graph"]]]
	,btTable=Values[inputReplicas[[All,"beta"]]],replicas=inputReplicas
	,groundStates=<|"minEnergy"->minEtoBeat,"minEstates"->{}|>,histories=<||>
	,maxGStateCount=500,replicaKeysOrderedByBeta,bt,minStates,candminE,measurements
	,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart},


replicaKeysOrderedByBeta=Keys[Sort[replicas[[All,"beta"]]]];

histories=<|Table[i-><|"beta"->replicas[[Key[i],Key["beta"]]],"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["beta"]]],{i,numRep}];

swap[offset_Integer]:=(*attempts swaps between replicas according to weight exp(-(delta b)(delta E)). *)
Module[{thisReplInd,nextReplInd,deltaBetaDeltaE,accept,tempThisReplicaState},

Do[
	thisReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl-1,numRep]+1]];
	nextReplInd=replicaKeysOrderedByBeta[[Mod[offset+repl,numRep]+1]];

	deltaBetaDeltaE= 
		(replicas[[Key[thisReplInd],Key["beta"]]]-replicas[[Key[nextReplInd],Key["beta"]]])
		*(replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]-replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]);

	(*Increase the number of swap tries by one for both replicas*)

	histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;
	histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;


	accept=0;
	If[deltaBetaDeltaE>0,accept = 1,
		If[RandomReal[]<Exp[deltaBetaDeltaE],accept =1];
	];


	If[accept==1, (*Do swap of replicas*)

		(*Update states: store the state of the high temp replica temporarily*)

		tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
		replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
		replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

		(*Update histories*)

		histories[[Key[thisReplInd],Key["history"],1]]=histories[[Key[nextReplInd],Key["history"],-1]];
		histories[[Key[nextReplInd],Key["history"],1]]=histories[[Key[thisReplInd],Key["history"],-1]];

		histories[[Key[thisReplInd],Key["history"]]]=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
		histories[[Key[nextReplInd],Key["history"]]]=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

		(*Increase the number of swap accepts by one for both replicas*)
		histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
		histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

	];

,{repl,1,numRep-1,2}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
While[numsweeps<=NN,

If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

Tempoutput=Association[
	With[{locrepl=replicas,locMinEtoBeat=groundStates[[Key["minEnergy"]]]},
		ParallelTable[
			<|i->ECGrav`GraphSweepReplica[locrepl[[Key[i],Key["state"],Key["graph"]]]
			,locrepl[[Key[i],Key["beta"]]],hamiltonian[hparams]
			,locrepl[[Key[i],Key["corrT"]]],locMinEtoBeat,UnlabeledVerticesYes]|>,
		{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
	]
];

(*Update replicas*)
Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i]]][[2]],{i,numRep}];

(*Update minimum energy*)
candminE=Min[Tempoutput[[All,1,"minEnergy"]]];

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


swap[numsweeps];

If[Mod[numsweeps,1]==0,
	Table[
		bt=replicas[[Key[i],Key["beta"]]];
		Sow[Flatten[{numsweeps,bt,replicas[[Key[i],Key["state"],Key["energy"]]]
			,Through[obs [replicas[[Key[i],Key["state"],Key["graph"]]]]]}],bt]
	,{i,numRep}];
];

numsweeps++;

]][[2]];

chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;Length[btTable]]]}]|>;

Do[histories[[Key[i],Key["history"]]]=
	(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]])
	,{i,numRep}];
(*remove the -1's in the initiation and select every nth point so that the total length is no more than numberOfDataPoints*)

result={groundStates,chart,replicas,histories};

Remove[vCount,groundStates,maxGStateCount,replicas,replicaKeysOrderedByBeta,histories,numRep,bt,minStates,candminE,measurements,numsweeps,Tempoutput,swap,printCase,repNumSweeps,chart];

result

];


(* ::Item::Closed:: *)
(*Overload GraphParallelTempering with one or more  external fields at a fixed beta and input external field  matrix*)


(*Primary pattern*)
ECGrav`GraphParallelTempering[seedGraph_List, bt_Real,hamiltonian_[hparams___],
	delH_[delHparams___],replicaSwapEdgesLabels_List,
	conjugateObs_/;MatchQ[conjugateObs,{__Function}],
	obs_/;MatchQ[obs,{__Function}],
	NN_Integer,UnlabeledVerticesYes_Integer,minEtoBeat_Real]:=

(*************************************)
(***  Last updated on: 04/04/2026  ***)
(*************************************)
(*Notes: ,*)
(*
Implements Parallel tempering algorithm on graph models at a fixed inverse temperature beta and at several different values of external field parameter. It does not determine the external field schedule here, it assumes the input external field list comes from a scheduler such as  determined by Constant Entropy Increase (CEI) or Constant Thermodynamic Length (CTL) temperature schedule. It equilibriates, computes correlation time, and then applies measurements. It does temperature swaps during the measurement step. It is parallelized so that equilibriation, computation of correlation time, and sweeps during measurement are all done in parallel. Temperature swaps are done on the master kernel. ,

Depends on the functions: GraphSweepReplicas, CorrelationTime, GraphEquilibriate, GraphComputeCorrelationTime., 

Inputs are:, 
1. seedGraph - adjacency matrix of the input seed graph,
2. btTable - a list of inverse temperatures for the replicas,
3. minEtoBeat -  seed minimum energy at or below which to save graphs,
4. hamiltonian =  function that assigns graphs energy,
5. delH = function that gives the change in energy when a single edge is flipped,
6. replicaSwapEdgesLabels = A list with two elements, the replica swap graph edges 
	and an association mapping graph vertices to external field values 
7. conjugateObs = a list of functions corresponding to the conjugate variables to the external fields
8. obs = the observable quantity in question,
9. NN = number of independent sweeps to be carried out(so that actual number of sweeps is correlation time times NN),
10. numberOfDataPoints = number of data points of measurements of the observables to be returned.,
11. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with four objects:, 
1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. the final state of each replicas ,
4. an association of each replica and its beta history ,
*)

Module[{result,groundStates,histories,maxGStateCount=500,replicas,numRep,
	minStates,candminE,measurements,numsweeps,Tempoutput,
	ChooseRandomIndependentEdgeSet,Swap,printCase,repNumSweeps,chart},

groundStates=<||>;
histories=<||>;
numRep=Length[replicaSwapEdgesLabels[[2]]];

PrintTemporary["Running PT for graph with numRep ",numRep," at beta ",bt,
	" replicaSwapEdgesLabels ",replicaSwapEdgesLabels];

(*(******************************)
(**      Equilibriate          **)
(******************************)*)

Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,bt,
		Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
		Apply[delH,replicaSwapEdgesLabels[[2,Key[i]]]],UnlabeledVerticesYes]
	,{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

groundStates=Tempoutput[[All,1]];

(*(******************************)
(*  Compute Correlation times *)
(******************************)*)

Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphComputeCorrelationTime[
		replicas[[Key[i],Key["state"],Key["graph"]]],bt,
		Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
		Apply[delH,replicaSwapEdgesLabels[[2,Key[i]]]],
		replicas[[Key[i],Key["eqlT"]]],
		groundStates[[Key[i],Key["minEnergy"]]],conjugateObs,UnlabeledVerticesYes]
		,{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];


(*Update replicas*)
Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];

(* Extract the minimum energy and corresponding states *)
Do[
	candminE=Tempoutput[[Key[i],1,"minEnergy"]];


	If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
		groundStates[[Key[i],Key["minEnergy"]]]=candminE;
		groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
		If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
			groundStates[[Key[i],Key["minEstates"]]]=Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
		];
	]

,{i,numRep}];


histories=<|Table[
			i-><|"externalField"->replicas[[Key[i],Key["externalField"]]],
			"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>
			,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["externalField"]]],{i,numRep}];

(*(****************************)
(*  Define the replica swap *)
(****************************)*)

(*Helper function to choose a set of non-intersecting edges*)
ChooseRandomIndependentEdgeSet[edgeList_List]:=
Module[{output={},remaining=edgeList,randomEdge},
	While[remaining!={},
		randomEdge=RandomChoice[remaining];
		output=Join[output,{randomEdge}];
		remaining=Select[remaining,!(MemberQ[#,randomEdge[[1]]]||MemberQ[#,randomEdge[[2]]])&]
	];
	output
];

Swap[]:=(*attempts swaps between replicas according to weight exp(-beta(delta H)(delta E)). *)
Module[{replicaSwapMatchings,thisReplInd,nextReplInd,betaDelHDelConj,accept,
	tempThisReplicaState},


	replicaSwapMatchings=ChooseRandomIndependentEdgeSet[replicaSwapEdgesLabels[[1]]];


	Do[
		thisReplInd=edgeIndex[[1]];
		nextReplInd=edgeIndex[[2]];

		betaDelHDelConj=bt*(replicas[[Key[thisReplInd],Key["externalField"]]]
							-replicas[[Key[nextReplInd],Key["externalField"]]])
							. (Through[conjugateObs[replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]]]
								-Through[conjugateObs[replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]]]
							);

		(*Increase the number of swap tries by one for both replicas*)

		histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;
		histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;

		accept=0;
		If[betaDelHDelConj>0,accept = 1,
			If[RandomReal[]<Exp[betaDelHDelConj],accept =1];
		];


		If[accept==1, (*Do swap of replicas*)

			(*Update states: store the state of the high temp replica temporarily*)

			tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
			replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
			replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

			(*Update the energies to reflect the new external fiel dvalues*)
			replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]},replicas[[Key[thisReplInd],Key["externalField"]]]]];
			replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]},replicas[[Key[nextReplInd],Key["externalField"]]]]];

			(*Update histories*)

			histories[[Key[thisReplInd],Key["history"],1]]
				=histories[[Key[nextReplInd],Key["history"],-1]];
			histories[[Key[nextReplInd],Key["history"],1]]
				=histories[[Key[thisReplInd],Key["history"],-1]];

			histories[[Key[thisReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
			histories[[Key[nextReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

			(*Increase the number of swap accepts by one for both replicas*)
			histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
			histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

		];


	,{edgeIndex,replicaSwapMatchings}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)

numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
	While[numsweeps<=NN,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		Tempoutput=Association[
			ParallelTable[
				<|i->ECGrav`GraphSweepReplica[
						replicas[[Key[i],Key["state"],Key["graph"]]],bt,
						Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
						Apply[delH,replicaSwapEdgesLabels[[2,Key[i]]]],
						replicas[[Key[i],Key["corrT"]]],
						groundStates[[Key[i],Key["minEnergy"]]],
						UnlabeledVerticesYes]
				|>,

			{i,numRep},
			DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
		];

		(*Update replicas*)
		Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i],2]],{i,numRep}];

		(*Update minimum energy*)
		Do[
			candminE=Tempoutput[[Key[i],1,"minEnergy"]];

			If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
				groundStates[[Key[i],Key["minEnergy"]]]=candminE;
				groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
				If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
					groundStates[[Key[i],Key["minEstates"]]]
						=Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
				];
			];

		,{i,numRep}];

		Swap[];

		If[Mod[numsweeps,1]==0,
			Table[
				Sow[{numsweeps,replicaSwapEdgesLabels[[2,Key[i]]],
					replicas[[Key[i],Key["state"],Key["energy"]]],
					Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					replicas[[Key[i],Key["state"],Key["energy"]]]
						-replicaSwapEdgesLabels[[2,Key[i]]]
							. Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					Through[obs[replicas[[Key[i],Key["state"],Key["graph"]]]]]},i]
			,{i,numRep}];
		];

		numsweeps++;

	]
][[2]];


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;numRep]]}]|>;

(*remove the -1's in the initiation*)
Do[If[histories[[Key[i],Key["swapAccept"]]]<Length[histories[[Key[i],Key["history"]]]],
	histories[[Key[i],Key["history"]]]
		=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]])]
,{i,numRep}];

result={ groundStates,chart,replicas,histories};

Remove[groundStates,maxGStateCount,replicas,histories,numRep,minStates,candminE,measurements,numsweeps,Tempoutput,ChooseRandomIndependentEdgeSet,Swap,printCase,repNumSweeps,chart];

result

];


(*Overload pattern without delH *)
ECGrav`GraphParallelTempering[seedGraph_List, bt_Real,hamiltonian_[hparams___],
	replicaSwapEdgesLabels_List,
	conjugateObs_/;MatchQ[conjugateObs,{__Function}],
	obs_/;MatchQ[obs,{__Function}],
	NN_Integer,UnlabeledVerticesYes_Integer,minEtoBeat_Real]:=

(*************************************)
(***  Last updated on: 03/05/2026  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the weep is done with labeled or unlabeled graphs,
*3. 02/16/2026 Update - made code more efficient by exporting the computeMBF and cleaning up the swap functions*)
(*
Implements Parallel tempering algorithm on graph models at a fixed inverse temperature beta and at several different values of external field parameter. It does not determine the external field schedule here, it assumes the input external field list comes from a scheduler such as  determined by Constant Entropy Increase (CEI) or Constant Thermodynamic Length (CTL) temperature schedule. It equilibriates, computes correlation time, and then applies measurements. It does temperature swaps during the measurement step. It is parallelized so that equilibriation, computation of correlation time, and sweeps during measurement are all done in parallel. Temperature swaps are done on the master kernel. ,

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

Outputs a list with four objects:, 
1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. the final state of each replicas ,
4. an association of each replica and its beta history ,
*)

Module[{result,groundStates,histories,maxGStateCount=500,replicas,numRep,
	minStates,candminE,measurements,numsweeps,Tempoutput,
	ChooseRandomIndependentEdgeSet,Swap,printCase,repNumSweeps,chart},

groundStates=<||>;
histories=<||>;
numRep=Length[replicaSwapEdgesLabels[[2]]];

PrintTemporary["Running PT for graph with numRep ",numRep," at beta ",bt,
	" replicaSwapEdgesLabels ",replicaSwapEdgesLabels];

(*(******************************)
(**      Equilibriate          **)
(******************************)*)

Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphEquilibriate[seedGraph,bt,
		Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
		UnlabeledVerticesYes]
	,{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];

(*Prepare replicas*)
replicas=Tempoutput[[All,2]];

groundStates=Tempoutput[[All,1]];

(*(******************************)
(*  Compute Correlation times *)
(******************************)*)

Tempoutput=Association[
	ParallelTable[i->ECGrav`GraphComputeCorrelationTime[
		replicas[[Key[i],Key["state"],Key["graph"]]],bt,
		Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
		replicas[[Key[i],Key["eqlT"]]],
		groundStates[[Key[i],Key["minEnergy"]]],conjugateObs,UnlabeledVerticesYes]
		,{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
];


(*Update replicas*)
Do[replicas[[Key[i]]]=Tempoutput[[Key[i],2]],{i,numRep}];

(* Extract the minimum energy and corresponding states *)
Do[
	candminE=Tempoutput[[Key[i],1,"minEnergy"]];


	If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
		groundStates[[Key[i],Key["minEnergy"]]]=candminE;
		groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
		If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
			groundStates[[Key[i],Key["minEstates"]]]=Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
		];
	]

,{i,numRep}];


histories=<|Table[
			i-><|"externalField"->replicas[[Key[i],Key["externalField"]]],
			"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>
			,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["externalField"]]],{i,numRep}];

(*(****************************)
(*  Define the replica swap *)
(****************************)*)

(*Helper function to choose a set of non-intersecting edges*)
ChooseRandomIndependentEdgeSet[edgeList_List]:=
Module[{output={},remaining=edgeList,randomEdge},
	While[remaining!={},
		randomEdge=RandomChoice[remaining];
		output=Join[output,{randomEdge}];
		remaining=Select[remaining,!(MemberQ[#,randomEdge[[1]]]||MemberQ[#,randomEdge[[2]]])&]
	];
	output
];

Swap[]:=(*attempts swaps between replicas according to weight exp(-beta(delta H)(delta E)). *)
Module[{replicaSwapMatchings,thisReplInd,nextReplInd,betaDelHDelConj,accept,
	tempThisReplicaState},


	replicaSwapMatchings=ChooseRandomIndependentEdgeSet[replicaSwapEdgesLabels[[1]]];


	Do[
		thisReplInd=edgeIndex[[1]];
		nextReplInd=edgeIndex[[2]];

		betaDelHDelConj=bt*(replicas[[Key[thisReplInd],Key["externalField"]]]
							-replicas[[Key[nextReplInd],Key["externalField"]]])
							. (Through[conjugateObs[replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]]]
								-Through[conjugateObs[replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]]]
							);

		(*Increase the number of swap tries by one for both replicas*)

		histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;
		histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;

		accept=0;
		If[betaDelHDelConj>0,accept = 1,
			If[RandomReal[]<Exp[betaDelHDelConj],accept =1];
		];


		If[accept==1, (*Do swap of replicas*)

			(*Update states: store the state of the high temp replica temporarily*)

			tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
			replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
			replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

			(*Update the energies to reflect the new external fiel dvalues*)
			replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]},replicas[[Key[thisReplInd],Key["externalField"]]]]];
			replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]},replicas[[Key[nextReplInd],Key["externalField"]]]]];

			(*Update histories*)

			histories[[Key[thisReplInd],Key["history"],1]]
				=histories[[Key[nextReplInd],Key["history"],-1]];
			histories[[Key[nextReplInd],Key["history"],1]]
				=histories[[Key[thisReplInd],Key["history"],-1]];

			histories[[Key[thisReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
			histories[[Key[nextReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

			(*Increase the number of swap accepts by one for both replicas*)
			histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
			histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

		];


	,{edgeIndex,replicaSwapMatchings}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)

numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
	While[numsweeps<=NN,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		Tempoutput=Association[
			ParallelTable[
				<|i->ECGrav`GraphSweepReplica[
						replicas[[Key[i],Key["state"],Key["graph"]]],bt,
						Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
						replicas[[Key[i],Key["corrT"]]],
						groundStates[[Key[i],Key["minEnergy"]]],
						UnlabeledVerticesYes]
				|>,

			{i,numRep},
			DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]
		];

		(*Update replicas*)
		Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i],2]],{i,numRep}];

		(*Update minimum energy*)
		Do[
			candminE=Tempoutput[[Key[i],1,"minEnergy"]];

			If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
				groundStates[[Key[i],Key["minEnergy"]]]=candminE;
				groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
				If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
					groundStates[[Key[i],Key["minEstates"]]]
						=Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
				];
			];

		,{i,numRep}];

		Swap[];

		If[Mod[numsweeps,1]==0,
			Table[
				Sow[{numsweeps,replicaSwapEdgesLabels[[2,Key[i]]],
					replicas[[Key[i],Key["state"],Key["energy"]]],
					Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					replicas[[Key[i],Key["state"],Key["energy"]]]
						-replicaSwapEdgesLabels[[2,Key[i]]]
							. Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					Through[obs[replicas[[Key[i],Key["state"],Key["graph"]]]]]},i]
			,{i,numRep}];
		];

		numsweeps++;

	]
][[2]];


chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;numRep]]}]|>;

(*remove the -1's in the initiation*)
Do[If[histories[[Key[i],Key["swapAccept"]]]<Length[histories[[Key[i],Key["history"]]]],
	histories[[Key[i],Key["history"]]]
		=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]])]
,{i,numRep}];

result={ groundStates,chart,replicas,histories};

Remove[groundStates,maxGStateCount,replicas,histories,numRep,minStates,candminE,measurements,numsweeps,Tempoutput,ChooseRandomIndependentEdgeSet,Swap,printCase,repNumSweeps,chart];

result

];


(* ::Item::Closed:: *)
(*Overload GraphParallelTempering with one or more external fields with input of previous external field replica as input*)


ECGrav`GraphParallelTempering[
	inputReplicas_Association/;(Length[DeleteDuplicates[inputReplicas[[All,"beta"]]]]==1&&Length[DeleteDuplicates[inputReplicas[[All,"externalField"]]]]>1),
	hamiltonian_[hparams___],delH_[delHparams___],replicaSwapEdgesLabels_List,
	conjugateObs_/;MatchQ[conjugateObs,{__Function}],
	obs_/;MatchQ[obs,{__Function}],NN_Integer,UnlabeledVerticesYes_Integer,
	minEtoBeat_Real]/;(Length[inputReplicas[[1,Key["externalField"]]]]==Length[replicaSwapEdgesLabels[[2,1]]]==Length[conjugateObs]):=

(*************************************)
(***  Last updated on: 03/03/2026  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the sweep is done with labeled or unlabeled graphs,
*3. 02/16/2026 Update - made code more efficient by exporting the computeMBF and cleaning up the swap functions*)
(*
This overload takes as input an association that is already equilibriated and with correlation times determined and takes measurements applies Parallel tempering algorithm. It sweeps, swaps, and takes measurements. ,

Depends on the functions: GraphSweepReplicas, CorrelationTime, GraphEquilibriate, GraphComputeCorrelationTime., 

Inputs are:, 
1. inputReplicas - association of replicas with ,
2. minEtoBeat -  seed minimum energy at or below which to save graphs,
3. hamiltonian =  function that assigns graphs energy,
4. delH = function that gives the change in energy when a single edge is flipped,
5. replicaSwapEdgesLabels = A list with two elements, the replica swap graph edges 
	and an association mapping graph vertices to external field values 
6. conjugateObs = a list of functions corresponding to the conjugate variables to the external fields
7. obs = the observable quantity in question,
8. NN = number of independent sweeps to be carried out(so that actual number of sweeps is correlation time times NN),
9. UnlabeledVerticesYes = 0 means no selection probability to make the graphs unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.,

Outputs a list with four objects:, 1. the lowest energy and corresponding states found from running the whole program, 
2. an association of each temperature with a table of energies, magnetizations, and meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = (NN * corrT times)/numberOfDataPoints.,
3. the final state of each replicas ,
4. an association of each replica and its beta history ,

*)

Module[{result,numRep=Length[inputReplicas],bt=inputReplicas[[1,"beta"]],
	replicas=inputReplicas,groundStates,histories=<||>,
	maxGStateCount=500,minStates,candminE,measurements,numsweeps,Tempoutput,
	ChooseRandomIndependentEdgeSet,Swap,printCase,repNumSweeps,chart},


groundStates=
	<|Table[
		candminE=Apply[hamiltonian,Join[{replicas[[Key[i],Key["state"],Key["graph"]]]},
						replicas[[Key[i],Key["externalField"]]]]];
		i->If[candminE<= minEtoBeat,
				<|"minEnergy"->candminE,"minEstates"->{replicas[[Key[i],Key["state"],Key["graph"]]]}|>,
				<|"minEnergy"->minEtoBeat,"minEstates"->{}|>
			],
{i,Keys[replicas]}]|>;

histories=<|Table[i-><|"externalField"->replicas[[Key[i],Key["externalField"]]],
			"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["externalField"]]],
	{i,numRep}];

(*(****************************)
(*  Define the replica swap *)
(****************************)*)

(*Helper function to choose a set of non-intersecting edges*)
ChooseRandomIndependentEdgeSet[edgeList_List]:=
Module[{output={},remaining=edgeList,randomEdge},
	While[remaining!={},
		randomEdge=RandomChoice[remaining];
		output=Join[output,{randomEdge}];
		remaining=Select[remaining,!(MemberQ[#,randomEdge[[1]]]||MemberQ[#,randomEdge[[2]]])&]
	];
output
];

Swap[]:=(*attempts swaps between replicas according to weight exp(-beta(delta H)(delta E)). *)
Module[{replicaSwapMatchings,thisReplInd,nextReplInd,betaDelHDelConj,accept,
	tempThisReplicaState},

	replicaSwapMatchings=ChooseRandomIndependentEdgeSet[replicaSwapEdgesLabels[[1]]];
	
	Do[
		thisReplInd=edgeIndex[[1]];
		nextReplInd=edgeIndex[[2]];
		
		betaDelHDelConj=
			bt*(replicas[[Key[thisReplInd],Key["externalField"]]]
				-replicas[[Key[nextReplInd],Key["externalField"]]])
				. (Through[conjugateObs[replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]]]
					-Through[conjugateObs[replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]]]
				);

		(*Increase the number of swap tries by one for both replicas*)

		histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;
		histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;

		accept=0;
		If[betaDelHDelConj>0,accept = 1,
			If[RandomReal[]<Exp[betaDelHDelConj],accept =1];
		];


		If[accept==1, (*Accept swap of replicas*)

			(*Update states: store the state of the high temp replica temporarily*)

			tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
			replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
			replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

			(*Update the energies to reflect the new external fiel dvalues*)
			replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]},replicas[[Key[thisReplInd],Key["externalField"]]]]];
			replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]},replicas[[Key[nextReplInd],Key["externalField"]]]]];

			(*Update histories*)

			histories[[Key[thisReplInd],Key["history"],1]]
				=histories[[Key[nextReplInd],Key["history"],-1]];
			histories[[Key[nextReplInd],Key["history"],1]]
				=histories[[Key[thisReplInd],Key["history"],-1]];

			histories[[Key[thisReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
			histories[[Key[nextReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

			(*Increase the number of swap accepts by one for both replicas*)
			histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
			histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

		];

	,{edgeIndex,replicaSwapMatchings}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
	While[numsweeps<=NN,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		Tempoutput=Association[
			ParallelTable[
				<|i->ECGrav`GraphSweepReplica[replicas[[Key[i],Key["state"],Key["graph"]]],
					bt,
					Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
					Apply[delH,replicaSwapEdgesLabels[[2,Key[i]]]],
					replicas[[Key[i],Key["corrT"]]],
					groundStates[[Key[i],Key["minEnergy"]]],
					UnlabeledVerticesYes]
				|>,
			{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]

		];


		(*Update replicas*)
		Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i],2]],{i,numRep}];

		(*Update minimum energy*)
		Do[
			candminE=Tempoutput[[Key[i],1,"minEnergy"]];

			If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
				groundStates[[Key[i],Key["minEnergy"]]]=candminE;
				groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
				If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
					groundStates[[Key[i],Key["minEstates"]]]=Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
				];
			];
		,{i,numRep}];

		Swap[];

		If[Mod[numsweeps,1]==0,
			Table[
				Sow[{numsweeps,
					replicaSwapEdgesLabels[[2,Key[i]]],
					replicas[[Key[i],Key["state"],Key["energy"]]],
					Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					replicas[[Key[i],Key["state"],Key["energy"]]]
						-replicaSwapEdgesLabels[[2,Key[i]]] . Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					Through[obs[replicas[[Key[i],Key["state"],Key["graph"]]]]]},i]
			,{i,numRep}];
		];

		numsweeps++;
	]
][[2]];

chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;numRep]]}]|>;


Do[
	If[histories[[Key[i],Key["swapAccept"]]]<Length[histories[[Key[i],Key["history"]]]],
	histories[[Key[i],Key["history"]]]
		=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]])
	],
{i,numRep}];(*remove the -1's in the initiation*)

result={ groundStates,chart,replicas,histories};

Remove[numRep,bt,replicas,groundStates,histories,maxGStateCount,minStates,candminE,
	measurements,numsweeps,Tempoutput,ChooseRandomIndependentEdgeSet,Swap,printCase,
	repNumSweeps,chart];

result

];


(*Overload with no delH *)
ECGrav`GraphParallelTempering[
	inputReplicas_Association/;(Length[DeleteDuplicates[inputReplicas[[All,"beta"]]]]==1&&Length[DeleteDuplicates[inputReplicas[[All,"externalField"]]]]>1),
	hamiltonian_[hparams___],replicaSwapEdgesLabels_List,
	conjugateObs_/;MatchQ[conjugateObs,{__Function}],
	obs_/;MatchQ[obs,{__Function}],NN_Integer,UnlabeledVerticesYes_Integer,
	minEtoBeat_Real]/;(Length[inputReplicas[[1,Key["externalField"]]]]==Length[replicaSwapEdgesLabels[[2,1]]]==Length[conjugateObs]):=

(*************************************)
(***  Last updated on: 03/03/2026  ***)
(*************************************)
(*Notes: ,
*1. 01/21/2023 Update - Memory leak in the parallelization is fixed.,
*2. 1/30/2025 Update - enabeled one to chose whether the sweep is done with labeled or unlabeled graphs,
*3. 02/16/2026 Update - made code more efficient by exporting the computeMBF and cleaning up the swap functions*)
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
3. the final state of each replicas ,
4. an association of each replica and its beta history ,

*)

Module[{result,numRep=Length[inputReplicas],bt=inputReplicas[[1,"beta"]],
	replicas=inputReplicas,groundStates,histories=<||>,
	maxGStateCount=500,minStates,candminE,measurements,numsweeps,Tempoutput,
	ChooseRandomIndependentEdgeSet,Swap,printCase,repNumSweeps,chart},


groundStates=
	<|Table[
		candminE=Apply[hamiltonian,Join[{replicas[[Key[i],Key["state"],Key["graph"]]]},
						replicas[[Key[i],Key["externalField"]]]]];
		i->If[candminE<= minEtoBeat,
				<|"minEnergy"->candminE,"minEstates"->{replicas[[Key[i],Key["state"],Key["graph"]]]}|>,
				<|"minEnergy"->minEtoBeat,"minEstates"->{}|>
			],
{i,Keys[replicas]}]|>;

histories=<|Table[i-><|"externalField"->replicas[[Key[i],Key["externalField"]]],
			"swapAccept"->0,"swapTry"->0, "history"->Table[-1.0,{NN}]|>,{i,numRep}]|>;

Do[histories[[Key[i],Key["history"],-1]]=replicas[[Key[i],Key["externalField"]]],
	{i,numRep}];

(*(****************************)
(*  Define the replica swap *)
(****************************)*)

(*Helper function to choose a set of non-intersecting edges*)
ChooseRandomIndependentEdgeSet[edgeList_List]:=
Module[{output={},remaining=edgeList,randomEdge},
	While[remaining!={},
		randomEdge=RandomChoice[remaining];
		output=Join[output,{randomEdge}];
		remaining=Select[remaining,!(MemberQ[#,randomEdge[[1]]]||MemberQ[#,randomEdge[[2]]])&]
	];
output
];

Swap[]:=(*attempts swaps between replicas according to weight exp(-beta(delta H)(delta E)). *)
Module[{replicaSwapMatchings,thisReplInd,nextReplInd,betaDelHDelConj,accept,
	tempThisReplicaState},

	replicaSwapMatchings=ChooseRandomIndependentEdgeSet[replicaSwapEdgesLabels[[1]]];
	
	Do[
		thisReplInd=edgeIndex[[1]];
		nextReplInd=edgeIndex[[2]];
		
		betaDelHDelConj=
			bt*(replicas[[Key[thisReplInd],Key["externalField"]]]
				-replicas[[Key[nextReplInd],Key["externalField"]]])
				. (Through[conjugateObs[replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]]]
					-Through[conjugateObs[replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]]]
				);

		(*Increase the number of swap tries by one for both replicas*)

		histories[[Key[nextReplInd],Key["swapTry"]]]+=1.0;
		histories[[Key[thisReplInd],Key["swapTry"]]]+=1.0;

		accept=0;
		If[betaDelHDelConj>0,accept = 1,
			If[RandomReal[]<Exp[betaDelHDelConj],accept =1];
		];


		If[accept==1, (*Accept swap of replicas*)

			(*Update states: store the state of the high temp replica temporarily*)

			tempThisReplicaState=replicas[[Key[thisReplInd],Key["state"]]];
			replicas[[Key[thisReplInd],Key["state"]]]=replicas[[Key[nextReplInd],Key["state"]]];
			replicas[[Key[nextReplInd],Key["state"]]]=tempThisReplicaState;

			(*Update the energies to reflect the new external fiel dvalues*)
			replicas[[Key[thisReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[thisReplInd],Key["state"],Key["graph"]]]},replicas[[Key[thisReplInd],Key["externalField"]]]]];
			replicas[[Key[nextReplInd],Key["state"],Key["energy"]]]
				=Apply[hamiltonian,Join[{replicas[[Key[nextReplInd],Key["state"],Key["graph"]]]},replicas[[Key[nextReplInd],Key["externalField"]]]]];

			(*Update histories*)

			histories[[Key[thisReplInd],Key["history"],1]]
				=histories[[Key[nextReplInd],Key["history"],-1]];
			histories[[Key[nextReplInd],Key["history"],1]]
				=histories[[Key[thisReplInd],Key["history"],-1]];

			histories[[Key[thisReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[thisReplInd],Key["history"]]]];
			histories[[Key[nextReplInd],Key["history"]]]
				=RotateLeft[histories[[Key[nextReplInd],Key["history"]]]];

			(*Increase the number of swap accepts by one for both replicas*)
			histories[[Key[thisReplInd],Key["swapAccept"]]]+=1.0;
			histories[[Key[nextReplInd],Key["swapAccept"]]]+=1.0 ;

		];

	,{edgeIndex,replicaSwapMatchings}];

];


(*(******************************************
* Start taking measurements and swapping **
******************************************)*)


numsweeps=1;

printCase=Floor[(NN*1.0)/5.0];


measurements=Reap[
	While[numsweeps<=NN,

		If[Mod[numsweeps,printCase]==0,PrintTemporary[" sweepno ",numsweeps]];

		Tempoutput=Association[
			ParallelTable[
				<|i->ECGrav`GraphSweepReplica[replicas[[Key[i],Key["state"],Key["graph"]]],
					bt,
					Apply[hamiltonian,replicaSwapEdgesLabels[[2,Key[i]]]],
					replicas[[Key[i],Key["corrT"]]],
					groundStates[[Key[i],Key["minEnergy"]]],
					UnlabeledVerticesYes]
				|>,
			{i,numRep},DistributedContexts->{$Context, "ECGrav`MCSims`Private`"}]

		];


		(*Update replicas*)
		Do[replicas[[Key[i],Key["state"]]]=Tempoutput[[Key[i],2]],{i,numRep}];

		(*Update minimum energy*)
		Do[
			candminE=Tempoutput[[Key[i],1,"minEnergy"]];

			If[candminE<groundStates[[Key[i],Key["minEnergy"]]],
				groundStates[[Key[i],Key["minEnergy"]]]=candminE;
				groundStates[[Key[i],Key["minEstates"]]]=Tempoutput[[Key[i],1,"minEstates"]],
				If[candminE==groundStates[[Key[i],Key["minEnergy"]]]&&Length[groundStates[[Key[i],Key["minEstates"]]]]<=maxGStateCount,
					groundStates[[Key[i],Key["minEstates"]]]=Union[groundStates[[Key[i],Key["minEstates"]]],Tempoutput[[Key[i],1,"minEstates"]]];
				];
			];
		,{i,numRep}];

		Swap[];

		If[Mod[numsweeps,1]==0,
			Table[
				Sow[{numsweeps,
					replicaSwapEdgesLabels[[2,Key[i]]],
					replicas[[Key[i],Key["state"],Key["energy"]]],
					Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					replicas[[Key[i],Key["state"],Key["energy"]]]
						-replicaSwapEdgesLabels[[2,Key[i]]] . Through[conjugateObs[replicas[[Key[i],Key["state"],Key["graph"]]]]],
					Through[obs[replicas[[Key[i],Key["state"],Key["graph"]]]]]},i]
			,{i,numRep}];
		];

		numsweeps++;
	]
][[2]];

chart=<|Table[i[[1,2]]->i,{i,measurements[[1;;numRep]]}]|>;


Do[
	If[histories[[Key[i],Key["swapAccept"]]]<Length[histories[[Key[i],Key["history"]]]],
	histories[[Key[i],Key["history"]]]
		=(histories[[Key[i],Key["history"],-(histories[[Key[i],Key["swapAccept"]]]+1);;-1]])
	],
{i,numRep}];(*remove the -1's in the initiation*)

result={ groundStates,chart,replicas,histories};

Remove[numRep,bt,replicas,groundStates,histories,maxGStateCount,minStates,candminE,
	measurements,numsweeps,Tempoutput,ChooseRandomIndependentEdgeSet,Swap,printCase,
	repNumSweeps,chart];

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
