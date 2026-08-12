(* ::Package:: *)

(* ::Input:: *)
(*(*:Name: ECGrav - Emergent Combinatorial Gravity*)*)
(*(*:Author: Kassahun Betre*)*)
(*(*:Contributors: Khai Luong*)*)
(*(*:Date: 08/06/2025*)*)
(*(*: To install package:,*)
(*   1. Save as ECGrav.wl, *)
(*   2. Evaluate Notebook, *)
(*   3. File->Install->Choose: Type = Package, Source=From File, Select ECGrav.wl,     Instal Name = ECGrav  *)*)
(*(*: To load package evaluate Needs["ECGrav`"] in a new notebook*)*)
(**)


(* ::Title:: *)
(*Begin Main ECGrav Package*)


BeginPackage["ECGrav`"];


Unprotect @@ Names["ECGrav`*"];

(* Everything except the user-settable package settings below. Re-loading the package in a
	live session -- which is also what happens when a parallel subkernel loads its own copy --
	must not silently discard a setting the session has already made. *)
(* Names["ECGrav`*"] yields SHORT names, so the exclusion must be written without the context
	prefix or it silently matches nothing and ClearAll wipes the setting anyway. *)
ClearAll @@ DeleteCases[Names["ECGrav`*"], "$ECGravMaxEquilibriationSweeps"];


(* ::Title:: *)
(*ECGrav Package Settings*)


(* :Usage Messages: *)

(* :Usage Messages: *)

ECGrav::usage="ECGrav carries the package-wide diagnostics that are shared rather than
	duplicated on each function: see ECGrav::badham and ECGrav::baddelh, issued by the Monte
	Carlo drivers when an energy callback cannot be called the way this package calls it.";

(* :Error Messages: *)

(* The Monte Carlo functions take their energy callbacks as inert head[params] expressions,
	bound by hamiltonian_[hparams___] and delH_[delHparams___]. Those patterns match ANY
	expression that has a head, so a misspelt symbol, a curried pure Function, or a wrapper of
	the wrong shape all pass argument checking unrejected and only surface much later as a Part
	or Function error naming nothing involved. No pattern can separate them -- a misspelt symbol
	is structurally identical to a correct one -- so the callbacks are checked by evaluating
	them once and requiring a number back. *)

ECGrav::badham=
	"The hamiltonian argument `1` is not usable: hamiltonian[graph, params] evaluated to `2` "<>
	"rather than a number. Pass an inert head[params] whose head is defined over "<>
	"(adjacencyMatrix, params...), for example h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed "<>
	"as h[-1.,0.]. A head with no definition, or a curried pure Function such as "<>
	"HIsing[#,-1.,0.]&, matches the argument pattern but cannot be called this way.";

ECGrav::baddelh=
	"The delH argument `1` is not usable: delH[graph, params, i, j] evaluated to `2` rather "<>
	"than a number. Pass an inert head[params] whose head is defined over "<>
	"(adjacencyMatrix, params..., i, j), for example "<>
	"dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.]. "<>
	"A curried pure Function such as delHIsing[#1,-1.,0.,#2,#3]& matches the argument pattern "<>
	"but binds delH to Function itself, and the four-argument call then fails.";


$ECGravMaxEquilibriationSweeps::usage="$ECGravMaxEquilibriationSweeps is the sweep budget every
	equilibriation routine works within: GraphEquilibriate and
	RandomPureSimplicialComplexMCMCEquilibriate each give up once they have run this many
	sweeps without meeting their convergence criterion. Exhausting the budget is reported
	through the corresponding ::noconv message, and the routine then returns
	\"converged\" -> False together with an eqlT that is a lower bound rather than an
	estimate. Raise it for systems that mix slowly. Setting it by plain assignment is enough:
	the parallel drivers push the current value out to the subkernels before launching work
	that can reach an equilibriator, so a subkernel's own copy of the package does not
	override the setting made in the master session.";

(* Guarded rather than assigned outright: a subkernel that already received the master's value,
	or a session that set it before re-loading, keeps its setting instead of being reset to the
	default. Paired with the DeleteCases above, which keeps ClearAll off this symbol. *)
If[!ValueQ[$ECGravMaxEquilibriationSweeps],
	$ECGravMaxEquilibriationSweeps=25000];


(* ::Title::Closed:: *)
(*ECGrav Public Functions Involving MCSims*)


(* ::Chapter:: *)
(*Helper Functions*)


(* ::Section::Closed:: *)
(*Aggregating Data*)


(* ::Item::Closed:: *)
(*CorrelationTime*)


(* :Usage Mesages: *)

CorrelationTime::usage="CorrelationTime[t,tbl] computes the autocorrelation at 
time t of the list of values given by tbl.
Inputs are:
1. t= Integer,\[IndentingNewLine]2. tbl = List of observable values such as magnetization or energy,\[IndentingNewLine]It returns a real number";
(* :Error Mesages: *)

CorrelationTime::argerr="An integer is expected at position 1, and a list of 
numbers at position 2.";


(* ::Item::Closed:: *)
(*EmpCorrelationTime*)


(* :Usage Mesages: *)

EmpCorrelationTime::usage="EmpCorrelationTime[t,tbl] computes the empirical autocorrelation
at time t of the list of values given by tbl based on the lecture by Evertz, 2020.
Inputs are:
1. t= Integer,\[IndentingNewLine]2. tbl = List of observable values such as magnetization or energy,\[IndentingNewLine]It returns a real number";
(* :Error Mesages: *)

EmpCorrelationTime::argerr="An integer is expected at position 1, and a list of 
numbers at position 2.";


(* ::Item::Closed:: *)
(*ErrorBootstrap*)


(* :Usage Mesages: *)

ErrorBootstrap::usage="ErrorBootstrap[formula,data] gives the uncertainty in the 
value of formula[data] for all i in data computed using the 
bootstrap or resampling method (Section 3.4.3 of Newman & Barkema).
Inputs are:
1. formula = formula that gives a number when applied on the data, e.g. Mean[data],\[IndentingNewLine]2. data = List of values,\[IndentingNewLine]It returns a real number";
(* :Error Mesages: *)

ErrorBootstrap::argerr="A formula is expected at position 1, and a list of 
numbers at position 2.";


(* ::Item::Closed:: *)
(*SpecificHeat*)


(* :Usage Mesages: *)

SpecificHeat::usage="SpecificHeat[energyTable,NN,beta] gives the specific heat 
per site where for the list of energy measurements energyTable given NN 
number of sites, and inverse temperature beta).
Inputs are:
1. energyTable = a list of energy values, \[IndentingNewLine]2. NN = Integer is teh number of sites,
3. beta = inverse temperature\[IndentingNewLine]It returns a real number equal to the specific heat per site.";
(* :Error Mesages: *)

SpecificHeat::argerr="A list of reals is expected at position 1, an integer at position
2, and a real number at position 3.";


(* ::Item::Closed:: *)
(*DSpecificHeat*)


(* :Usage Mesages: *)

DSpecificHeat::usage="DSpecificHeat[energyTable,NN,beta] gives the derivative of the 
specific heat per site wrt inverse temperature where for the list of energy measurements 
energyTable given NN number of sites, and inverse temperature beta).
Inputs are:
1. energyTable = a list of energy values, \[IndentingNewLine]2. NN = Integer is teh number of sites,
3. beta = inverse temperature\[IndentingNewLine]It returns a real number equal to the D[specific heat per site,beta].";
(* :Error Mesages: *)

DSpecificHeat::argerr="A list of reals is expected at position 1, an integer at position
2, and a real number at position 3.";


(* ::Item::Closed:: *)
(*Susceptibility*)


(* :Usage Mesages: *)

Susceptibility::usage="Susceptibility[obsTable,NN,beta] gives the susceptibility of 
the observable (e.g. magnetization) per site for the list of observable values 
obsTable, number of sites NN, and inverse temperature beta).
Inputs are:
1. obsTable = List, a list of observable values, \[IndentingNewLine]2. NN = Integer, the number of sites,
3. beta = Real, inverse temperature\[IndentingNewLine]It returns a real number equal to the susceptibility per site.";
(* :Error Mesages: *)

Susceptibility::argerr="A list of reals is expected at position 1, an integer at 
position 2, and a real number at position 3.";


(* ::Item::Closed:: *)
(*LogSumExp*)


(* :Usage Mesages: *)

LogSumExp::usage="LogSumExp[lst] Computes the log of the sum of exponentials of the 
values in the list lst.
Inputs are:
1. lst = List, a list of values, \[IndentingNewLine]It returns a real number equal to the Log[Sum[Exp[i],{i,lst}]].";
(* :Error Mesages: *)

LogSumExp::argerr="A list of reals is expected at position 1.";


(* ::Item::Closed:: *)
(*ComputeMinusBetaTimesFreeEnergy*)


(* :Usage Mesages: *)

ComputeMinusBetaTimesFreeEnergy::usage="ComputeMinusBetaTimesFreeEnergy[dat_Association] 
Computes the value of -beta*free energy at the values of the inverse temperature beta
given as the keys in the association dat.
Overloads: 
	ComputeMinusBetaTimesFreeEnergy[dat_Association, bt_Real] computes 
	-beta*free energy at the values of the external field values
	given as the keys in the association dat at the inverse temperature bt.
	
Inputs is:,
1. dat - an association of external field values (J) and the corresponding 
   list of values of the conjugate field (O) measured at these external field values.
   (i.e., the Hamiltonian is assumed to have the term +J*O)
   e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
2. bt - Real = the inverse temperature

Returns an association with the external field values (J_i) as keys and 
-(betaFixed)(F[betaFixed,J_i] as the values.
";

(* :Error Mesages: *)

ComputeMinusBetaTimesFreeEnergy::argerr="Input has to be of the form 
	ComputeMinusBetaTimesFreeEnergy[dat_Association]
	or 
	ComputeMinusBetaTimesFreeEnergy[dat_Association, bt_Real]";


(* ::Item::Closed:: *)
(*NegativeBetaTimesFreeEnergy*)


(* :Usage Mesages: *)

NegativeBetaTimesFreeEnergy::usage="NegativeBetaTimesFreeEnergy[bf_Real,minusBetaF_Association,energyMeasurements_Association] 
Computes the value of -beta*(free energy) at the value of the inverse temperature bf.
Overloads: 
	NegativeBetaTimesFreeEnergy[betaFixed_Real,targetExtField_Real,
		minusBetaF_Association,conjugateExtFieldMeasurements_Association] computes 
	-beta*free energy at the value of the inverse temperature betaFixed and target 
	external field targetExtField.
	
Inputs are:
1. bf = Real, inverse temperature
2. targetExtField - Real value for the external field at which the extrapolated value 
	of -beta*freeEnergy is wanted,
3. minusBetaF = Association, an association of inverse temperatures and the 
   corresponding value of -beta*(free energy) computed from the simulation at
   each beta, e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,
4. energyMeasurements (or conjugateExtFieldMeasurements) = an association of inverse 
	temperature (or external field )
   and the corresponding list of energies (or conjugate observable values) measured at 
	that beta (or external value)., 
   e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta (external field) values which are the keys for both associations have 
to be equal as sets! Also, the lengths of the lists of energy measurements have to be 
equal for all betas/external fields.\[IndentingNewLine]It returns a real number equal to the -(bf)F[bf], or -(betaFixed)F[betaFixed,targetExtField]";
(* :Error Mesages: *)

NegativeBetaTimesFreeEnergy::argerr="Input has to be of the form 
	NegativeBetaTimesFreeEnergy[bf_Real,minusBetaF_Association,energyMeasurements_Association]
	or 
	NegativeBetaTimesFreeEnergy[betaFixed_Real,targetExtField_Real,
		minusBetaF_Association,conjugateExtFieldMeasurements_Association]";


(* ::Item::Closed:: *)
(*InternalEnergy*)


(* :Usage Mesages: *)

InternalEnergy::usage="InternalEnergy[bf,minusBetaF,energyMeasurements] 
Computes the value of the interval energy U[bf] at the value of the inverse temperature bf.
Inputs are:
1. bf = Real, inverse temperature
2. minusBetaF = Association, an association of inverse temperatures and the 
   corresponding value of -beta*(free energy) computed from the simulation at
   each beta, e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,
3. energyMeasurements = energyMeasurements - an association of inverse temperature 
   and the corresponding list of energies measured at that beta., 
   e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta values which are the keys for both associations have to be equal as 
sets! Also, the lengths of the lists of energy measurements have to be equal for all 
betas.\[IndentingNewLine]It returns a real number equal to U[bf].";
(* :Error Mesages: *)

InternalEnergy::argerr="A Real is expected at position 1, an association
at position 2, and an association at 3.";


(* ::Item::Closed:: *)
(*ExtrapolatedExpectationValue*)


(* :Usage Mesages: *)

ExtrapolatedExpectationValue::usage="ExtrapolatedExpectationValue[bf_Real,
	minusBetaF_Association,energyMeasurements_Association, 
	measuredObservableValues_Association] 
Computes the extrapolated expectation value of the observable at the target value of 
the inverse temperature bf.
Overloads: 
	ExtrapolatedExpectationValue[betaFixed_Real,targetExtField_Real,
		minusBetaF_Association,conjugateExtFieldMeasurements_Association,
		measuredObservableValues_Association] 
	computes the extrapolated value of the observable O at the fixed value of the 
	inverse temperature betaFixed and at the value of the external field targetExtField. 
	The microscopic Hamiltonian is assumed to have a term of the form J*E where J is tunable 
	external field parameter and E is the conjugate observable..
	
Inputs are:,
1. betaFixed - Real, the fixed value of the  inverse temperature,\[IndentingNewLine]2. targetExtField - Real, the value of the target external field,\[IndentingNewLine]3. minusBetaF - an association of the external field H and the corresponding value of -beta*free energy(H) e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,\[IndentingNewLine]4. conjugateExtFieldMeasurements - an association of the external field values J and the corresponding list of energies (list of E's without the J multiplying them)measured at that J. e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>,\[IndentingNewLine]5. measuredObservableValues - an association of the external field values J and the corresponding list of meaurements of the observable O measured at those external field values. e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>\[IndentingNewLine]Note, the J values which are the keys for all the associations have to be equal as sets! Also, the lengths of the lists of values have to be equal for all J values,\[IndentingNewLine]
Returns a real number equal to  the extrapolated expectation value of the observable at 
the target value of the inverse temperature bf (or the target external field parameter 
and fixed temperature).
";

(* :Error Mesages: *)

ExtrapolatedExpectationValue::argerr="Input has to be of the form 
	ExtrapolatedExpectationValue[bf_Real,
	minusBetaF_Association, energyMeasurements_Association, measuredObservableValues_Association] 
	or 
	ExtrapolatedExpectationValue[betaFixed_Real,targetExtField_Real,
		minusBetaF_Association,conjugateExtFieldMeasurements_Association,
		measuredObservableValues_Association]";


(* ::Item::Closed:: *)
(*CvOverT*)


(* :Usage Mesages: *)

CvOverT::usage="CvOverT[bf,minusBetaF,energyMeasurements] 
Computes the value of the specific heat divided by temperature (Cv/T) at the value of 
the inverse temperature bf. Note, it is not computing Cv/T per site!
Inputs are:
1. bf = Real, inverse temperature
2. minusBetaF = Association, an association of inverse temperatures and the 
   corresponding value of -beta*(free energy) computed from the simulation at
   each beta, e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,
3. energyMeasurements = energyMeasurements - an association of inverse temperature 
   and the corresponding list of energies measured at that beta., 
   e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>
Note, the beta values which are the keys for both associations have to be equal as 
sets! Also, the lengths of the lists of energy measurements have to be equal for all 
betas.\[IndentingNewLine]It returns a real number equal to Cv/T at bf.";
(* :Error Mesages: *)

CvOverT::argerr="A Real is expected at position 1, an association
at position 2, and an association at 3.";


(* ::Item::Closed:: *)
(*Constrained Probability of the Conjugate Field*)


(* :Usage Mesages: *)

ConstrainedProbConjugateField::usage="ConstrainedProbConjugateField[betaFixed_Real,targetExtField_List,
	minusBetaF_Association,
	conjugateExtFieldMeasurements_Association]
	Computes the probability density of the magnetization or conjugate field
	P(beta,M) which depends implicitely on the external field H to which M is the 
	conjugate field..
	
Inputs are:,
1. betaFixed - Real, the fixed value of the  inverse temperature,\[IndentingNewLine]2. targetExtField - List, a vector of external-field values at which the density is evaluated,\[IndentingNewLine]3. minusBetaF - an association of the external field H and the corresponding value of
	-beta*free energy(H) e.g. <|0.1 -> -2.3, 2.5 -> 34.2|>,\[IndentingNewLine]4. conjugateExtFieldMeasurements - an association of the external field values J and
	the corresponding list of energies (list of E's without the J multiplying them)
	measured at that J. e.g. <|0.1 -> {1.1,2.3,5.2}, 2.5 -> {-2.0,-4.3,-20.1}|>,
\[IndentingNewLine]Note, the J values which are the keys for all the associations have to be equal as 
sets! Also, the lengths of the lists of values have to be equal for all J values,\[IndentingNewLine]
Returns a list containing a SmoothKernelDistribution of the conjugate field (together with plot ranges), which can be plotted.
";

(* :Error Mesages: *)

ConstrainedProbConjugateField::argerr="Input has to be of the form 
	ConstrainedProbConjugateField[betaFixed_Real,targetExtField_List,
	minusBetaF_Association,conjugateExtFieldMeasurements_Association] ";


(* ::Chapter:: *)
(*Hamiltonians*)


(* ::Section::Closed:: *)
(*Graph Hamiltonians*)


(* ::Item::Closed:: *)
(*HIsing*)


(* :Usage Mesages: *)

HIsing::usage="HIsing[Am,J,L] = J/2 sum_{i!=j}A^2_{ij} + L/2 sum_{i!=j}A_{ij}.
Inputs are:
1. Am = List, adjacency matrix of a graph\[IndentingNewLine]2. J = Real, coupling constant,
3. L = Real, akin to external magnetic field\[IndentingNewLine]It returns a real number";
(* :Error Mesages: *)

HIsing::argerr="An adjacency matrix is expected at position 1, a real number at 
position 2, and a real number at position 3.";


(* :Usage Mesages: *)

delHIsing::usage="delHIsing[Am,J,L,a,b] computes HIsing[Amnew] - HIsing[Am] where 
Amnew if found by toggling Am at 
row a and col b.
Inputs are:
1. Am = List, adjacency matrix of a graph\[IndentingNewLine]2. J = Real, coupling constant,
3. L = Real, akin to external magnetic field,
4. a = Integer, row number
5. b = Integer, column number\[IndentingNewLine]It returns a real number";

(* :Error Mesages: *)
delHIsing::argerr="An adjacency matrix is expected at position 1, a real number at 
position 2, a real number at position 3, an integer at position 4, an integer at 
position 5.";


(* ::Item::Closed:: *)
(*HWeightedFaceCounts*)


(* :Usage Mesages: *)

HWeightedFaceCounts::usage="HWeightedFaceCounts[Am,J1,J2,J3,J4,J5] = 
J1*(# of vertices) + J2*(number of edges)+... + J5*(number of 5-cliques).
Inputs are:
1. Am = List, adjacency matrix of a graph
2. J1 = Real, coupling constant,
3. J2 = Real, coupling constant,
4. J3 = Real, coupling constant,
5. J4 = Real, coupling constant,
6. J5 = Real, coupling constant\[IndentingNewLine]It returns a real number";
(* :Error Mesages: *)

HWeightedFaceCounts::argerr="An adjacency matrix is expected at position 1, and real 
numbers at positions 2 through 6.";


(* ::Item::Closed:: *)
(*HEdgeDeg*)


(* :Usage Mesages: *)

HEdgeDeg::usage="HEdgeDeg[Am,J,D1,D2] = 
J/(2N)*(Tr(A^4-A(2D1K + 2D2I)A)+D1^2*n(n-1)), where K is the 
adjacency matrix of the complete graph.
Inputs are:
1. Am = List, adjacency matrix of a graph\[IndentingNewLine]2. J = Real, coupling constant,
3. D1 = Real, coupling constant,
4. D2 = Real, coupling constant,\[IndentingNewLine]It returns a real number";
(* :Error Mesages: *)

HEdgeDeg::argerr="An adjacency matrix is expected at position 1, and real numbers
at positions 2, 3, 4.";


(* :Usage Mesages: *)

delHEdgeDeg::usage="delHEdgeDeg[Am,J, D1, D2, Amsq,a,b] computes 
HEdgeDeg[Amnew,J,D1,D2]] - HEdgeDeg[Am,J,D1,D2] where 
Amnew is found by toggling Am at row a and col b.
Inputs are:
1. Am = List, adjacency matrix of a graph\[IndentingNewLine]2. J = Real, coupling constant,
3. D1 = Real, akin to external magnetic field,
4. D2 = Integer, row number
5. Amsq = Am.Am, optional input that is not needed because of the overload
6. a = Integer, row number
7. b = Integer, column number\[IndentingNewLine]It returns a real number";

(* :Error Mesages: *)
delHEdgeDeg::argerr="An adjacency matrix is expected at position 1, real numbers at 
positions 2,3,4, optional square symmetric integer matrix equalt to the square
of the adjacency matrix at position 5, an integer at position 6, an integer at 
position 7.";


(* ::Item::Closed:: *)
(*HLaplacian*)


(* :Usage Mesages: *)

HLaplacian::usage="HLaplacian[Amat,J1]=Sum_{ij} deg(i)L_{ij}deg(j), where 
L is the graph Laplacian, L = Diag[deg] - Amat.
Inputs are:
1. Am = List, adjacency matrix of a graph\[IndentingNewLine]2. J = Real, coupling constant,\[IndentingNewLine]It returns a real number";
(* :Error Mesages: *)

HLaplacian::argerr="An adjacency matrix is expected at position 1, a real number at 
position 2.";


(* ::Item::Closed:: *)
(*H1dCombManifold*)


(* :Usage Mesages: *)

H1dCombManifold::usage="H1dCombManifold[Amat,J]. Gives zero on 1D geometric graphs
(cycle and path graphs)
Inputs are:
1. Amat = List, adjacency matrix of a graph
2. J = Real, coupling constant,
It returns a real number";
(* :Error Mesages: *)

H1dCombManifold::argerr="An adjacency matrix is expected at position 1, a real number at 
position 2.";


(* ::Item::Closed:: *)
(*H2dCombManifold*)


(* :Usage Mesages: *)

H2dCombManifold::usage="H2dCombManifold[Amat,J]. Gives zero on 2d combinatorial 
	manifolds (combinatorial manifolds homeomorphic to the sphere S^2 or the disk B^2)
Inputs are:
1. Amat = List, adjacency matrix of a graph
2. J = Real, coupling constant,
It returns a real number";
(* :Error Mesages: *)

H2dCombManifold::argerr="An adjacency matrix is expected at position 1, a real number at 
position 2.";


(* :Usage Mesages: *)

delH2dCombManifold::usage="delH2dCombManifold[Amat_List,J_,a_Integer,b_Integer] computes 
	H2dCombManifold[Amnew] - H2dCombManifold[Am] where Amnew if found by toggling Amat 
	at row a and col b.
Inputs are:
1. Am = List, adjacency matrix of a graph\[IndentingNewLine]2. J = Real, coupling constant,
3. a = Integer, row number
4. b = Integer, column number\[IndentingNewLine]It returns a real number";

(* :Error Mesages: *)
delH2dCombManifold::argerr="Input must be of the form 
	delH2dCombManifold[Amat_List,J_,a_Integer,b_Integer].";


(* ::Item::Closed:: *)
(*H2dPseudoCombManifold*)


(* :Usage Mesages: *)

H2dPseudoCombManifold::usage="H2dPseudoCombManifold[Amat,J]. Gives zero on 2d 
	pseudo-combinatorial manifolds (graphs that 2 pure where the degree of each
	edge is atmost 2 and are 2-path connected)
Inputs are:
1. Amat = List, adjacency matrix of a graph
2. J = Real, coupling constant,
It returns a real number";
(* :Error Mesages: *)

H2dPseudoCombManifold::argerr="An adjacency matrix is expected at position 1, a real number at 
position 2.";


(* ::Chapter:: *)
(*Ground State Searches*)


(* ::Section::Closed:: *)
(*Exact Low Energy State Searches*)


(* ::Item::Closed:: *)
(*LowEnergyStates*)


(* :Usage Mesages: *)

LowEnergyStates::usage="LowEnergyStates[states_List,H_,klevel_Integer] Calculates the 
energies of the states in the list using the hamiltonian function H. Outputs an 
ordered association with the lowest klevel energies as the keys and the states 
with those energies as the values.
Inputs are:
1. states = List, a list of adjacency matrix of a graphs
2. H = function, the Hamiltonian,  
3. klevel = Integer, the number of lowest energy levels to return,

Outputs an association with the lowest klevel energies as keys and lists of states 
with the corresponding energies as the values.";
(* :Error Mesages: *)

LowEnergyStates::argerr="Input should be of the form 
	LowEnergyStates[states_List,H_,klevel_Integer].";


(* ::Section::Closed:: *)
(*Gradient Descent*)


(* ::Item::Closed:: *)
(*GradDescent*)


(* :Usage Mesages: *)

GradDescent::usage="GradDescent[seedAmat,delH[delHparams],cutoff] runs the 
gradient algorithm on the input adjacency matrix seedAmat
Inputs are:
1. seedAmat = List, adjacency matrix of a graph
2. delH = an inert head[params] for delta E when one edge is toggled,
3. cutoff -  number of sweeps,

Outputs an adjacency matrix.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";
(* :Error Mesages: *)

GradDescent::argerr="An adjacency matrix is expected at position 1, a hamiltonian 
at position 2, a formula for delta Hamiltonian in position 3, and integer at position 4.";


(* ::Section::Closed:: *)
(*Stochastic Gradient Descent*)


(* ::Item::Closed:: *)
(*SGradDescent*)


(* :Usage Mesages: *)

SGradDescent::usage="SGradDescent[seedAmat, hamiltonian[hparams],delH[delHparams],cutoff] runs 
the stochastic gradient descent algorithm with softmax at parameter beta on the input 
adjacency matrix seedAmat.
Inputs are:
1. seedAmat = List, adjacency matrix of a graph\[IndentingNewLine]2. hamiltonian = formula for that gives the energy of the graph, coupling constant,
3. delH - a formula for delta E (when one edge is flipped to expedite computation),  
4. beta - inverse temperature,
5. NN -  number of sweeps,\[IndentingNewLine]
Outputs an association with three elements,  
1. the minimum energy visited throughout the search,
2. states with that minimum energy. If multiple states have degenerate minimum energy, 
they will all be included.,
3. The last state visited.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

SGradDescent::argerr="An adjacency matrix is expected at position 1, a hamiltonian 
at position 2, an optional formula for delta Hamiltonian in position 3, a real at 
position 4 (positin 3 if optional delta Hamiltonian formula not given), 
and integer at position 5 (positin 4 if optional delta Hamiltonian formula not given).";


(* ::Section::Closed:: *)
(*Simulated Annealing*)


(* ::Item::Closed:: *)
(*SimulatedAnnealing*)


(* :Usage Mesages: *)

SimulatedAnnealing::usage="SimulatedAnnealing[seedAmat, hamiltonian,delH,betai,betaf,roundLength,NN] 
runs the simulated annealing ground state search algorithm ground state search from initial to final 
inverse temperatures given by betai and betaf on the input adjacency matrix seedAmat.
Inputs are:
1. seedAmat = List, adjacency matrix of a graph\[IndentingNewLine]2. hamiltonian = formula for that gives the energy of the graph, coupling constant,
3. delH - optional formula for delta E (when one edge is flipped to expedite computation),  
4. betai - Real, starting low inverse temperature,
5. betaf - Real, the final high inverse temperature,
6. roundLength - Integer, related to how many rounds of MC steps to do at each 
   temperaturenumber before advancing to the next temperature,
7. NN -  number of sweeps,\[IndentingNewLine]
Outputs an association with five elements,
1. the minimum energy visited throughout the search,
2. the second minimum energy,
3. the 'ground states', the states with that minimum energy. If multiple states have 
   degenerate minimum energy, they will all be included.,
4. the 'first excited states', i.e., the states correspondign to the second minimum 
    energy,
5. The last state visited.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

SimulatedAnnealing::argerr="An adjacency matrix is expected at position 1, a hamiltonian 
at position 2, an optional formula for delta Hamiltonian in position 3, a real at 
position 4 (position 3 if optional delta Hamiltonian formula not given), a real at 
position 5 (position 4 if optional delta Hamiltonian formula not given),
and integers at the last two positions.";


(* ::Chapter:: *)
(*Exact Expectation Value Calculations*)


(* ::Section::Closed:: *)
(*Expectation Value*)


(* ::Item::Closed:: *)
(*ExactExpectationValue*)


(* :Usage Mesages: *)

ExactExpectationValue::usage="ExactExpectationValue[ensemble_List,degenValues_List,obsfunction_,Hamiltonian_,beta_]
Calculates the expectation value of the observable function(s). Energies in Boltzman weights are computed using the Hamiltonian. degeneracy is vector of functions that give the different weights for a given graph. Beta is inverse temperature. \[IndentingNewLine]Usage: ExactExpectationValue[{amat1, amat2},{{1.0,1.0},{1.0/3.0,1.0/5.0}},{AvgDeg, AvgDim},HIsing[#,-1.0,0.0]&,10.0]. \[IndentingNewLine]
Inputs are:
1. states = List, a list of adjacency matrix of a graphs
2. degenValues = List of list, a list of numbers representing the degeneracy or other 
	weight to be applied to each graph in the state. Each list within the list must have 
	the same length as states.
3. obsfunction = a list of functions of observables whose expectation values in the canonical
	ensemble is to be computed. In the overload this can be replaced by a list of list of 
	observable values.
4. Hamiltonian = a one-argument function of an adjacency matrix (e.g. HIsing[#,-1.0,0.0]&),  

5. beta = Real or a variable, the inverse temperature

Outputs a list of lists with as many items as the length of the degeneracy list. 
Each item in the output list has in itself the partition function, the free energy, 
and the expectation values of the observable functions, all of which are computed with 
the corresponding degeneracy weight.";
(* :Error Mesages: *)

ExactExpectationValue::argerr="Input should be of the form
	ExactExpectationValue[ensemble_List,degenValues_List,obsfunction_,Hamiltonian_,beta_],
	or ExactExpectationValue[ensemble_List,degenValues_List,obsValues_List,Hamiltonian_,beta_].";


(* ::Chapter:: *)
(*Metropolis MC*)


(* ::Section::Closed:: *)
(*Single Flip Metropolis*)


(* ::Item::Closed:: *)
(*GraphMetropolis*)


(* :Usage Mesages: *)

GraphMetropolis::usage="GraphMetropolis[seedAmat,beta,H, observables,maxSweep] 
runs Metropolic MC simulation on the seed input graph seedAmat at inverse temperature beta, 
with Hamiltonian H for maxSweep MC sweeps.
Inputs are:
1. seedAmat= seed adjacency matrix,\[IndentingNewLine]2. beta = inverse temperature,\[IndentingNewLine]3. hamiltonian = the hamiltonian which can depend on arbitrary parameters,\[IndentingNewLine]4. observables = a list of cuntions on graph adjacency matrices e.g. magnetization, average degree etc.,\[IndentingNewLine]5. maxSweep= the total number of MC sweeps where one sweep is N(N-1)/2 MC flip attempts.,
The function returns a lsit with two elements:
1. The adjacency matrix of the last graph visited in the MC run,
2. A list of values of the observables applied to each intermediate graph after 
1 MC sweep.";

(* :Error Mesages: *)

GraphMetropolis::argerr="A graph adjacency matrix was expected at position 1, 
a real number at position 2, a hamiltonian at position 3, a list of functions
at positions 4, and an integer at positions 5.";


(* ::Chapter:: *)
(*Parallel Tempering*)


(* ::Section::Closed:: *)
(*Parallel Tempering*)


(* ::Subsection::Closed:: *)
(*GraphSweepReplica*)


(* ::Item::Closed:: *)
(*GraphSweepReplica*)


(* :Usage Mesages: *)

GraphSweepReplica::usage="GraphSweepReplica[seedGraph_List,beta_Real,
	hamiltonian_[hparams___],delH_[delHparams___],
	NN_Integer,minEToBeat_Real,UnlabeledVerticesYes_Integer], overload
GraphSweepReplica[seedGraph_List,beta_Real,hamiltonian_[hparams___],NN_Integer,minEToBeat_Real,
	UnlabeledVerticesYes_Integer] performs NN MCMC sweeps on a seed graph state 
given by the adjacency matrix seedGraph
Inputs are:
1. seedGraph = List, a seed graph as an adjacency matrix, \[IndentingNewLine]2. beta = Real, inverse temperature beta, \[IndentingNewLine]3. hamiltonian[hparams] = a formula, the hamiltonian together with its inpur parameters, \[IndentingNewLine]4. delH[delHparam] = optional formula for delta E (when one edge is flipped to expedite computation),  \[IndentingNewLine]5. NN = Integer, number of sweeps NN,\[IndentingNewLine]6. minEToBeat = Real, a value for energy such that the lowest energy states with 
   energy lower than minEToBeat will be saved (for ground state search),\[IndentingNewLine]7. UnlabeledVerticesYes = Integer, 0 means no selection probability to make the graphs
   unlabeled, UnlabeledVerticesYes = 1 means graphs are unlabeled.   \[IndentingNewLine]\[IndentingNewLine]Outputs a list with two associations,\[IndentingNewLine]1. the minimum energy visited throughout the sweep states with that energy. 
    If multiple states have degenerate minimum energy, they will all be included. 
    It saves only non-isomorphic graphs. ,\[IndentingNewLine]2. The second association is the temperature, final graph, energy, and magnetization 
    at the end of the sweeps.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

GraphSweepReplica::argerr="A graph adjacency matrix is expected at position 1, 
a real number at position 2, a hamiltonian at position 3, an optional formula for 
the change in energy when one edge is toggled at position 4, integer at position 5, 
a real number at position 6, an integer at position 7.";


(* ::Subsection::Closed:: *)
(*GraphEquilibriate*)


(* ::Item::Closed:: *)
(*GraphEquilibriate*)


(* :Usage Mesages: *)

GraphEquilibriate::usage="GraphEquilibriate[seedGraph_List, beta_Real, 
	hamiltonian_[hparams___],delH_[delHparams___], UnlabeledVerticesYes_Integer] 
overload
	GraphEquilibriate[seedGraph_List, beta_Real, hamiltonian_[hparams___],
		UnlabeledVerticesYes_Integer] equilibriates an input graph configuration to the 
	temperature beta.,\[IndentingNewLine]Inputs are:,\[IndentingNewLine]1. seedGraph = List, adjacency matrix of the seed graph,\[IndentingNewLine]2. beta = Real, inverse temperature,\[IndentingNewLine]3. hamiltonian = formula, the hamiltonian, \[IndentingNewLine]4. delH = optional formula for delta E (when one edge is flipped to expedite computation),  \[IndentingNewLine]5. UnlabeledVerticesYes = Integer, 0 means no selection probability to make the 
   graphs unlabeled so graphs are labeled, 
   UnlabeledVerticesYes = 1 means graphs are unlabeled.,\[IndentingNewLine]\[IndentingNewLine]Depends on the function GraphSweepReplicas., \[IndentingNewLine]\[IndentingNewLine]Outputs a list with two associations:, ;, \[IndentingNewLine]The second is the final equilibriated state  \[IndentingNewLine]\[IndentingNewLine]Outputs a list with two associations,\[IndentingNewLine]1. the minimum energy visited throughout the equilibriation and the states with that 
    energy. If multiple states have degenerate minimum energy, they will all be included, 
    but if multiple identical adjacency graphs are found, only one unique adjacency graph 
    is kept. \[IndentingNewLine]2. The second association has the inverse temperature, equilibriation time, and the final 
   equilibriated state, which itself is an association which includes the adjacency matrix, 
   magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

GraphEquilibriate::argerr="A graph adjacency matrix is expected at position 1, 
a real number at position 2, a hamiltonian at position 3, an optional formula for 
the change in energy when one edge is toggled at position 4, integer at position 5.";

(* Built by concatenation rather than as one literal spanning several source lines: embedded
	newlines and tabs make a message render as misaligned columns once the `n` slots are
	filled in. *)
GraphEquilibriate::noconv=
	"Equilibriation did not converge within the budget of `1` sweeps. At exit the mean "<>
	"squared pairwise difference in `2` across the tracks -- the observable that agreed "<>
	"least well -- was `3`, against `4`, the value that quantity is expected to take once "<>
	"the tracks are equilibriated; convergence requires the former to be below the latter. "<>
	"Returning eqlT -> `5`, which is a lower bound rather than an estimate, and a state that "<>
	"may not be equilibriated. Raise $ECGravMaxEquilibriationSweeps, or start from a better "<>
	"seed graph.";

(* Concatenated for the same reason as ::noconv above. *)
GraphEquilibriate::nosignal=
	"None of the monitored observables (`1`) took more than one value over the whole "<>
	"comparison window, on any track, so there was nothing for the convergence test to "<>
	"read and equilibriation could not be checked. Returning eqlT -> `2` and converged -> "<>
	"True, since a chain whose observables never move is stationary in them, but note that "<>
	"this is the absence of evidence rather than evidence of mixing: if the chain is stuck, "<>
	"this is what that looks like. Pass observables that distinguish the states you care "<>
	"about, or check the returned state by hand.";


(* ::Subsection::Closed:: *)
(*GraphComputeCorrelationTime*)


(* ::Item::Closed:: *)
(*GraphComputeCorrelationTime*)


(* :Usage Mesages: *)

GraphComputeCorrelationTime::usage="GraphComputeCorrelationTime[seedGraph_List,beta_Real,
	hamiltonian_[hparams___], delH_[delHparams___], eqlT_Integer, minEToBeat_Real,
	EnergyOrMag_Integer,UnlabeledVerticesYes_Integer]
overload GraphComputeCorrelationTime[seedGraph_List,beta_Real,
	hamiltonian_[hparams___], delH_[delHparams___], eqlT_Integer, minEToBeat_Real,
	operators_,UnlabeledVerticesYes_Integer] and their versions with delH 
	takes an equilibriated graph with equilibriation time and graph configurations
	and computes the correlation time.,\[IndentingNewLine]
Depends on the functions: GraphSweepReplicas, CorrelationTime.,

Inputs are:,\[IndentingNewLine]1. seedGraph = List, adjacency matrix of the seed graph,\[IndentingNewLine]2. beta = Real, inverse temperature,\[IndentingNewLine]3. hamiltonian = formula, the hamiltonian, \[IndentingNewLine]4. delH = optional formula for delta E (when one edge is flipped to expedite computation),  \[IndentingNewLine]5. eqlT = Integer, equilibriation time.,
6. minEToBeat = Real, energy value such that if a state (or states ) are found with energy 
   lower than this, then the lowest such energy and configurations with that energy 
   will be saved.,\[IndentingNewLine]7. EnergyOrMag = Integer, a variable to specify whether energy or magnetization is used to 
    compute correlation time. If EnergyOrMag = 0, energy is used; if EnergyOrMag = 1, 
    magnetization is used. OR, 
7. operators = a list of functions whose correlation time is to be computed,\[IndentingNewLine]8. UnlabeledVerticesYes = Integer, 0 means no selection probability to make the 
   graphs unlabeled so graphs are labeled, 
   UnlabeledVerticesYes = 1 means graphs are unlabeled.,\[IndentingNewLine]\[IndentingNewLine]Outputs a list with two associations,\[IndentingNewLine]1. the minimum energy visited throughout the equilibriation and the states with that 
    energy. If multiple states have degenerate minimum energy, they will all be included, 
    but if multiple identical adjacency graphs are found, only one unique adjacency graph 
    is kept. \[IndentingNewLine]2. The second association has the inverse temperature, equilibriation time, 
    correlation timeand the final state visited, which itself is an association which 
    includes the adjacency matrix, 
   magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

GraphComputeCorrelationTime::argerr="A graph adjacency matrix is expected at position 1,
a real number at position 2, a hamiltonian at position 3, an optional formula for
the change in energy when one edge is toggled at position 4, and Integer at position 5,
a Real numbers at positions 6, and integers at positions 7 and 8. The integer at position
8 has to be 1 or 0 for unlabeled ot labeled graphs respectively";

GraphComputeCorrelationTime::stuck="All measured values of `1` are stuck at `2`.
	Autocorrelation time can not be computed; a default value of 2 will be assigned.";

GraphComputeCorrelationTime::alldefault="Autocorrelation time could not be computed
	since all measured values of the observables are stuck at the values `1`.
	Exiting with default values corrT -> 2, corrTValues -> `2`.";


(* ::Subsection::Closed:: *)
(*GraphMultiHistogram*)


(* ::Item::Closed:: *)
(*GraphMultiHistogram*)


(* :Usage Mesages: *)

GraphMultiHistogram::usage="GraphMultiHistogram[seedGraph_List,betaLow_Real,betaHigh_Real,
	hamiltonian_[hparams___],delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer] 
implements the Multiple Histogram Method for the graph models to get a smooth plot of the 
observables in the list of functions obs as a function of inverse temperature 
ranging from betaLow to betaHigh. It equilibriates, computes correlation times, 
takes NN independent measurements, then computes the free energies that will be used 
for the extrapolation.

The overload GraphMultiHistogram[seedGraph_List,betaLow_Real,betaHigh_Real,
	hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer] allows delH to be optional.

The overloads GraphMultiHistogram[seedGraph_List,betaTable_List,
	hamiltonian_[hparams___],delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer] and 
	GraphMultiHistogram[seedGraph_List,betaTable_List,
	hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer] take in a list of 
inverse temperatures in betaTable instead of just high and low.\[IndentingNewLine]
The overloads GraphMultiHistogram[seedGraph_List,beta_Real,
	hamiltonian_[hparams___], delH_[delHparams___], externalFieldTable_List, conjugateObs_,
	obs_, NN_Integer,UnlabeledVerticesYes_Integer] and 
	GraphMultiHistogram[seedGraph_List,beta_Real,
	hamiltonian_[hparams___],externalFieldTable_List,conjugateObs_, obs_,NN_Integer,
	UnlabeledVerticesYes_Integer] take in a single value
	of inverse temperatures and a list of external field parameters in externalFieldTable 
	to do multihistogram estrapolation with other external field parameters.

Depend on the functions: GraphSweepReplicas, GraphEquilibriate, 
	GraphComputeCorrelationTime, CorrelationTime.,

Inputs are:,\[IndentingNewLine]1. seedGraph = List, adjacency matrix of the seed graph,\[IndentingNewLine]2. betaLow = Real, lower bound of inverse temperature,
3. betaHigh = Real, upper bound of inverse temperature,
   Note: the overload takes a list of inverse temperature values instead of betaLow and 
   betaHigh\[IndentingNewLine]4. hamiltonian = formula, the hamiltonian, \[IndentingNewLine]5. delH = an optional formula for delta E (when one edge is flipped to expedite computation),  \[IndentingNewLine]6. obs = List, a list of formulas of observables that act on graph adjacency matrix and 
    output a number,
7. NN = Integer, number of independent measurements,
8. UnlabeledVerticesYes = Integer, 0 means no selection probability to make the 
   graphs unlabeled so graphs are labeled, 
   UnlabeledVerticesYes = 1 means graphs are unlabeled.,\[IndentingNewLine]\[IndentingNewLine]Outputs a list with three objects,\[IndentingNewLine]1. an association of the the minimum energy found throughout run and the states with that 
    energy. If multiple states have degenerate minimum energy, they will all be included, 
    but if multiple identical adjacency graphs are found, only one unique adjacency graph 
    is kept. 
2. an association with inverse temperatures (or external field parameters) as keys and 
   -(beta)*(free energy), energy values, (conjugate observable values), and observable 
	values at the corresponding value of the inverse 
   temperature beta (or the external field parameters).\[IndentingNewLine]3. the replicas in the last step, i.e., an association with inverse temperatures 
	(external field parameters) as keys 
   and graph states as values where a state is itself an association which includes the 
   adjacency matrix, magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

GraphMultiHistogram::argerr="A graph adjacency matrix is expected at position 1, real 
 numbers at positions 2 and 3 (alternatively a list of real numbers at position 2), 
 a hamiltonian at position 4 (or 3), an optional formula for the change in energy when one edge 
 is toggled at position 5 (or 4), a list of formulas corresponding to observables at 
 position 6 (or 5 or 4), an Integer at position 7 (or 6 or 5), and an Integer (0 or 1 only for 
 unlabeled ot labeled graphs respectively) at 
 position 8 (or 7 or 6)."


(* ::Subsection::Closed:: *)
(*GraphCEITempSchedule*)


(* ::Item::Closed:: *)
(*GraphCEITempSchedule*)


(* :Usage Mesages: *)

GraphCEITempSchedule::usage="GraphCEITempSchedule[seedGraph_List,betaLow_Real,betaHigh_Real,
	hamiltonian_[hparams___],delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer] 
finds a temperature schedule to be used for parallel tempering using the Constant
Entropy Increase (CEI) approach based on the paper by Sobo et. al. 
https://doi.org/10.1063/1.2907846. 

The overload GraphCEITempSchedule[seedGraph_List,betaLow_Real,betaHigh_Real,
hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer]  allows delH to be optional.\[IndentingNewLine]
The overloads GraphCEITempSchedule[seedGraph_List,betaTable_List,
		hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer] and 
	GraphCEITempSchedule[seedGraph_List,betaTable_List,
		hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer] take in a list of 
inverse temperatures in betaTable instead of just high and low.\[IndentingNewLine]
Depends on the functions: GraphSweepReplicas, GraphMultiHistogram (which also depends on
	GraphEquilibriate, GraphComputeCorrelationTime), CorrelationTime.,

Inputs are:,\[IndentingNewLine]1. seedGraph = List, adjacency matrix of the seed graph,\[IndentingNewLine]2. betaLow = Real, lower bound of inverse temperature,
3. betaHigh = Real, upper bound of inverse temperature,
   Note: the overload takes a list of inverse temperature values instead of betaLow and 
   betaHigh\[IndentingNewLine]4. hamiltonian[hparams] = formula, the hamiltonian, \[IndentingNewLine]5. delH[delHparams] = an optional formula for delta E (when one edge is flipped to expedite computation),  \[IndentingNewLine]6. obs = List, a list of formulas of observables that act on graph adjacency matrix and 
    output a number,
7. NN = Integer, number of independent measurements,
8. UnlabeledVerticesYes = Integer, 0 means no selection probability to make the 
   graphs unlabeled so graphs are labeled, 
   UnlabeledVerticesYes = 1 means graphs are unlabeled.,\[IndentingNewLine]\[IndentingNewLine]Outputs a list with four objects,\[IndentingNewLine]1. an association of the the minimum energy found throughout run and the states with that 
    energy. If multiple states have degenerate minimum energy, they will all be included, 
    but if multiple identical adjacency graphs are found, only one unique adjacency graph 
    is kept. 
2. temperature-entropy pairs,i.e., the temperature schedule,
3. an association with inverse temperatures as keys and -(beta)*(free energy), 
   energy values, and observable values at the corresponding value of the inverse 
   temperature beta.\[IndentingNewLine]4. the replicas in the last step, i.e., an association with inverse temperatures as keys 
   and graph states as values where a state is itself an association which includes the 
   adjacency matrix, magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

GraphCEITempSchedule::argerr="A graph adjacency matrix is expected at position 1, real 
 numbers at positions 2 and 3 (alternatively a list of real numbers at position 2), 
 a hamiltonian at position 4 (or 3), an optional formula for the change in energy when one edge 
 is toggled at position 5 (or 4), a list of formulas corresponding to observables at 
 position 6 (or 5 or 4), an Integer at position 7 (or 6 or 5), and an Integer (0 or 1 only for 
 unlabeled ot labeled graphs respectively) at 
 position 8 (or 7 or 6)."


(* ::Subsection::Closed:: *)
(*GraphCTLSchedule*)


(* :Usage Mesages: *)

GraphCTLSchedule::usage="GraphCTLSchedule[seedGraph_List,btTable_List,hamiltonian_[hparams___],
	delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer,
	numReplicas_Integer:$KernelCount] 
A program to find an inverse temperature schedule for parallel tempering that keeps 
acceptance ratios constant based on the constant thermodynamic length approach. 

The overload GraphCTLSchedule[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	delH_[delHparams___],externalFieldTable_List,conjugateObs_,obs_,NN_Integer,
	UnlabeledVerticesYes_Integer,numReplicas_Integer:$KernelCount] finds find a 
schedule of external field values at a fixed inverse temperature for use in parallel 
tempering that keeps acceptance ratios constant based on the constant thermodynamic 
length approach.

Overloads GraphCTLSchedule[betaMin_Real,betaMax_Real,minusbetaF_Association,energyValues_Association,
	numReplicas_Integer:$KernelCount] and 
	GraphCTLSchedule[bt_Real,hMin_Real,hMax_Real,minusbetaF_Association,conjugateFieldMeasurements_Association,
	numReplicas_Integer:$KernelCount] finds schedules for beta or external field using
	outputs of previous runs of GraphMultiHistogram or GraphCTLSchedule

remaining overloads without delH allow inputing a formula for delH to be optional.\[IndentingNewLine]\[IndentingNewLine]
Depends on the functions: GraphSweepReplicas, GraphMultiHistogram (which also depends on
	GraphEquilibriate, GraphComputeCorrelationTime), CorrelationTime.,

Inputs are:,\[IndentingNewLine]1. seedGraph = List, adjacency matrix of the seed graph,\[IndentingNewLine]2. betaLow = Real, lower bound of inverse temperature,
3. betaHigh = Real, upper bound of inverse temperature,
   Note: the overload takes a list of inverse temperature values instead of betaLow and 
   betaHigh\[IndentingNewLine]4. hamiltonian[hparams] = formula, the hamiltonian, \[IndentingNewLine]5. delH[delHparams] = an optional formula for delta E (when one edge is flipped to expedite computation),  \[IndentingNewLine]6. obs = List, a list of formulas of observables that act on graph adjacency matrix and 
    output a number,
7. NN = Integer, number of independent measurements,
8. UnlabeledVerticesYes = Integer, 0 means no selection probability to make the 
   graphs unlabeled so graphs are labeled, 
   UnlabeledVerticesYes = 1 means graphs are unlabeled.,\[IndentingNewLine]\[IndentingNewLine]Outputs a list with four objects,\[IndentingNewLine]1. an association of the the minimum energy found throughout run and the states with that 
    energy. If multiple states have degenerate minimum energy, they will all be included, 
    but if multiple identical adjacency graphs are found, only one unique adjacency graph 
    is kept. 
2. temperature-entropy pairs,i.e., the temperature schedule,
3. an association with inverse temperatures as keys and -(beta)*(free energy), 
   energy values, and observable values at the corresponding value of the inverse 
   temperature beta.\[IndentingNewLine]4. the replicas in the last step, i.e., an association with inverse temperatures as keys 
   and graph states as values where a state is itself an association which includes the 
   adjacency matrix, magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

GraphCTLSchedule::argerr="Input has to be of the form:
GraphCTLSchedule[seedGraph_List,btTable_List,hamiltonian_[hparams___],
	delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer,
	numReplicas_Integer:$KernelCount],
GraphCTLSchedule[seedGraph_List,btTable_List,hamiltonian_[hparams___],obs_,
	NN_Integer,UnlabeledVerticesYes_Integer, numReplicas_Integer:$KernelCount],
GraphCTLSchedule[betaMin_Real,betaMax_Real,minusbetaF_Association,energyValues_Association,
	numReplicas_Integer:$KernelCount],
GraphCTLSchedule[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	delH_[delHparams___],externalFieldTable_List,conjugateObs_,obs_,NN_Integer,
	UnlabeledVerticesYes_Integer,numReplicas_Integer:$KernelCount],
GraphCTLSchedule[seedGraph_List,bt_Real,hamiltonian_[hparams___],
	delH_[delHparams___],externalFieldTable_List,conjugateObs_,obs_,NN_Integer,
	UnlabeledVerticesYes_Integer,numReplicas_Integer:$KernelCount],
GraphCTLSchedule[bt_Real,hMin_Real,hMax_Real,minusbetaF_Association,conjugateFieldMeasurements_Association,
	numReplicas_Integer:$KernelCount]."


(* ::Subsection::Closed:: *)
(*GraphParallelTempering*)


(* ::Item::Closed:: *)
(*GraphParallelTempering*)


(* :Usage Mesages: *)

GraphParallelTempering::usage="GraphParallelTempering[seedGraph_List, btTable_List, 
	minEtoBeat_Real, hamiltonian_[hparams___],delH_[delHparams___],obs_,
	EnergyOrMag_Integer,NN_Integer, UnlabeledVerticesYes_Integer] 
Implements Parallel tempering algorithm on graph models at several different 
temperatures. It does not determine the temperature schedule here, it assumes the 
input beta list comes from a scheduler such as Constant Entropy Increase (CEI) based 
on the paper by Sobo et. al. https://doi.org/10.1063/1.2907846, or Constant Thermodynamic 
Length (CTL) temperature schedule. It equilibriates, computes correlation time, and 
then applies measurements. It does temperature swaps during the measurement step. 
It is parallelized so that equilibriation, computation of correlation time, and sweeps 
during measurement are all done in parallel. Temperature swaps are done on the master 
kernel. 

The overload GraphParallelTempering[seedGraph_List, btTable_List, minEtoBeat_Real,
          hamiltonian_[hparams___],obs_,EnergyOrMag_Integer,NN_Integer,
		UnlabeledVerticesYes_Integer]
	allows delH to be an optional formula.

The overloads GraphParallelTempering[inputReplicas_Association,minEtoBeat_Real,hamiltonian_[hparams___],delH_[delHparams___],obs_,
	EnergyOrMag_Integer,NN_Integer,UnlabeledVerticesYes_Integer] and
	GraphParallelTempering[inputReplicas_Association,minEtoBeat_Real,hamiltonian_,obs_,
	EnergyOrMag_Integer,NN_Integer,UnlabeledVerticesYes_Integer] take as an input an 
association of replicas with keys as integer numbers of the replica and values assiciations listing
the temperature, equilibriation time, correlation time, and equilibriated state of each
replica. In other words, the overload takes as an input, the 4th entry in the output of 
GraphParallelTempering. \[IndentingNewLine]
The other overloads implement parallel tempering for an external field schedule at a
fixed beta. 

Depends on the functions: GraphSweepReplicas, CorrelationTime, GraphEquilibriate, 
GraphComputeCorrelationTime. 

Inputs are:,\[IndentingNewLine]1. seedGraph = List, adjacency matrix of the seed graph,\[IndentingNewLine]2. btTable - a list of inverse temperatures for the replicas, ideally returned from the
   temperature schedule found by GraphCEITempSchedule
	Note, inputs 1 and 2 maybe replaced by the output association inputReplicas of a previous run
3. minEToBeat = Real, a value for energy such that the lowest energy states with 
   energy lower than minEToBeat will be saved (for ground state search),\[IndentingNewLine]4. hamiltonian = formula, the hamiltonian, \[IndentingNewLine]5. delH = a formula for delta E (when one edge is flipped to expedite computation),  \[IndentingNewLine]6. obs = List, a list of formulas of observables that act on graph adjacency matrix and 
    output a number,
7. EnergyOrMag = an integer (0 or 1) to specify whether to use energy (= 0) 
   or magnetization ( = 1) for computing correlation time,
8. NN = Integer, number of independent measurements,
9. UnlabeledVerticesYes = Integer, 0 means no selection probability to make the 
   graphs unlabeled so graphs are labeled, 
   UnlabeledVerticesYes = 1 means graphs are unlabeled.,\[IndentingNewLine]\[IndentingNewLine]Outputs a list with four objects,\[IndentingNewLine]1. an association of the the minimum energy found throughout run and the states with that 
    energy. If multiple states have degenerate minimum energy, they will all be included, 
    but if multiple identical adjacency graphs are found, only one unique adjacency graph 
    is kept. 
2. an association of each temperature with a table of energies, magnetizations, and 
   meausrements collected every t = numIndependentMeasurements/numberOfDataPoints = 
   (NN * corrT times)/numberOfDataPoints.,\[IndentingNewLine]3. the replicas in the last step, i.e., an association with inverse temperatures as keys 
   and graph states as values where a state is itself an association which includes the 
   adjacency matrix, magnetization, and energy, 
   i.e., <|'state'-><|'graph;->curAmat,`energy' \[Rule]hamiltonian[curAmat],
                      'mag'\[Rule]Total[Flatten[seedGraph]]*1.0/(vCount(vCount-1))|>.
4. an association of each replica and its history.

Both energy callbacks are passed as inert head[params] expressions and are called as
hamiltonian[adjacencyMatrix,hparams] and delH[adjacencyMatrix,delHparams,i,j]. Define the
head over all of those arguments and pass it partially applied, e.g.
h[am_List,j_Real,l_Real]:=HIsing[am,j,l] passed as h[-1.,0.], and
dH[am_List,j_Real,l_Real,a_Integer,b_Integer]:=delHIsing[am,j,l,a,b] passed as dH[-1.,0.].
A curried pure Function such as HIsing[#,-1.,0.]& matches the argument pattern but cannot
be called this way and is rejected with ECGrav::badham or ECGrav::baddelh.";

(* :Error Mesages: *)

GraphParallelTempering::argerr="Input must be of the form: 
	GraphParallelTempering[seedGraph_List, btTable_List,minEtoBeat_Real,hamiltonian_[hparams___],delH_[delHparams___],obs_,EnergyOrMag_Integer,NN_Integer,UnlabeledVerticesYes_Integer]
, GraphParallelTempering[seedGraph_List, btTable_List,minEtoBeat_Real,hamiltonian_[hparams___],obs_,EnergyOrMag_Integer,NN_Integer,UnlabeledVerticesYes_Integer]
, GraphParallelTempering[seedGraph_List, bt_Real,minEtoBeat_Real,hamiltonian_[hparams___],delH_[delHparams___],externalFieldTable_List,conjugateObs_,obs_,EnergyOrMag_Integer,NN_Integer,UnlabeledVerticesYes_Integer]
, GraphParallelTempering[seedGraph_List, bt_Real,minEtoBeat_Real,hamiltonian_[hparams___],externalFieldTable_List,conjugateObs_,obs_,EnergyOrMag_Integer,NN_Integer,UnlabeledVerticesYes_Integer]
, GraphParallelTempering[inputReplicas_Association,minEtoBeat_Real,hamiltonian_[hparams___],delH_[delHparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer]
, GraphParallelTempering[inputReplicas_Association,minEtoBeat_Real,hamiltonian_[hparams___],obs_,NN_Integer,UnlabeledVerticesYes_Integer]
, GraphParallelTempering[inputReplicas_Association,minEtoBeat_Real,hamiltonian_[hparams___],delH_[delHparams___],conjugateObs_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]
, GraphParallelTempering[inputReplicas_Association/;(L,minEtoBeat_Real,hamiltonian_[hparams___],conjugateObs_,obs_,NN_Integer,UnlabeledVerticesYes_Integer]"


(* ::Title:: *)
(*ECGrav Public Functions Involving PureComplexes*)


(* ::Chapter:: *)
(*Helper Functions*)


(* ::Section::Closed:: *)
(*General*)


(* ::Subsection::Closed:: *)
(*Basic Simplicial Complex Constructions*)


(* ::Item::Closed:: *)
(*FacetIncidenceMatrix*)


(* :Usage Messages: *)

FacetIncidenceMatrix::usage ="FacetIncidenceMatrix[facets_List], 
	overload FacetIncidenceMatrix[g_Graph] returns the 
	facet incidence matrix of the simplicial complex given by a list of its facets, or
	the graph g. ";

(* :Error Messages: *)

FacetIncidenceMatrix::argerr="Input should be a list of lists (i.e., list of facets) 
	or a graph object";


(* ::Item::Closed:: *)
(*FacetAdjacencyMatrix*)


(* :Usage Messages: *)

FacetAdjacencyMatrix::usage ="FacetAdjacencyMatrix[facets_List], 
	overload FacetIncidenceMatrix[g_Graph] returns the 
	facet adjacency matrix of the simplicial complex given by a list of its facets, or
	the graph g. ";

(* :Error Messages: *)

FacetAdjacencyMatrix::argerr="Input should be a list of lists (i.e., list of facets) 
	or a graph object";


(* ::Item::Closed:: *)
(*GraphFromCliques*)


(* :Usage Messages: *)

GraphFromCliques::usage="GraphFromCliques[facets] returns the graph or clique complex 
completion of the complex given by the list of its facets.";

(* :Error Messages: *)

GraphFromCliques::argerr="Input should be a list of lists (i.e., list of masimal cliques).";


(* ::Item::Closed:: *)
(*GraphFromFacetIncidence*)


(* :Usage Messages: *)

GraphFromFacetIncidence::usage="GraphFromFacetIncidence[incidenceMat]
	returns the graph or clique complex completion of the complex given a facet 
	incidence matrix.";

(* :Error Messages: *)

GraphFromFacetIncidence::argerr="Input should be a {0,1} matrix (the facet-incidence)";


(* ::Item::Closed:: *)
(*CliquesFromFacetIncidence*)


(* :Usage Messages: *)

CliquesFromFacetIncidence::usage="CliquesFromFacetIncidence[incidenceMat]
	returns the list of cliques of the clique complex completion of the complex given 
	as facet incidence matrix";

(* :Error Messages: *)

CliquesFromFacetIncidence::argerr="Input should be a {0,1} matrix (the facet-incidence)";


(* ::Item::Closed:: *)
(*ComplexFromFacetLabeledVertexList*)


(* :Usage Messages: *)

ComplexFromFacetLabeledVertexList::usage="
	ComplexFromFacetLabeledVertexList[facetLabeledVertices] returns the complex as 
	a normal list of facets given an input as list of vertices of the complex labeled 
	by the facets they belong to (i.e., a facet-labeled complex).";

(* :Error Messages: *)

ComplexFromFacetLabeledVertexList::argerr="Input should be a list of lists 
	(i.e., list of vertices given as a list of facets they belong to)";


(* ::Item::Closed:: *)
(*FacetLabeledVertexListFromComplex*)


(* :Usage Messages: *)

FacetLabeledVertexListFromComplex::usage="
	FacetLabeledVertexListFromComplex[facetsLst] returns the complex given a list of 
	vertices labeled by the facets they belong to (i.e., a facet-labeled complex)";

(* :Error Messages: *)

FacetLabeledVertexListFromComplex::argerr="Input should be a list of lists (i.e., list of facets)";


(* ::Item::Closed:: *)
(*PureComplexQ*)


(* :Usage Messages: *)

PureComplexQ::usage = "PureComplexQ[facets_List]returns True if the complex given as a list 
	of facets is pure.";

(* :Error Messages: *)

PureComplexQ::argerr="Input should be a list of lists (i.e., list of facets)";


(* ::Item::Closed:: *)
(*PureGraphQ*)


(* :Usage Messages: *)

PureGraphQ::usage = "PureGraphQ[g_Graph] returns True if the graph is pure.";

(* :Error Messages: *)

PureGraphQ::argerr="Input should be a graph object.";


(* ::Item::Closed:: *)
(*CliqueComplexQ*)


(* :Usage Messages: *)

CliqueComplexQ::usage = "CliqueComplexQ[facetsLst_List] returns True if the complex
	given as a list of facets is a clique complex.";

(* :Error Messages: *)

CliqueComplexQ::argerr="Input should be a list of lists (the facets).";


(* ::Subsection::Closed:: *)
(*Spheres, Balls, Links, Stars*)


(* ::Item::Closed:: *)
(*Sph*)


(* :Usage Messages: *)

Sph::usage="Sph[Amat,i,r]; overloads Sph[g,i,r] returns the 
	sphere of radius r, at vertex i, in the graph given by Amat or g. 
	If r is not specified it gives the unit sphere with r=1.";

(* :Error Messages: *)

Sph::argerr="Adjacency matrix or a graph object is expected at
	position 1 and integers at positions 2 and/or 3.";

Sph::vtxnotfound="The vertex `1` is not in the graph. Returning {}.";


(* ::Item::Closed:: *)
(*Bll*)


(* :Usage Messages: *)

Bll::usage="Bll[Amat,i,r]; overloads Bll[g,i,r] returns the 
	ball of radius r, at vertex i, in the graph given by Amat or g. 
	If r is not specified it gives the unit ball with r=1.";

(* :Error Messages: *)

Bll::argerr="Adjacency matrix or a graph object is expected at
	position 1 and integers at positions 2 and/or 3.";

Bll::vtxnotfound="The vertex `1` is not in the graph. Returning {}.";


(* ::Item::Closed:: *)
(*Lnk*)


(* :Usage Messages: *)

Lnk::usage="Lnk[facetsLst,face], returns the Link of the face in the complex given as 
	a list facetsLst of facets.";

(* :Error Messages: *)

Lnk::argerr="A list of facets is expected at 
	position 1 a list at position 2.";


(* ::Item::Closed:: *)
(*Str*)


(* :Usage Messages: *)

Str::usage="Str[facetsLst,face], returns the closed star Str the face in the complex 
	given as a list facetsLst of facets.";

(* :Error Messages: *)

Str::argerr="A list of facets is expected at 
	position 1 a list at position 2.";


(* ::Subsection::Closed:: *)
(*Basic Graph Observables*)


(* ::Item::Closed:: *)
(*Deg*)


(* :Usage Messages: *)

Deg::usage="Deg[Amat,i] gives the degree of vertex i in the graph with adjacency matrix Amat.";

(* :Error Messages: *)

Deg::argerr="An adjacency matrix is expected at 
	position 1.";


(* ::Item::Closed:: *)
(*AvgDeg*)


AvgDeg::usage="AvgDeg[Amat_List], AvgDeg[g_Graph] gives the average degree of a 
	graph given as an adjacency 
	matrix Amat or a graph g";

(* :Error Messages: *)

AvgDeg::argerr="Input must be of the form AvgDeg[Amat_List] of AvgDeg[g_Graph].";


(* ::Item::Closed:: *)
(*FacetDeg*)


(* :Usage Messages: *)

FacetDeg::usage="FacetDeg[facetsLst,face] gives the facet degree of the face
	given the complex as a list facetsLst of facets."; 				

(* :Error Messages: *)

FacetDeg::argerr="A list (of facets) is expected at position 1 an integer (a vertex) 
	at position 2 .";


(* ::Item::Closed:: *)
(*HyperDeg*)


(* :Usage Messages: *)
 				

HyperDeg::usage="HyperDeg[g_Graph,clq_List];
	overload HyperDeg[Amat_List,clq_List], HyperDeg[facetsLst_List,clq_List],
	computes the degree of the clique (face) clq given as a list of vertices in the 
	complex given as a graph g, adjacency matrix Amat, or a list facetsLst of facets. 
	It checks whether or not the input is a clique (face) in the graph/complex.";

(* :Error Messages: *)

HyperDeg::argerr="A graph, an adjacency matrix, or a list of facets is expected at 
	position 1 and a list at position 2.";


(* ::Item::Closed:: *)
(*FVector*)


(* :Usage Messages: *)
 				

FVector::usage="FVector[g_Graph];
	overload FVector[Amat_List], FVector[facetsLst_List],
	Computes the f-vector of the graph/complex which is the vector {f0,f1,...} of 
	the number of vertices, number of edges, etc.";

(* :Error Messages: *)

FVector::argerr="A graph, an adjacency matrix, or a list of facets is expected as as
	input.";


(* ::Item::Closed:: *)
(*GVolume*)


(* :Usage Messages: *)

GVolume::usage="GVolume[Amat_List,Dmat_List,r_Integer,i_Integer]; overloads
	GVolume[g_Graph,Dmat_List,r_Integer,i_Integer], 
	GVolume[Amat_List,r_Integer,i_Integer]
	GVolume[g_Graph,r_Integer,i_Integer] Computes volume of 
	the ball (number of vertices within) in the graph with adjacency matrix Amat, 
	centered at node i and having radius r, where the distance between nodes is given 
	by Dmat. GVolume[Amat_List,r_Integer,i_Integer], GVolume[Amat_List,r_Integer,i_Integer]
	do the computation with the distance matrix computed from the graph.";
	

(* :Error Messages: *)

GVolume::argerr="Input must be of the form GVolume[Amat_List,Dmat_List,r_Integer,i_Integer]
	, GVolume[g_Graph,Dmat_List,r_Integer,i_Integer], 
	GVolume[Amat_List,r_Integer,i_Integer], GVolume[g_Graph,r_Integer,i_Integer].";


(* ::Item::Closed:: *)
(*KFaceDistance*)


(* :Usage Messages: *)

KFaceDistance::usage="KFaceDistance[facetsLst_List,face1_List,face2_List] 
	Given a simplicial complex as a list of facets and two equal cardinality faces, 
	face1 and face2 of cardinality k in the graph, it computes the minimum k-face 
	distance between them where the path steps on only (k+1)-dimensional faces joined 
	along k-dimensional subfaces. It  uses the built in Mathematica GraphDistance 
	function on a graph constructed such that the vertices are all k-faces and edges 
	between them are 1 if two faces are contained in a bigger face and 0 otherwise.";

(* :Error Messages: *)

KFaceDistance::argerr="Input must be of the form KFaceDistance[facetsLst_List,face1_List,face2_List].";


(* ::Item::Closed:: *)
(*KFaceDistanceMatrix*)


(* :Usage Messages: *)

KFaceDistanceMatrix::usage="KFaceDistanceMatrix[facetsLst_List,k_Integer]; Given a 
	simplicial complex as a list of facets dimensionality k, it computes the 
	k-face distance matrix between every pair of faces of cardinality facedim, where each path 
	steps on only (k+1)-dimensional faces.  It  uses the built in Mathematica 
	GraphDistanceMatrix function on a graph constructed such that the vertices 
	are all clqdim-faces and edges between them are 1 if two faces are contained in a 
	bigger clique and 0 otherwise.";

(* :Error Messages: *)

KFaceDistanceMatrix::argerr="Input must be of the form KFaceDistanceMatrix[facetsLst_List,k_Integer].";


(* ::Item::Closed:: *)
(*KpathConnectedComponents*)


(* :Usage Messages: *)

KpathConnectedComponents::usage="KpathConnectedComponents[facetsLst_List,k_Integer]; 
	Overload KpathConnectedComponents[g_Graph,k_Integer]: Given a complex as a list 
	facetsLst of facets or as a graph g,  it computes the 
	number of components connected by (k+1)-paths where those k-faces that form a 
	component connected by a path of (k+1)-dimensional faces will form a component.  
	It  uses the built in Mathematica KVertexConnectedComponents function on a 
	graph constructed such that the vertices are all k-faces and edges between 
	them are 1 if two faces are contained in a bigger clique and 0 otherwise.";

(* :Error Messages: *)

KpathConnectedComponents::argerr="Input must be of the form 
	KpathConnectedComponents[facetsLst_List,k_Integer] or 
	KpathConnectedComponents[g_Graph,k_Integer].";


(* ::Item::Closed:: *)
(*ConnectedComplexComponents*)


(* :Usage Messages: *)

ConnectedComplexComponents::usage="ConnectedComplexComponents[facelst_List] 
	gives a list of simplicial sub-complexes of the complex facelst which are 
	connected. E.g.ConnectedComplexComponents[{{1,2},{3,4}}] = {{{1,2}},{{3,4}}}";

(* :Error Messages: *)

ConnectedComplexComponents::argerr="Input mys be of he form 
	ConnectedComplexComponents[facelst_List].";


(* ::Item::Closed:: *)
(*FractionInLargestComponent*)


(* :Usage Messages: *)

FractionInLargestComponent::usage="FractionInLargestComponent[g_Graph], 
	FractionInLargestComponent[Amat_List] Given a graph, it computes the ratio 
	of the number of vertices in the largest connected component to the total 
	number of vertices in the graph.";

(* :Error Messages: *)

FractionInLargestComponent::argerr="Input must be of the form 
	FractionInLargestComponent[g_Graph] or FractionInLargestComponent[Amat_List].";


(* ::Item::Closed:: *)
(*FractionInLargestKPathComponent*)


(* :Usage Messages: *)

FractionInLargestKPathComponent::usage="FractionInLargestKPathComponent[g_Graph,k_Integer] 
	overloads FractionInLargestKPathComponent[g_Graph,k_Integer], 
	FractionInLargestKPathComponent[amat_List,k_Integer], 
	FractionInLargestKPathComponent[amat_List,k_Integer]: Given a graph, it computes the 
	ratio of the number of vertices in the largest 
	k-path connected component to the total number of vertices in the complex. If k is 
	not specified, it uses k = size of the largest maximal clique.";

(* :Error Messages: *)

FractionInLargestKPathComponent::argerr="Input must be of the form 
	FractionInLargestKPathComponent[facetsLst_List,k_Integer], 
	FractionInLargestKPathComponent[g_Graph,k_Integer], 
	FractionInLargestKPathComponent[facetsLst_List], or 
	FractionInLargestKPathComponent[g_Graph].";


(* ::Item::Closed:: *)
(*CliqueOrder*)


(* :Usage Messages: *)
	
CliqueOrder::usage="CliqueOrder[g_Graph] Given a graph, it gives the facet order of the graph, 
	i.e., the number of maximal cliques of the graph.";


(* :Error Messages: *)

CliqueOrder::argerr="Input must be of the form CliqueOrder[g_Graph].";


(* ::Item::Closed:: *)
(*FacetOrder*)


(* :Usage Messages: *)

FacetOrder::usage="FacetOrder[facetsLst_List] gives the facet order of the complex
fiven as a list of facets. ";

(* :Error Messages: *)

FacetOrder::argerr="Input must be of the form FacetOrder[facetsLst_List].";


(* ::Item::Closed:: *)
(*EulerChi*)


(* :Usage Messages: *)
	
EulerChi::usage="EulerChi[facetsList_List]; overloads EulerChi[Amat_List],
	EulerChi[g_Graph] gives the Euler Characteristic of the 
	complex given as a list of facets facetsLst, a graph g, .";

(* :Error Messages: *)

EulerChi::argerr="Input must be of the form EulerChi[facetsList_List].";


(* ::Item::Closed:: *)
(*RmsPurity*)


(* :Usage Messages: *)

RmsPurity::usage="RmsPurity[facetsLst_List]; overload RmsPurity[g_Graph], 
	RmsPurity[Amat_List] computes the 
	mean purity = 1/|F| sum_{fi,fj}(dim(f_i)-dim(fj))^2 of the complex facetsLst,
	graph g or graph Amat, 
	where fi, fj run over all facets.
	If there's just one facet it outputs 0.";

(* :Error Messages: *)

RmsPurity::argerr="Input must be of the form RmsPurity[g_Graph] or RmsPurity[facetsLst_List].";


(* ::Item::Closed:: *)
(*Branchedness*)


(* :Usage Messages: *)

Branchedness::usage="Branchedness[g_Graph,n_Integer]; overloads
	Branchedness[g_Graph], Branchedness[Amat_List,n_Integer],
	Branchedness[Amat] Given a graph and a dimension n, 
	it computes the rms value of the degree of the n-facets, minus 3/2. This quantity will be 
	zero for pseudo manifolds if n is the purity of the graph minus one since all codimension 1 
	faces have degree 1 or 2. ";


(* :Error Messages: *)

Branchedness::argerr="Input must be of the form Branchedness[g_Graph,n_Integer].";


(* ::Subsection::Closed:: *)
(*Graph Dimensions*)


(* ::Item::Closed:: *)
(*AvgKDim*)


(* :Usage Messages: *)

AvgKDim::usage = "AvgKDim[g_Graph] Computes the inductive (Knill) dimension 
	of a graph as defined by Oliver Knill.";

(* :Error Messages: *)

AvgKDim::argerr="A graph adjacency matrix or graph was expected.";


(* ::Item::Closed:: *)
(*SpectralDim*)


(* :Usage Messages: *)

SpectralDim::usage = "SpectralDim[Amat_List,s_Integer]; overload
	SpectralDim[g_Graph,s_Integer]Spectral dimension for all nodes at step s.";

(* :Error Messages: *)

SpectralDim::argerr="A graph adjacency matrix or graph was expected at position 1 and an integer greater than 1 at position 2.";


(* ::Item::Closed:: *)
(*HausdorffDim*)


(* :Usage Messages: *)

HausdorffDim::usage = " HausdorffDim[Amat_List,Dmat_List, s_Integer]; overloads
	HausdorffDim[Amat_List, s_Integer], 
	HausdorffDim[g_Graph,Dmat_List, s_Integer],
	HausdorffDim[g_Graph, s_Integer] computes the Hausdorff
	dimension for the first node at step s when the distance between nodes is 
	given by the matrix Dmat. If Dmat is not given, it computes and uses the 
	unweighted graph distance matrix";

(* :Error Messages: *)

HausdorffDim::argerr="A graph adjacency matrix or graph was expected at position 1, an optional square matrix at position 2, and an integer greater than 1 at position 3.";


(* ::Subsection::Closed:: *)
(*Pure and Geometric Complexes*)


(* ::Item::Closed:: *)
(*ClosedDGraphQ*)


(* :Usage Messages: *)

ClosedDGraphQ::usage = "ClosedDGraphQ[g_Graph] returns True if the clique complex of the graph
g is a combinatorial manifold without boundary (a closed d-graph): the unit sphere of every
vertex is recursively a closed d-graph. This is the topological manifold condition only; a
graph whose clique complex is a torus returns True. Use DSphereQ to test for a sphere.";

(* :Error Messages: *)

ClosedDGraphQ::argerr="Input should be a graph.";


(* ::Item::Closed:: *)
(*DSphereQ*)


(* :Usage Messages: *)

DSphereQ::usage = "DSphereQ[g_Graph] returns True if the clique complex of the graph g is a
combinatorial manifold homeomorphic to a sphere: a closed d-graph (see ClosedDGraphQ) whose
clique complex additionally has the Euler characteristic of a sphere of the same dimension
(chi = 1 + (-1)^d). For surfaces this is exact (only the 2-sphere returns True; a torus
returns False); in dimension >= 3 it is a necessary homology-level condition only.";

(* :Error Messages: *)

DSphereQ::argerr="Input should be a graph.";


(* ::Item::Closed:: *)
(*DGraphQ*)


(* :Usage Messages: *)

DGraphQ::usage = "DGraphQ[g_Graph] returns true of the graph is a geometric 
	graph and false if not.";

(* :Error Messages: *)

DGraphQ::argerr="A graph was expected.";


(* ::Item::Closed:: *)
(*DGraphBoundary*)


(* :Usage Messages: *)

DGraphBoundary::usage = "DGraphBoundary[g_Graph] or DGraphBoundary[Amat_List] outputs 
	the induced subgraph 
	which is the boundary of the graph g given as a graph object g or adjacency matrix 
	Amat. 
	The boudnary of a geometric graph is the induced subgraph over boundary 
	nodes, i.e. those vertices whose unit spheres are contractible, or paths.";

(* :Error Messages: *)

DGraphBoundary::argerr="Input must be a graph or an adjacency matrix.";


(* ::Item::Closed:: *)
(*CombinatorialManifoldQ*)


(* :Usage Messages: *)

CombinatorialManifoldQ::usage = "CombinatorialManifoldQ[facetsLst_List] returns True if 
	the complex is a combinatorial manifold.";

(* :Error Messages: *)

CombinatorialManifoldQ::argerr="Input should be a list of facets.";


(* ::Item::Closed:: *)
(*ClosedCombinatorialManifoldQ*)


(* :Usage Messages: *)

ClosedCombinatorialManifoldQ::usage = "ClosedCombinatorialManifoldQ[facetsLst_List] returns
True if the pure simplicial complex (given as a list of facets) is a combinatorial manifold
without boundary, i.e. a closed combinatorial manifold: connected, with every codimension-1
face contained in exactly two facets and every vertex link recursively a closed combinatorial
manifold. This is the topological manifold condition only; a torus, Klein bottle, or sphere
all return True. Use CombinatorialSphereQ to test whether the manifold is a sphere.";

(* :Error Messages: *)

ClosedCombinatorialManifoldQ::argerr="Input should be a list of facets.";


(* ::Item::Closed:: *)
(*CombinatorialSphereQ*)


(* :Usage Messages: *)

CombinatorialSphereQ::usage = "CombinatorialSphereQ[facetsLst_List] returns True if the pure
simplicial complex (given as a list of facets) is a combinatorial sphere: a closed combinatorial
manifold (see ClosedCombinatorialManifoldQ) that in addition has the Euler characteristic of a
sphere of the same dimension, chi = 1 + (-1)^d. For surfaces this is exact (only the 2-sphere
returns True; a torus, Klein bottle, or projective plane return False); in dimension >= 3 it is
a necessary homology-level condition only.";

(* :Error Messages: *)

CombinatorialSphereQ::argerr="Input should be a list of facets.";


(* ::Item::Closed:: *)
(*OrientableCombinatorialManifoldQ*)


(* :Usage Messages: *)

OrientableCombinatorialManifoldQ::usage = "OrientableCombinatorialManifoldQ[facetsLst_List] 
returns True if the complex assumed to be a pseudo combinatorial manifold is 
orientable.";

(* :Error Messages: *)

OrientableCombinatorialManifoldQ::argerr="Input should be a list of facets.";


(* ::Item::Closed:: *)
(*CombinatorialBoundary*)


(* :Usage Messages: *)

CombinatorialBoundary::usage = "CombinatorialBoundary[facetsLst_List] returns the 
	boundary of the complex assuming it is a combinatorial manifold.";

(* :Error Messages: *)

CombinatorialBoundary::argerr="Input should be a list of facets representing a
	combinatorial manifold.";

CombinatorialBoundary::notmanifold="Attempting to find the combinatorial boundary of a
	complex that is not a combinatorial manifold. Returning $Failed.";


(* ::Item::Closed:: *)
(*CountHoles*)


(* :Usage Messages: *)

CountHoles::usage = "CountHoles[g_Graph,k_Integer] gives the number of k-holes,
	i.e., the k-th Betti number of a graph. Uses Mathematica's built in 
	ResourceFunction[``BettiNumbers``]";

(* :Error Messages: *)

CountHoles::argerr="A list of lists or graph is expected at position 1 and 
	an optional integer at position 2.";


(* ::Item::Closed:: *)
(*Visualize2DComplex*)


(* :Usage Messages: *)

Visualize2DComplex::usage = "Visualize2DComplex[facetsLst_List] gives 3D rendering of 
	a 2D pure complex.";

(* :Error Messages: *)

Visualize2DComplex::argerr="Input should be a list of facets representing a 
	combinatorial manifold.";


(* ::Subsection::Closed:: *)
(*Counting Pure Complexes*)


(* ::Item::Closed:: *)
(*RankComb*)


(* :Usage Messages: *)

RankComb::usage = " RankComb[set_List,numLabels_Integer]: \[IndentingNewLine]Given a sorted set which is a subset of the set {0,1,,...,numLabels-1}, 
it assigns a unique integer to the set between 0 and numLabels choose 
length(set).\[IndentingNewLine]E.g., Rank[{0,1},3] = 0, and  Rank[{1,2},3]=3.";

RankComb::argerr="Input should be of the form RankComb[set_List,numLabels_Integer].";


(* ::Item::Closed:: *)
(*UnrankComb*)


(* :Usage Messages: *)

UnrankComb::usage = "UnrankComb[l_Integer,numLabels_Integer,setSize_Integer]: 
Given an integer l between 0 and numLabels choose setSize-1,
it assigns a unique sorted set which is a sorted setSize-subset of 
{0,1,,...,numLabels-1}.\[IndentingNewLine]E.g., UnrankComb[0,3,2] = {0,1}, and  UnrankComb[2,3,2] = {1,2}.";

UnrankComb::argerr="Input should be of teh form 
	UnrankComb[l_Integer,numLabels_Integer,setSize_Integer].";


(* ::Item::Closed:: *)
(*NumVertexLabeledPureComplexes*)


(* :Usage Messages: *)

NumVertexLabeledPureComplexes::usage = "NumVertexLabeledPureComplexes[p_Integer, M_Integer, n_Integer]
gives the number of vertex labeled pure simplicial complexes of purity p, facet
order M, and number of vertices n -- that is, the number of M-element sets of p-subsets
of [n] whose union is all of [n]. It is zero unless p <= n <= p M. Computes one whole
row of n at a time, working up from M = 1, and memoizes the rows for the session.
NumVertexLabeledPureComplexes[p_Integer, M_Integer] gives the number of pure simplicial complexes
of purity p and facet order M, summed over the vertex count.
M = 0 is the empty complex: NumVertexLabeledPureComplexes[p, 0, n] is 1 at n = 0 and 0 otherwise,
and NumVertexLabeledPureComplexes[p, 0] is 1. Negative M gives 0.
See NumFacetLabeledPureComplexes for the count with the facets labelled and the vertices not,
and NumUnlabeledPureComplexes for the count with neither labelled.";

NumVertexLabeledPureComplexes::argerr="Input should be of the form NumVertexLabeledPureComplexes[p_Integer, M_Integer, n_Integer]
	or NumVertexLabeledPureComplexes[p_Integer, M_Integer].";


(* ::Item::Closed:: *)
(*NumPureComplexes*)


(* :Usage Messages: *)

(* Former name of NumVertexLabeledPureComplexes, kept as an alias. The argument order is
   unchanged -- (purity, facet order, vertex count), with the facet order now written M rather
   than q -- so old calls keep returning what they always did. *)

NumPureComplexes::usage = "NumPureComplexes is the former name of NumVertexLabeledPureComplexes
and is kept as an alias for it, with the same arguments in the same order. Prefer
NumVertexLabeledPureComplexes, which says which of the two labellings it counts.";


(* ::Item::Closed:: *)
(*NumFacetLabeledPureComplexes*)


(* :Usage Messages: *)

NumFacetLabeledPureComplexes::usage = "NumFacetLabeledPureComplexes[p_Integer, M_Integer, n_Integer]
gives the number of facet labeled pure simplicial complexes of purity p, facet order M, and
number of vertices n: M labelled, pairwise distinct p-subsets of an n-element vertex set whose
union is the whole set, counted up to permutation of the vertices, which are not labelled.
Equivalently the number of separating (p, n, M) incidence tableaux -- multisets of n nonempty
subsets of [M] in which every label occurs in exactly p of them and no two labels occupy the
same set of rows.
NumFacetLabeledPureComplexes[p_Integer, M_Integer] sums that over the vertex count.
It is zero unless p <= n <= p M and Binomial[n, p] >= M.
See NumVertexLabeledPureComplexes for the count with the vertices labelled instead, and
NumUnlabeledPureComplexes for the count with neither labelled.";

NumFacetLabeledPureComplexes::argerr="Input should be of the form NumFacetLabeledPureComplexes[p_Integer, M_Integer, n_Integer]
	or NumFacetLabeledPureComplexes[p_Integer, M_Integer].";


(* ::Item::Closed:: *)
(*NumUnlabeledPureComplexes*)


(* :Usage Messages: *)

NumUnlabeledPureComplexes::usage = "NumUnlabeledPureComplexes[p_Integer, M_Integer, n_Integer]
gives the number of fully unlabeled pure simplicial complexes of purity p, facet order M, and
number of vertices n: M pairwise distinct p-subsets of an n-element vertex set whose union is the
whole set, counted up to relabelling the vertices and permuting the facets. These are the
isomorphism classes, so this is the count the two labelled ones both refine.
NumUnlabeledPureComplexes[p_Integer, M_Integer] sums that over the vertex count.
It is zero unless p <= n <= p M and Binomial[n, p] >= M.
It is computed by Burnside over the vertex permutations and so sums over the cycle types of
S_n: cost grows like PartitionsP[n], making n, rather than p or M, the limiting argument.
See NumVertexLabeledPureComplexes and NumFacetLabeledPureComplexes for the two counts that keep
one of the two labellings.";

NumUnlabeledPureComplexes::argerr="Input should be of the form NumUnlabeledPureComplexes[p_Integer, M_Integer, n_Integer]
	or NumUnlabeledPureComplexes[p_Integer, M_Integer].";


(* ::Section::Closed:: *)
(*Choosing Isomorphism Classes*)


(* ::Subsection::Closed:: *)
(*Graph Isomorphism Classes*)


(* ::Item::Closed:: *)
(*GraphIsContained*)


(* :Usage Messages: *)

GraphIsContained::usage="GraphIsContained[glist_List,g_Graph]; overload 
	GraphIsContained[AmatList_List,Amat_List], checks whether or not the graph 
	(either graph object g or adjacency matrix Amat)is isomorphic to any graphs in glist 
	(respectively Amatlist), returning True after the first occurance of isomorphism and 
	False otherwise. The elements of glist and AmatList must be exclusively graphs or 
	exclusively graph adjacency matrices respectively.";

GraphIsContained::argerr="Input should be of the form 
	GraphIsContained[glist_List,g_Graph] or GraphIsContained[AmatList_List,Amat_List].";


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicGraphs*)


(* :Usage Messages: *)

ChooseNonIsomorphicGraphs::usage="Given a list, multiple lists, or list of sublists GList, either exclusively of graphs or 
	exclusively of graph adjacency matrices, ChooseNonIsomorphicGraphs[GList_List] returns a list 
	of non-isomorphic graphs. When GList is a list of graphs/graph adjacency matrices,
	ChooseNonIsomorphicGraphs[GList_List] returns the largest list of non-isomorphic graphs 
	from GList. When GList is multiple lists or a list of sublists of graphs/graph adjacency 
	matrices, ChooseNonIsomorphicGraphs[GList_List] returns a list where any two graphs 
	originating from different sublists are non-isomorphic, but two graphs originating from the 
	same sublist are not checked for isomorphism.";
	
ChooseNonIsomorphicGraphs::argerr="Input should be a sequence, a list, a sequence of lists, 
	or list of lists of graphs only or adjacency matrices only.";


(* ::Subsection::Closed:: *)
(*Clique Complex Isomorphism Classes*)


(* ::Item::Closed:: *)
(*IsContainedClqComp*)


(* :Usage Messages: *)

IsContainedClqComp::usage="Given a list of clique complexes cList (each given as a list of 
	maximal cliques), and single other clique complex cmpx, 
	IsContainedClqComp[cList_List, cmpx_List] checks whether cmpx is isomorphic to any complexes 
	in cList, returning True after the first occurance of isomorphism and False otherwise.";

IsContainedClqComp::argerr="Input should be of the form 
	IsContainedClqComp[cList_List,cmpx_List], where cList is a list of list of lists (a list of 
	clique complexes each given as a list of facets), and 
	cmpx is a list of lists (facets).";


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicClqComplexes*)


(* :Usage Messages: *)

ChooseNonIsomorphicClqComplexes::usage="Given a sequence, a list, multiple lists, or 
	a list of list (cList) of clique complexes (each given as a list of maximal cliques),
	ChooseNonIsomorphicClqComplexes[cList_List] returns a list of non-isomorphic clique complexes. 
	When cList is a list of clique complexes, ChooseNonIsomorphicClqComplexes[cList] returns the 
	largest list of non-isomorphic complexes from cList. 
	When cList is a sequence of lists or a list of lists of clique complexes, 
	ChooseNonIsomorphicClqComplexes[cList] returns a single list where 
	any two complexes originating from different sublists are non-isomorphic, but two complexes 
	originating from the same sublist are not checked for isomorphism.";

ChooseNonIsomorphicClqComplexes::argerr="Input should be a sequence, a list, a sequence of lists
	or a list of lists of clique complexes only.";


(* ::Subsection::Closed:: *)
(*Simplicial Complex Isomorphism Classes*)


(* ::Item::Closed:: *)
(*IsomorphicSimplicialComplexQ*)


(* :Usage Messages: *)

IsomorphicSimplicialComplexQ::usage="IsomorphicSimplicialComplexQ[c1_List, c2_List] checks whether the 
	pure simplicial complexes c1 and c2 (each given as a list of facets) are isomorphic under 
	vertex labeling, returing True if isomorphic.";

IsomorphicSimplicialComplexQ::argerr="Input should be of the form 
	IsomorphicSimplicialComplexQ[c1_List, c2_List], where c1 and c2 are each a list of list (the list
	of facets of the complexes).";


(* ::Item::Closed:: *)
(*IsContainedSimplicialComp*)


(* :Usage Messages: *)

IsContainedSimplicialComp::usage="Given a list of simplicial complexes xList (each given 
	as a list of facets), and single other simplicial comples cmpx, 
	IsContainedSimplicialComp[xList_List, cmpx_List] checks whether cmpx is isomorphic to any 
	complexes in xList, returning True after the first occurance of isomorphism and False 
	otherwise.";

IsContainedSimplicialComp::argerr="Input should be of the form 
	IsContainedSimplicialComp[xList_List, cmpx_List], where xList is a list of list of lists 
	(a list of simplicial complexes each given as a list of facets), and 
	cmpx is a list of lists (facets).";
	


(* ::Item::Closed:: *)
(*ChooseNonIsomorphicSimplicialComplexes*)


(* :Usage Messages: *)

ChooseNonIsomorphicSimplicialComplexes::usage="Given a sequence, a list, multiple lists, or 
	a list of list (xList) of simplicial complexes (each given as a list of facets),
	ChooseNonIsomorphicSimplicialComplexes[xList_List] returns a list of non-isomorphic 
	simplicial complexes. 
	When xList is a list of simplicial complexes, ChooseNonIsomorphicClqComplexes[xList] returns 
	the largest list of non-isomorphic complexes from xList. 
	When xList is a sequence of lists or a list of lists of simplicial complexes, 
	ChooseNonIsomorphicSimplicialComplexes[xList] returns a single list where 
	any two complexes originating from different sublists are non-isomorphic, but two complexes 
	originating from the same sublist are not checked for isomorphism.";

ChooseNonIsomorphicSimplicialComplexes::argerr="Input should be a sequence, a list, 
	a sequence of lists or a list of lists of simplicial complexes only.";
	


(* ::Section::Closed:: *)
(*Automorphism Groups and Orders*)


(* ::Subsection::Closed:: *)
(*Simplicial Complex Automorphism and Facet Automorphism*)


(* ::Item::Closed:: *)
(*SimplicialComplexAutomorphismGroupOrder*)


(* :Usage Messages: *)

SimplicialComplexAutomorphismGroupOrder::usage="SimplicialComplexAutomorphismGroupOrder[facetsLst_List]
	returns the automorphism group order of the complex (given as a list facetsLst of facets).";
	
SimplicialComplexAutomorphismGroupOrder::argerr="Input should be of the form 
	SimplicialComplexAutomorphismGroupOrder[facetsLst_List].";


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroup*)


(* :Usage Messages: *)

PureComplexAutomorphismGroup::usage="PureComplexAutomorphismGroup[facetsLst_List] returns 
	a list {a,b} where a is the automorphism group of the canonically relabeled version of the
	input complex (given as a list facetsLst of facets), and b is the vertex relabeling rule
	from facetsLst to the canonically labeled version. The canonical relabeling sends the vertices
	to 1..n, so a acts on those labels; b maps them back. a is given by generators, not by its
	elements.";
	
PureComplexAutomorphismGroup::argerr="Input should be of the form 
	PureComplexAutomorphismGroup[facetsLst_List].";


(* ::Item::Closed:: *)
(*PureComplexAutomorphismGroupOrder*)


(* :Usage Messages: *)

PureComplexAutomorphismGroupOrder::usage="PureComplexAutomorphismGroupOrder[facetsLst_List] 
	computes the automorphism group order of the pure complex facetsLst.";
	
PureComplexAutomorphismGroupOrder::argerr="Input should be of the form 
	PureComplexAutomorphismGroupOrder[facetsLst_List].";


(* ::Item::Closed:: *)
(*PureComplexFacetStabilizerGroupOrder*)


(* :Usage Messages: *)

PureComplexFacetStabilizerGroupOrder::usage="PureComplexFacetStabilizerGroupOrder[facetsLst_List]
	computes the order of the facet stabilizer group of the pure complex facetsLst.";

PureComplexFacetStabilizerGroupOrder::argerr="Input should be of the form 
	PureComplexFacetStabilizerGroupOrder[facetsLst_List].";


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroupOrder*)


(* :Usage Messages: *)

PureComplexFacetAutomorphismGroupOrder::usage="
	PureComplexFacetAutomorphismGroupOrder[facetsLst_List] computes the order of the 
	facet automorphism group of the pure complex given as a list facetsLst of facets.";
	
PureComplexFacetAutomorphismGroupOrder::argerr="Input should be of the form 
	PureComplexFacetAutomorphismGroupOrder[facetsLst_List].";


(* ::Item::Closed:: *)
(*PureComplexFacetAutomorphismGroup*)


(* :Usage Messages: *)

PureComplexFacetAutomorphismGroup::usage="PureComplexFacetAutomorphismGroup[facetsLst_List]
	returns the facet automorphism group of the pure simplicial complex given as a list facetsLst
	of facets, that is the image of its automorphism group in the symmetric group on the facets.
	The permutations returned act on the FACET INDICES 1..Length[facetsLst], where index j means
	facetsLst[[j]] -- the position in the list as supplied, NOT the lexicographic rank of the
	facet. For {{1,2,3},{1,2,4}} facet 1 is {1,2,3}; for {{1,2,4},{1,2,3}} facet 1 is {1,2,4}.
	The group is given by generators rather than by its elements, which is often the only
	possibility: a 20-facet complex with pairwise interchangeable facets has the full symmetric
	group on 20 points here, of order 2432902008176640000.";
	
PureComplexFacetAutomorphismGroup::argerr="Input should be of the form 
	PureComplexFacetAutomorphismGroup[facetsLst_List].";


(* ::Subsection::Closed:: *)
(*Clique Complex Automorphism and Facet Automorphism*)


(* ::Item::Closed:: *)
(*CliqueFacetStabilixerGroupOrder*)


(* :Usage Messages: *)

CliqueFacetStabilizerGroupOrder::usage="CliqueFacetStabilizerGroupOrder[facetsLst_List] 
	returns the order of the facet stabilizer group of the clique complex given as a list facetsLst
	of maximal cliques. It does not check whether or not the complex is a clique complex, it assumes
	it is.";

CliqueFacetStabilizerGroupOrder::argerr="Input should be of the form 
	CliqueFacetStabilizerGroupOrder[facetsLst_List].";


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroupOrder*)


(* :Usage Messages: *)

CliqueFacetAutomorphismGroupOrder::usage="CliqueFacetAutomorphismGroupOrder[facetsLst_List] 
	returns the facet automorphism group order of the clique complex given as a list facetsLst of 
	maximal cliques. It does not check whether or not the input complex is a clique complex.";
	
CliqueFacetAutomorphismGroupOrder::argerr="Input should be of the form 
	CliqueFacetAutomorphismGroupOrder[facetsLst_List].";


(* ::Item::Closed:: *)
(*CliqueFacetAutomorphismGroup*)


(* :Usage Messages: *)

CliqueFacetAutomorphismGroup::usage="CliqueFacetAutomorphismGroup[facetsLst_List] returns the
	facet automorphism group of the clique complex given as a list facetsLst of maximal cliques,
	that is the image of its automorphism group in the symmetric group on the facets. It does not
	check whether or not the complex is a clique complex.
	The permutations returned act on the FACET INDICES 1..Length[facetsLst], where index j means
	facetsLst[[j]] -- the position in the list as supplied, NOT the lexicographic rank of the
	facet. For {{1,2,3},{1,2,4}} facet 1 is {1,2,3}; for {{1,2,4},{1,2,3}} facet 1 is {1,2,4}.
	The group is given by generators rather than by its elements, which is often the only
	possibility: a 20-facet complex with pairwise interchangeable facets has the full symmetric
	group on 20 points here, of order 2432902008176640000.";
	
CliqueFacetAutomorphismGroup::argerr="Input should be of the form 
	CliqueFacetAutomorphismGroup[facetsLst_List].";


(* ::Chapter:: *)
(*Generate Pure Simplicial Complexes*)


(* ::Section::Closed:: *)
(*Generate All Pure Simplicial Complexes*)


(* ::Subsection::Closed:: *)
(*Generate all vertex labeled pure simplicial complexes *)


(* ::Item::Closed:: *)
(*GenerateAllVertexLabeledPureSimplicialComplexes*)


(* :Usage Messages: *)

GenerateAllVertexLabeledPureSimplicialComplexes::usage="
	GenerateAllVertexLabeledPureSimplicialComplexes[{p_Integer, q_Integer, n_Integer}]; overload
	GenerateAllVertexLabeledPureSimplicialComplexes[p_Integer, q_Integer] returns a 
	list of all vertex-labeled p-pure simplicial complexes of facet-order q. If the number of 
	vertices (n) is specified, it generates all p-pure complexes with q facets and n vertices. 
	Otherwise it generates all p-pure complexes with q facets and all vertex counts between nmin
	and n=p*q. Vertices are labeled from the set [n]={1,2,...,n}.";

(* :Error Messages: *)

GenerateAllVertexLabeledPureSimplicialComplexes::argerr="Input should be of the form 
	GenerateAllVertexLabeledPureSimplicialComplexes[{p_Integer, q_Integer, n_Integer}] or 
	GenerateAllVertexLabeledPureSimplicialComplexes[{p_Integer, q_Integer}].";


(* ::Subsection::Closed:: *)
(*Generate all facet-labeled pure simplicial complexes*)


(* ::Item::Closed:: *)
(*GenerateAllFacetLabeledPureSimplicialComplexes*)


(* :Usage Messages: *)

GenerateAllFacetLabeledPureSimplicialComplexes::usage="
	GenerateAllFacetLabeledPureSimplicialComplexes[{p_Integer,q_Integer}] returns a list 
	of all facet-labeled p-pure simplicial complexes of facet-order q, where each 
	simplicial complex is represented as a list of its facets.";

(* :Error Messages: *)

GenerateAllFacetLabeledPureSimplicialComplexes::argerr="Input should be of the form  
	GenerateAllFacetLabeledPureSimplicialComplexes[{p_Integer, q_Integer}].";


(* ::Section:: *)
(*Generate A Random Pure Simplicial Complex*)


(* ::Subsection::Closed:: *)
(*Random [purity*facetOrder]-labeled Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomPQLabeledPureSimplicialComplex*)


(* :Usage Messages: *)
 
RandomPQLabeledPureSimplicialComplex::usage="
	RandomPQLabeledPureSimplicialComplex[{p_Integer,q_Integer}]; 
	generates a random simplicial complex with purity p and facet-order q, whose vertices are
	labeled from the set [pq]={1,2,...p*q}.

	RandomPQLabeledPureSimplicialComplex[{p_Integer,q_Integer}, numSamples_Integer] generates
	a list with numSamples random simplicial complex with purity p and facet-order q, whose 
	vertices are labeled from the set [pq]={1,2,...p*q}
	";

(* :Error Messages: *)

RandomPQLabeledPureSimplicialComplex::argerr="Input should be of the form
	RandomPQLabeledPureSimplicialComplex[{p_Integer,q_Integer}] or 
	RandomPQLabeledPureSimplicialComplex[{p_Integer,q_Integer}, numSamples_Integer].";


(* ::Subsection::Closed:: *)
(*Random Vertex labeled Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomVertexLabeledPureSimplicialComplex*)


(* :Usage Messages: *)
 
RandomVertexLabeledPureSimplicialComplex::usage="
	RandomVertexLabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}, numSamples_Integer]; overload 
	RandomVertexLabeledPureSimplicialComplex[{p_Integer,q_Integer}, numSamples_Integer]
	generates numSamples random vertex-labeled pure simplicial complexes with purity p, 
	facet-order q, and vertex count n (if specified). If numSamples is not specified it returns
	a single random complex. ";

(* :Error Messages: *)	

RandomVertexLabeledPureSimplicialComplex::argerr="Input should be of the form
	RandomVertexLabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}, numSamples_Integer] or
	RandomVertexLabeledPureSimplicialComplex[{p_Integer,q_Integer}, numSamples_Integer], or
	RandomVertexLabeledPureSimplicialComplex[{p_Integer,q_Integer}], or
	RandomVertexLabeledPureSimplicialComplex[{p_Integer,q_Integer}].";

RandomVertexLabeledPureSimplicialComplex::empty="There is no pure simplicial complex of purity
	`1`, facet order `2` and vertex count `3`, so there is nothing to sample. Purity and facet
	order must be at least 1, and the vertex count must satisfy p <= n <= p q with at least q
	distinct p-subsets available.";


(* ::Subsection::Closed:: *)
(*Random Uniform Unlabeled Pure Simplicial Complex by rejection from vertex-labeled*)


(* ::Item::Closed:: *)
(*RandomUniformUnlabeledPureSimplicialComplex*)


(* :Usage Messages: *)
 
RandomUniformUnlabeledPureSimplicialComplex::usage="
	RandomUniformUnlabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}, numSamples_Integer], overload;
	RandomUniformUnlabeledPureSimplicialComplex[{p_Integer,q_Integer}, numSamples_Integer] generates
	numSamples random unlabeled pure simplicial complexes with purity p,
	facet-order q, and if specified vertex count n. If numSamples is not specified it returns
	a single random complex.
Neither the vertices nor the facets are labelled, so these are the isomorphism classes: the
distribution is uniform over the classes counted by NumUnlabeledPureComplexes, and is produced
without rejection. With the vertex count free it is drawn from NumUnlabeledPureComplexes too.";

(* :Error Messages: *)	

RandomUniformUnlabeledPureSimplicialComplex::argerr="Input should be of the form
	RandomUniformUnlabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}, numSamples_Integer] or
	RandomUniformUnlabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}].";

RandomUniformUnlabeledPureSimplicialComplex::empty="There is no pure simplicial complex of purity
	`1`, facet order `2` and vertex count `3`, so there is nothing to sample. Purity and facet
	order must be at least 1, and the vertex count must satisfy p <= n <= p q with at least q
	distinct p-subsets available.";


(* ::Subsection::Closed:: *)
(*Random Uniform Facet-labeled Pure Simplicial Complex by Burnside pair sampling*)


(* ::Item::Closed:: *)
(*RandomUniformFacetLabeledPureSimplicialComplex*)


(* :Usage Messages: *)
 
RandomUniformFacetLabeledPureSimplicialComplex::usage="
	RandomUniformFacetLabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}, numSamples_Integer], overload;
	RandomUniformFacetLabeledPureSimplicialComplex[{p_Integer,q_Integer}, numSamples_Integer] generates numSamples random FACET-LABELED pure simplicial complexes with purity p,
	facet-order q, and vertex count n if specified -- q labeled, pairwise distinct p-subsets of an
	n-element vertex set covering it, with the vertices unlabeled. If numSamples is not specified
	it returns a single random complex.
	The distribution is exactly uniform over those complexes, the number of which is
	NumFacetLabeledPureComplexes[p, q, n], and is produced without rejection.";

(* :Error Messages: *)	

RandomUniformFacetLabeledPureSimplicialComplex::argerr="Input should be of the form
	RandomUniformFacetLabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}, numSamples_Integer] or
	RandomUniformFacetLabeledPureSimplicialComplex[{p_Integer,q_Integer,n_Integer}].";

RandomUniformFacetLabeledPureSimplicialComplex::empty="There is no pure simplicial complex of
	purity `1`, facet order `2` and vertex count `3`, so there is nothing to sample. Purity and
	facet order must be at least 1, and the vertex count must satisfy p <= n <= p q with at
	least q distinct p-subsets available.";


(* ::Subsection::Closed:: *)
(*Random facet labeled Pure Simplicial Complex ( not uniform)*)


(* ::Item:: *)
(*RandomFacetLabeledPureSimplicialComplex*)


(* :Usage Messages: *)

RandomFacetLabeledPureSimplicialComplex::usage="
	RandomFacetLabeledPureSimplicialComplex[{p_Integer, q_Integer}]; overload
	RandomFacetLabeledPureSimplicialComplex[{p_Integer, q_Integer}, numSamples_Integer] 
	generates a random 
	facet-labeled simplicial complex with purity p and facet-order q.";
	
(* :Error Messages: *)

RandomFacetLabeledPureSimplicialComplex::argerr="Integers were expected at positions 1 and 2.";


(* ::Subsection:: *)
(*MCMC Random Pure Simplicial Complex*)


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexIterMCMCSweep*)


(* :Usage Messages: *)

RandomPureSimplicialComplexMCMCSweep::usage="
	RandomPureSimplicialComplexMCMCSweep[seedComplex_List,labelingChoise_Integer,NN_Integer],
	overload RandomPureSimplicialComplexMCMCSweep[seedComplex_List,HoldNumVerticesFixed_Boole,labelingChoise_Integer,NN_Integer] 
	performs NN Markov Chain Monte Carlo sweeps (a sweep is purity*facetOrder steps) on a 
	seed pure simplicial complex seedComplex. labelingChoise is an integer in {0,1,2}. 0 = 
	vertex labeled, 1 = facet-labeled 2 = unlabeled. In the overload, HoldNumVerticesFixed=True 
	keeps the total numbe rof vertices fixed, otherwise it defaults to the primary function
	which does not hold number of vertices fixed";
	
(* :Error Messages: *)

RandomPureSimplicialComplexMCMCSweep::argerr="Input has to be of the form 
	RandomPureSimplicialComplexMCMCSweep[seedComplex_List,labelingChoise_Integer,NN_Integer].";


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexMCMCEquilibriate*)


(* :Usage Messages: *)

RandomPureSimplicialComplexMCMCEquilibriate::usage="
	RandomPureSimplicialComplexMCMCEquilibriate[seedComplex_List,labelingChoise_Integer] 
	runs enough MC sweeps to equilibriate the input seedComplex. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled. Equilibriation 
	is established when four different complexes (the input, the complex with the minumum 
	number of vertices, the complex with the maximum number of vertices, and another random 
	complex) all start to have close enough probabilities such that the fluctions of each is
	of the order of the fluctuation across them.";
	
(* :Error Messages: *)

RandomPureSimplicialComplexMCMCEquilibriate::argerr="Input has to be of the form 
	RandomPureSimplicialComplexMCMCEquilibriate[seedComplex_List,labelingChoise_Integer].";

(* Concatenated for the same reason as GraphEquilibriate::noconv above. *)
RandomPureSimplicialComplexMCMCEquilibriate::noconv=
	"Equilibriation did not converge within the budget of `1` sweeps. At exit the mean "<>
	"squared pairwise difference in `2` across the tracks -- the observable that agreed "<>
	"least well -- was `3`, against `4`, the value that quantity is expected to take once "<>
	"the tracks are equilibriated; convergence requires the former to be below the latter. "<>
	"Returning eqlT -> `5`, which is a lower bound rather than an estimate, and a state that "<>
	"may not be equilibriated. Raise $ECGravMaxEquilibriationSweeps, or start from a better "<>
	"seed complex.";

(* Concatenated for the same reason as ::noconv above. *)
RandomPureSimplicialComplexMCMCEquilibriate::nosignal=
	"None of the monitored observables (`1`) took more than one value over the whole "<>
	"comparison window, on any track, so there was nothing for the convergence test to "<>
	"read and equilibriation could not be checked. Returning eqlT -> `2` and converged -> "<>
	"True, since a chain whose observables never move is stationary in them, but note that "<>
	"this is the absence of evidence rather than evidence of mixing: if the chain is stuck, "<>
	"this is what that looks like. Pass observables that distinguish the states you care "<>
	"about, or check the returned state by hand.";


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexMCMCCorrelationTime*)


(* :Usage Messages: *)

RandomPureSimplicialComplexMCMCCorrelationTime::usage="
	RandomPureSimplicialComplexMCMCCorrelationTime[equilibriatedComplex_List,
	eqlT_Integer, operators_,labelingChoise_Integer] 
	takes as an input equilibriated complex, the equlibriation time, a list of functions, and 
	a choice of labeleing and determines the correltion time. Correlation time is the maximum 
	of the correlation times for energy = -Log[prob], the operator which counts the number
	of edges, and the input list of operators. labelingChoise is an 
	integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled.";
	
(* :Error Messages: *)

RandomPureSimplicialComplexMCMCCorrelationTime::argerr="Input has to be of the form
	RandomPureSimplicialComplexMCMCEquilibriate[equilibriatedComplex_List,
	eqlT_Integer, operators_,labelingChoise_Integer].";

RandomPureSimplicialComplexMCMCCorrelationTime::stuck="All measured values of `1` are
	stuck at `2`. Autocorrelation time can not be computed; a default value of 2 will
	be assigned.";

RandomPureSimplicialComplexMCMCCorrelationTime::alldefault="Autocorrelation time could
	not be computed since all measured values of the observables are stuck at the
	values `1`. Exiting with default values corrT -> 2, corrTValues -> `2`.";


(* ::Item::Closed:: *)
(*RandomPureSimplicialComplexMCMC*)


(* :Usage Messages: *)

RandomPureSimplicialComplexMCMC::usage="
	RandomPureSimplicialComplexMCMC[seedComplex_List, operators_,
	labelingChoise_Integer, numSamples_Integer] 
	takes as an input seed complex, a list of functions, a choice of labeleing and 
	number of samples and outputs a list of numSamples independent random pure simplicial 
	complexes from a uniform disribution. it first runs equilibriation run and determines correlation times
	labelingChoise is an integer in {0,1,2}. 0 = vertex labeled, 1 = facet-labeled 2 = unlabeled.";
	
(* :Error Messages: *)

RandomPureSimplicialComplexMCMC::argerr="Input has to be of the form 
	RandomPureSimplicialComplexMCMC[seedComplex_List, operators_,
	labelingChoise_Integer, numSamples_Integer].";

(* Concatenated for the same reason as the ::noconv messages. *)
RandomPureSimplicialComplexMCMC::noneql=
	"Equilibriation did not converge, so the chain has not been shown to have reached its "<>
	"stationary distribution and samples drawn from it now would not be the uniform draws "<>
	"this function promises. Aborting instead of sampling, and returning $Failed. See the "<>
	"RandomPureSimplicialComplexMCMCEquilibriate::noconv message just issued, then raise "<>
	"$ECGravMaxEquilibriationSweeps (currently `1`) or start from a better seed complex.";


(* ::Section::Closed:: *)
(*Generate A Random Pseudo-manifold*)


(* ::Subsection::Closed:: *)
(*Random Pseudo Manifold Through Successive Facet Addition*)


(* ::Item::Closed:: *)
(*RandomUnlabeledPseudoManifold*)


(* :Usage Messages: *)

RandomUnlabeledPseudoManifold::usage="
	RandomUnlabeledPseudoManifold[{p_Integer, q_Integer}];  
	generates a random ulabeled pseudo-manifold with purity p and facet-order q.";
	
(* :Error Messages: *)

RandomUnlabeledPseudoManifold::argerr="Integers were expected at positions 1 and 2.";


(* ::Subsection::Closed:: *)
(*Add a random facet to a PseudoManifold*)


(* ::Item:: *)
(*AddRandomUnlabeledFacetToPseudoManifold*)


(* :Usage Messages: *)

AddRandomUnlabeledFacetToPseudoManifold::usage="
	AddRandomUnlabeledFacetToPseudoManifold[facetsLst_List,apastingSites_List];  
	adds a random facet to the existing pseudo-manifold given as the list
	facetsLst of facets at pasting sites apastingSites which are co-dimension 1 faces 
	available for binding. The method picks out one representative from the orbit of the 
	binding sites under the automorphism group of the existing pseudo-manifold.";
	
(* :Error Messages: *)

AddRandomUnlabeledFacetToPseudoManifold::argerr="Input should be of the form
	AddRandomUnlabeledFacetToPseudoManifold[facetsLst_List,apastingSites_List].";


(* ::Title:: *)
(*ECGrav Private*)


Begin["`Private`"]; (* Begin private context *)


(*Needs["ECGrav`MCSims`"];*)


(* Load the subpackages *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "Subpackages", "MCSims.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "Subpackages", "PureComplexes.wl"}]];


(* ::Title:: *)
(*ECGrav Protect and End*)


End[] (* End private context *)

(* Protect any symbols you want to prevent modifications to *)
(*Protect[
    PublicFunction1,
    PublicFunction2,
    MCSimsFunction,
    PureComplexesFunction
];*)
Protect @@ Names["ECGrav`*"];

(* The equilibriation budget is a user setting, so it has to stay assignable after the
	blanket Protect above. *)
Unprotect[$ECGravMaxEquilibriationSweeps];

EndPackage[]
