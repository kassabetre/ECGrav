#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*TemperingPlots -- MBAR surfaces from GraphCTLSchedule / GraphParallelTempering runs*)

(*
	Load with

		Needs["ECGrav`"];
		Get["/path/to/ECGrav/TemperingPlots.wls"];

	Requires ECGrav 1.11.0 or later.  The whole pipeline is built on the private MBAR helpers
	(MBARWeightBasis and friends) introduced there.  On an older paclet
	ECGrav`Private`MBARWeightBasis has no definition, returns unevaluated, and every surface
	silently inherits the unevaluated head -- CTLSurfaceData checks for this and says so.

	WHAT IT DOES

	A run measures its conjugate observables at K bootstrap field points.  MBAR turns that into a
	smooth function of the target field, and HomogeneousHamiltonian.md sections 7.1-7.4 give the
	recipe one point at a time.  Evaluating a whole plotting grid that way rebuilds the MBAR basis
	at every point and re-walks all S samples once per observable per point.

	Everything here is computed in ONE batched pass instead:

	  * the basis is built once (section 7.1's advice, worth 5-10x on its own);
	  * the weights for the entire grid are one T x S matrix product, not T separate S-vectors;
	  * every surface -- means, variances, susceptibilities, Binder cumulants, the free energy,
	    the internal energy, the heat capacity -- is one matrix product against those weights,
	    because MBAR weights do not depend on the observable.

	The cost is therefore essentially independent of how many surfaces you ask for.

	COORDINATES.  In homogeneous c.O form the sampled field vector is c = beta {1, nt1, nt2, ...}
	and the physical hamiltonian is H_phys = O.v with v = {1, nt1, nt2, ...}.  Plot coordinates
	here are always the physical {beta, nt1, ...}; c is formed internally.  Plotting in c instead
	turns the bootstrap rectangle into a fan through the origin (section 7.1).

	UNITS.  The stored beta is 1.0 and the recorded energy in chart column 3 is c.O = beta*H_phys,
	so a naive Mean of that column is beta*<H_phys> and the textbook beta^2*Variance of it carries
	beta^2 too many.  Every surface here, and every bootstrap value, is in physical units.
*)

(*ECGrav` must be a declared dependency, not merely loaded: BeginPackage with one argument sets
	$ContextPath to {"ECGravPlots`", "System`"} and hides it, so ComputeMinusBetaTimesFreeEnergy
	would resolve to a fresh undefined symbol in this package's private context and return
	unevaluated -- silently, and with the same signature as calling a private helper that a too-old
	paclet never defined.*)
BeginPackage["ECGravPlots`", {"ECGrav`"}];

CTLSurfaceData::usage =
"CTLSurfaceData[res] bundles the return of GraphCTLSchedule (or of GraphMultiHistogram) into the \
association every other function here takes.  res may be the {schedule, hist} pair or hist itself. \
Option \"ObservableNames\" names the entries of chart column 6, in the order the obs list was \
given to the driver.";

BetaSurfaceData::usage =
"BetaSurfaceData[res] bundles a BETA-tempering run -- GraphParallelTempering or \
GraphMultiHistogram over a beta ladder -- into the same association the external-field entry \
points return, so every surface here works on it.  Beta tempering is the one-component \
homogeneous form, c = {beta} and O = {energy}, and this does that mapping.  Non-numeric \
observable columns are dropped with a message; \"ObservableNames\" names the ones that survive.";

TemperingSurfaceData::usage =
"TemperingSurfaceData[res] does the same for the return of GraphParallelTempering, whose chart \
carries no free energies; they are solved with ComputeMinusBetaTimesFreeEnergy.";

AddObservable::usage =
"AddObservable[data, name, f] registers a derived per-sample quantity.  f is applied to the \
association of already-registered flat sample vectors, e.g. \
AddObservable[data, \"mag\", #numEdges/78 &].  A plain length-S vector is accepted too.";

SurfaceGrid::usage =
"SurfaceGrid[{bMin, bMax, nB}, {ntMin, ntMax, nNt}] builds the regular (beta, nt0) grid the \
surfaces are evaluated on.  SurfaceGrid[points] wraps an explicit list of coordinates.";

MBARSurfaces::usage =
"MBARSurfaces[data, grid] evaluates every surface on grid in one batched pass and returns an \
association of names to value lists, one value per grid point.  Options: \"ESSFloor\" (Automatic \
= half a replica's worth of samples, None to disable masking), \"SiteCounts\" (association \
name -> NN for the per-site susceptibilities), \"ChunkSize\".";

BootstrapValues::usage =
"BootstrapValues[data] gives the same quantities at the K bootstrap field points from their own \
samples, with no reweighting -- the reference the surfaces must pass through.  Returns an \
association of names to lists of {beta, nt0, ..., value}.";

SurfacePlots::usage =
"SurfacePlots[data, grid, surfaces] returns an association of names to ListPlot3D of each surface \
with the bootstrap points overlaid, mirroring Show[contplot, listplot].";

SurfacePoints::usage =
"SurfacePoints[grid, values] pairs a grid with a value list into plottable triples, dropping \
masked (Indeterminate) points.";

SurfaceNames::usage =
"SurfaceNames[data] lists the surface names MBARSurfaces produces for this data bundle.";

MBARHelpersAvailableQ::usage =
"MBARHelpersAvailableQ[] is True when the loaded ECGrav paclet has the private MBAR helpers this \
pipeline needs (1.11.0 or later).  Testing the definitions is the reliable check: a version string \
compared with Less does not even evaluate.";

SwapAcceptance::usage =
"SwapAcceptance[res] gives swapAccept/swapTry per replica from a run's histories.";

PhysicalToField::usage =
"PhysicalToField[{beta, nt1, ...}] gives the sampled field vector c = beta {1, nt1, ...}.";

FieldToPhysical::usage =
"FieldToPhysical[c] inverts PhysicalToField.";

Begin["`Private`"];

(* ::Section:: *)
(*Coordinates*)

(*c = beta {1, nt1, ...} and H_phys = O.{1, nt1, ...}.  Both directions are needed: the plot axes
	are physical, the chart keys are c.*)
PhysicalToField[p_List] := N[p[[1]]*Prepend[Rest[p], 1.]];
FieldToPhysical[c_List] := N[Prepend[Rest[c]/c[[1]], c[[1]]]];
hamiltonianDirection[p_List] := N[Prepend[Rest[p], 1.]];

(*A single-field run keys its chart on bare reals rather than one-element lists -- MBARWeightBasis
	coerces them and the nested-Sum form it replaced tolerated it through Listable broadcasting.
	Accept both here so nothing downstream has to care which shape a run happens to use.*)
PhysicalToField[b_?NumericQ] := {N[b]};
FieldToPhysical[c_?NumericQ] := {N[c]};
hamiltonianDirection[b_?NumericQ] := {1.};
asFieldVectors[keys_List] := If[MatrixQ[keys], N[keys], List /@ N[keys]];

(* ::Section:: *)
(*Extraction*)

(* The multi-histogram return is {groundStates, chart, replicas, histories}; the chart is
	hist[[2]], keyed by the field vector c, and its "data" rows carry the six columns listed in
	HomogeneousHamiltonian.md 7.1.  A GraphCTLSchedule return wraps that as {schedule, hist}. *)
histOf[res_List] := Which[
	Length[res] === 2 && ListQ[res[[2]]] && Length[res[[2]]] === 4, res[[2]],
	Length[res] === 4 && AssociationQ[res[[2]]], res,
	True, $Failed
];

(* True/False observables (PureGraphQ and friends) have to become numbers before anything can be
	averaged. *)
numericColumn[v_List] := If[VectorQ[v, NumericQ], N[v], N[Boole[v]]];

CTLSurfaceData::obscount = "`1` observable names given for `2` columns in chart column 6.";
CTLSurfaceData::badres =
"Expected the {schedule, hist} return of GraphCTLSchedule or the 4-element multi-histogram itself.";
(*No backticks in the template: StringForm reads the pair around Private in a context path as a
	named slot and then complains about the numbered one.*)
CTLSurfaceData::nohelpers =
"The private MBAR weight-basis helper has no definition.  This pipeline needs ECGrav 1.11.0 or \
later; the loaded paclet is `1`.  Install the current paclet (or PacletDirectoryLoad the source \
directory) and restart the kernel.";

MBARHelpersAvailableQ[] := Length[DownValues[ECGrav`Private`MBARWeightBasis]] > 0;

checkHelpers[] :=
	If[Length[DownValues[ECGrav`Private`MBARWeightBasis]] === 0,
		Message[CTLSurfaceData::nohelpers,
			Quiet@Check[PacletObject["ECGrav"]["Version"], "unknown"]]; False,
		True];

(*bt must be Real: MBARWeightBasis takes ?NumericQ, but the public entry points pattern-match
	betaFixed_Real and an integer 1 falls through to their catch-all.*)
buildData[chart_Association, bt_, obsNames_] :=
Module[{keys, mbf, conj, om, nVars, counts, vals, nObs, col6},
	keys = Keys[chart];
	mbf = chart[[All, "minusBetaF"]];
	conj = chart[[All, "data", All, 4]];
	counts = Length /@ Lookup[conj, Key /@ keys];

	(*Flatten in the basis's key order -- Lookup[..., Key/@keys], never Values.  The basis stacks
		its samples that way and a mismatch is silent (HomogeneousHamiltonian.md 7.1).*)
	om = N[Join @@ Lookup[conj, Key /@ keys]];
	If[!MatrixQ[om], om = Transpose[{om}]];
	nVars = Last[Dimensions[om]];

	col6 = Join @@ Lookup[chart[[All, "data", All, 6]], Key /@ keys];
	nObs = If[Length[col6] > 0, Length[First[col6]], 0];
	If[ListQ[obsNames] && Length[obsNames] =!= nObs,
		Message[CTLSurfaceData::obscount, Length[obsNames], nObs]];

	vals = Join[
		Association@Table["O" <> ToString[i] -> om[[All, i]], {i, nVars}],
		Association@Table[
			If[ListQ[obsNames] && i <= Length[obsNames], obsNames[[i]], "obs" <> ToString[i]]
				-> numericColumn[col6[[All, i]]],
			{i, nObs}]];

	<|
		"beta" -> N[bt], "keys" -> asFieldVectors[keys], "chartKeys" -> keys, "counts" -> counts,
		"minusBetaF" -> mbf, "conj" -> conj,
		"nVars" -> nVars, "nSamples" -> Length[om],
		"obs" -> om, "values" -> vals,
		"chart" -> chart,
		"basis" -> ECGrav`Private`MBARWeightBasis[N[bt], mbf, conj]
	|>
];

Options[CTLSurfaceData] = {"ObservableNames" -> Automatic};

CTLSurfaceData[res_List, OptionsPattern[]] :=
Module[{hist = histOf[res], chart, bt},
	If[hist === $Failed, Message[CTLSurfaceData::badres]; Return[$Failed]];
	If[!checkHelpers[], Return[$Failed]];
	chart = hist[[2]];
	bt = N@First[Values[chart[[All, "beta"]]]];
	buildData[chart, bt, OptionValue["ObservableNames"]]
];

Options[TemperingSurfaceData] = {"ObservableNames" -> Automatic, "Beta" -> Automatic};

TemperingSurfaceData::badres =
"Expected the 4-element {groundStates, chart, replicas, histories} return of GraphParallelTempering.";

(*GraphParallelTempering's chart uses the SAME six columns but carries no free energies, so they
	are solved here -- sub-second since 1.12.0.*)
TemperingSurfaceData[res_List, OptionsPattern[]] :=
Module[{chart, bt, mbf, wrapped},
	If[Length[res] =!= 4 || !AssociationQ[res[[2]]],
		Message[TemperingSurfaceData::badres]; Return[$Failed]];
	If[!checkHelpers[], Return[$Failed]];
	chart = res[[2]];
	bt = N[OptionValue["Beta"] /. Automatic -> 1.0];
	mbf = ComputeMinusBetaTimesFreeEnergy[chart[[All, All, 4]], bt];
	wrapped = Association@Table[
		k -> <|"minusBetaF" -> mbf[[Key[k]]], "beta" -> bt, "externalField" -> k,
			"data" -> chart[[Key[k]]]|>,
		{k, Keys[chart]}];
	buildData[wrapped, bt, OptionValue["ObservableNames"]]
];

(*Beta tempering IS the one-component homogeneous form: the field is c = {beta}, the conjugate
	observable is the energy, and betaFixed = 1, so c.O = beta E.  Everything in this file then
	applies unchanged -- but the beta chart does not LOOK like the external-field one, and reading
	it as if it did is silently wrong: its column 4 is the observable list, not the conjugate
	observables, and it has no residual column at all.

	Two row layouts occur.  Since 1.13.0 a beta row is {sweep, beta, E, {obs...}, tag}; before that
	the observables were spliced in, giving {sweep, beta, E, obs...}.  Both are read.

	The chart is re-keyed on {beta}, a one-element vector rather than the bare Real the drivers use,
	because the coordinate helpers here take lists -- FieldToPhysical[c_List] does not match a Real.

	The free energies are always re-solved rather than taken from the run, because a run made before
	1.12.0 stored under-converged ones (HomogeneousHamiltonian.md, and CHANGELOG 1.12.0). It costs
	well under a second.*)

BetaSurfaceData::badres =
"Expected a beta-tempering return: the {schedule, hist} pair, or a driver return whose second \
element is the chart.";
BetaSurfaceData::dropped =
"Dropped `1` observable column(s) that are neither numeric nor True/False -- positions `2`. An \
observable returning a graph is the usual cause; it cannot be averaged.";
BetaSurfaceData::noobs = "No usable observable columns remain.";

Options[BetaSurfaceData] = {"ObservableNames" -> Automatic};

BetaSurfaceData[res_List, OptionsPattern[]] :=
Module[{h, chart, rows, keys, first, nested, obsOf, probe, nObs, keep, drop, names, given,
		reshaped},

	h = Which[
		Length[res] === 2 && ListQ[res[[2]]] && Length[res[[2]]] >= 4, res[[2]],
		Length[res] >= 2 && AssociationQ[res[[2]]], res,
		True, $Failed];
	If[h === $Failed, Message[BetaSurfaceData::badres]; Return[$Failed]];
	If[!checkHelpers[], Return[$Failed]];

	chart = h[[2]];
	rows = If[AssociationQ[First[Values[chart]]], chart[[All, "data"]], chart];
	keys = Keys[rows];
	first = First[First[Values[rows]]];

	(*Nested since 1.13.0, spliced before it.*)
	nested = ListQ[first[[4]]];
	obsOf[row_] := If[nested, row[[4]], row[[4 ;;]]];
	nObs = Length[obsOf[first]];

	(*Keep only what can be averaged, judged on the first row of every replica so a column that is
		a graph in one replica and a number in another cannot slip through.*)
	probe = Transpose[obsOf[First[#]] & /@ Values[rows]];
	keep = Select[Range[nObs], AllTrue[probe[[#]], NumericQ[#] || BooleanQ[#] &] &];
	drop = Complement[Range[nObs], keep];
	If[drop =!= {}, Message[BetaSurfaceData::dropped, Length[drop], drop]];
	If[keep === {}, Message[BetaSurfaceData::noobs]; Return[$Failed]];

	given = OptionValue["ObservableNames"];
	names = If[ListQ[given] && Length[given] === nObs, given[[keep]], given];

	reshaped = Association @ Table[
		With[{b = N[If[ListQ[k], First[k], k]]},
			{b} -> (With[{row = #}, {row[[1]], {b}, b*row[[3]], {N[row[[3]]]}, 0.,
				obsOf[row][[keep]]}] & /@ rows[[Key[k]]])],
		{k, keys}];

	TemperingSurfaceData[{<||>, reshaped, <||>, <||>}, "ObservableNames" -> names]
];

(* ::Section:: *)
(*Derived per-sample quantities*)

AddObservable::len = "Vector for `1` has length `2`, expected `3` samples.";

AddObservable[$Failed, ___] := $Failed;

AddObservable[data_Association, name_String, f_Function] :=
	AddObservable[data, name, f[data["values"]]];

AddObservable[data_Association, name_String, v_List] :=
	If[Length[v] =!= data["nSamples"],
		Message[AddObservable::len, name, Length[v], data["nSamples"]]; data,
		MapAt[Append[#, name -> numericColumn[v]] &, data, {Key["values"]}]];

SurfaceNames[data_Association] :=
	Join[{"ESS", "logZ", "betaF", "FPhys", "EPhys", "EHom", "Cv", "BinderEnergy"},
		Join @@ Table[{n, "Var" <> n, "chi" <> n, "Binder" <> n}, {n, Keys[data["values"]]}]];

(* ::Section:: *)
(*The grid*)

SurfaceGrid[ranges : {_, _, _Integer} ..] :=
	With[{axes = Table[N@Subdivide[r[[1]], r[[2]], r[[3]] - 1], {r, {ranges}}]},
		<|"points" -> Flatten[Outer[List, Sequence @@ axes], Length[axes] - 1],
		  "dimensions" -> ({ranges}[[All, 3]]),
		  "axes" -> axes|>];

SurfaceGrid[pts_List] := <|"points" -> N[pts], "dimensions" -> None, "axes" -> None|>;

(* ::Section:: *)
(*The one batched pass*)

(*The MBAR weight of sample s at target c is Exp[-(beta c.O_s + logBasis_s)] up to a per-target
	constant, and every quantity built from it is a ratio, so the constant cancels.  For T targets
	that is one T x S matrix -- C.Transpose[om] -- instead of T separate S-vectors, and each
	observable is then a single matrix product against it.  The per-row shift is the only thing
	the free energy needs.*)
batchStats[basis_Association, fields_List, valueMat_List] :=
Module[{om = basis["obs"], lb = basis["logBasis"], bt = basis["beta"], logUT, shifts, ut, tot},
	(*Built S x T first, because logBasis is per SAMPLE and Plus threads over the outer dimension:
		in that orientation each entry of logBasis meets a whole row and broadcasts, which it
		cannot do against a T x S matrix.  Transposed once, after which every reduction is
		row-wise over targets.*)
	logUT = Transpose[-(bt*(om . Transpose[N[fields]]) + lb)];
	shifts = Max /@ logUT;
	ut = Exp[Clip[logUT - shifts, {-700., 0.}]];
	tot = Total[ut, {2}];
	<|"logNorm" -> shifts + Log[tot],
	  "ess" -> tot^2/Total[ut^2, {2}],
	  "means" -> If[Length[valueMat] > 0,
			(ut . Transpose[valueMat])/tot,
			ConstantArray[{}, Length[fields]]]|>
];

(*Binder cumulant 1 - <x^4>/(3 <x^2>^2).  A constant-zero observable -- an all-False boolean
	column such as ClosedDGraphQ on a run that never closed -- has <x^2> = 0 and no cumulant; say
	so with Indeterminate rather than letting 0/0 raise Power::infy for every point.*)
SetAttributes[binder, Listable];
binder[fourth_, second_] := If[second == 0., Indeterminate, 1.0 - fourth/(3.0*second^2)];

(*Flattened outer powers of the conjugate observables.  <H^p> for H = O.v is then just
	Flatten[v^(x)p] . <O^(x)p>, with no multi-index bookkeeping; p = 1, 2, 4 covers the internal
	energy, the heat capacity and the energy Binder cumulant.*)
outerPowerColumns[om_List, p_Integer] :=
	Table[Times @@ Table[om[[All, i]], {i, t}], {t, Tuples[Range[Last@Dimensions[om]], p]}];

Options[MBARSurfaces] = {"ESSFloor" -> Automatic, "SiteCounts" -> <||>, "ChunkSize" -> 512};

MBARSurfaces[$Failed, ___] := $Failed;

MBARSurfaces[data_Association, grid_Association, opts : OptionsPattern[]] :=
	MBARSurfaces[data, grid["points"], opts];

MBARSurfaces[data_Association, targets_List, OptionsPattern[]] :=
Module[{om = data["obs"], vals = data["values"], names, valueMat, nV, n1, n2, n4,
		res, means, ess, logNorm, floor, siteCounts, betas, dirs,
		v1, v2, v4, m1, m2, m4, out, mask, fields},

	names = Keys[vals];
	nV = Length[names];
	n1 = data["nVars"]; n2 = n1^2; n4 = n1^4;

	(*One value matrix for the whole run: every registered observable with its square and fourth
		power, then the outer powers of the conjugates.  All of them ride the same weights, which
		is why adding surfaces is nearly free.*)
	valueMat = Join[
		Join @@ Table[{v, v^2, v^4}, {v, Values[vals]}],
		outerPowerColumns[om, 1], outerPowerColumns[om, 2], outerPowerColumns[om, 4]];

	fields = PhysicalToField /@ N[targets];
	res = batchStats[data["basis"], #, valueMat] & /@
			Partition[fields, UpTo[Max[1, OptionValue["ChunkSize"]]]];
	means   = Join @@ res[[All, "means"]];
	ess     = Join @@ res[[All, "ess"]];
	logNorm = Join @@ res[[All, "logNorm"]];

	betas = N[targets][[All, 1]];

	(*v = {1, nt1, ...} is the direction of the physical hamiltonian in observable space.  It moves
		with the plot coordinate; freezing it at the bootstrap nt0 is the trap in
		HomogeneousHamiltonian.md 7.3.*)
	dirs = hamiltonianDirection /@ N[targets];
	v1 = dirs;
	v2 = Flatten /@ (Outer[Times, #, #] & /@ dirs);
	v4 = Flatten /@ (Outer[Times, #, #, #, #] & /@ dirs);

	m1 = MapThread[Dot, {means[[All, 3 nV + 1 ;; 3 nV + n1]], v1}];
	m2 = MapThread[Dot, {means[[All, 3 nV + n1 + 1 ;; 3 nV + n1 + n2]], v2}];
	m4 = MapThread[Dot, {means[[All, 3 nV + n1 + n2 + 1 ;; 3 nV + n1 + n2 + n4]], v4}];

	floor = OptionValue["ESSFloor"] /. Automatic -> 0.5*Mean[N@data["counts"]];
	mask = If[floor === None, ConstantArray[True, Length[targets]], Thread[ess >= floor]];
	siteCounts = OptionValue["SiteCounts"];

	out = <|
		"ESS"          -> ess,
		"logZ"         -> logNorm,
		"betaF"        -> -logNorm,
		"FPhys"        -> -logNorm/betas,
		"EPhys"        -> m1,
		"EHom"         -> betas*m1,
		"Cv"           -> betas^2*(m2 - m1^2),
		"BinderEnergy" -> binder[m4, m2]
	|>;

	Do[
		With[{i = First@FirstPosition[names, nm], nn = N@Lookup[siteCounts, nm, 1]},
			With[{mu = means[[All, 3 i - 2]], sq = means[[All, 3 i - 1]], q = means[[All, 3 i]]},
				AssociateTo[out, {
					nm             -> mu,
					"Var" <> nm    -> sq - mu^2,
					"chi" <> nm    -> betas*nn*(sq - mu^2),
					"Binder" <> nm -> binder[q, sq]}]]],
		{nm, names}];

	(*Indeterminate outside the sampled region, so a plotter leaves a hole instead of drawing
		extrapolation.  ESS itself is never masked -- it is the diagnostic.*)
	If[floor =!= None,
		out = AssociationMap[
			If[# === "ESS", out[#], MapThread[If[#2, #1, Indeterminate] &, {out[#], mask}]] &,
			Keys[out]]];
	out
];

(* ::Section:: *)
(*The bootstrap points*)

Options[BootstrapValues] = {"SiteCounts" -> <||>};

(*The same quantities read straight off each replica's own samples: the reference the continuous
	surfaces must pass through.  Computed in PHYSICAL units -- chart column 3 is c.O = beta*H_phys,
	so Mean of that column carries one factor of beta too many, and
	SpecificHeat[column3, 1, beta] = beta^2 Variance[column3] carries two.*)
BootstrapValues[$Failed, ___] := $Failed;

BootstrapValues[data_Association, OptionsPattern[]] :=
Module[{keys = data["keys"], counts = data["counts"], vals = data["values"], om = data["obs"],
		offs, slices, coords, mbf, eph, siteCounts, trip, out},

	offs = Prepend[Accumulate[counts], 0];
	slices = Table[Span[offs[[i]] + 1, offs[[i + 1]]], {i, Length[keys]}];
	coords = FieldToPhysical /@ keys;
	mbf = Lookup[data["minusBetaF"], Key /@ data["chartKeys"]];
	eph = Table[om[[slices[[i]]]] . hamiltonianDirection[coords[[i]]], {i, Length[keys]}];
	siteCounts = OptionValue["SiteCounts"];

	trip[f_] := Table[Append[coords[[i]], f[i]], {i, Length[keys]}];

	(*No "ESS" here.  A replica's raw sample count is not the MBAR effective sample size, and
		overlaying it on the ESS surface invites reading it as one.*)
	out = <|
		"logZ"         -> trip[mbf[[#]] &],
		"betaF"        -> trip[-mbf[[#]] &],
		"FPhys"        -> trip[-mbf[[#]]/coords[[#, 1]] &],
		"EPhys"        -> trip[Mean[eph[[#]]] &],
		"EHom"         -> trip[coords[[#, 1]]*Mean[eph[[#]]] &],
		"Cv"           -> trip[coords[[#, 1]]^2*(Mean[eph[[#]]^2] - Mean[eph[[#]]]^2) &],
		"BinderEnergy" -> trip[binder[Mean[eph[[#]]^4], Mean[eph[[#]]^2]] &]
	|>;

	Do[
		With[{nn = N@Lookup[siteCounts, nm, 1], vec = vals[nm]},
			With[{xs = Table[vec[[slices[[i]]]], {i, Length[keys]}]},
				AssociateTo[out, {
					nm             -> trip[Mean[xs[[#]]] &],
					"Var" <> nm    -> trip[Mean[xs[[#]]^2] - Mean[xs[[#]]]^2 &],
					"chi" <> nm    -> trip[coords[[#, 1]]*nn*
											(Mean[xs[[#]]^2] - Mean[xs[[#]]]^2) &],
					"Binder" <> nm -> trip[binder[Mean[xs[[#]]^4], Mean[xs[[#]]^2]] &]}]]],
		{nm, Keys[vals]}];

	out
];

(* ::Section:: *)
(*Plots*)

SurfacePoints[grid_Association, values_List] := SurfacePoints[grid["points"], values];

SurfacePoints[targets_List, values_List] :=
	Select[MapThread[Append, {N[targets], values}], NumericQ[Last[#]] &];

(*A regular grid is plotted as an ARRAY rather than as triples.  ListPlot3D leaves a genuine hole
	wherever the array holds Indeterminate, which is the point of the ESS mask; handed triples
	instead it drops the masked points and triangulates straight across the gap, drawing exactly
	the invented surface the mask exists to suppress.  SurfaceGrid's points run beta-major, and
	ListPlot3D reads rows as y, hence the Transpose.*)
surfaceGraphic[grid_Association, values_List, nm_String, popts_List] :=
Module[{dims = grid["dimensions"], axes = grid["axes"]},
	If[ListQ[dims] && Length[dims] === 2,
		ListPlot3D[Transpose[Partition[values, Last[dims]]],
			DataRange -> MinMax /@ axes, PlotRange -> All, PlotLabel -> nm,
			AxesLabel -> {"\[Beta]", "nt0", nm}, Sequence @@ popts],
		ListPlot3D[SurfacePoints[grid, values],
			PlotRange -> All, PlotLabel -> nm,
			AxesLabel -> {"\[Beta]", "nt0", nm}, Sequence @@ popts]]
];

surfaceGraphic[targets_List, values_List, nm_String, popts_List] :=
	ListPlot3D[SurfacePoints[targets, values], PlotRange -> All, PlotLabel -> nm,
		AxesLabel -> {"\[Beta]", "nt0", nm}, Sequence @@ popts];

Options[SurfacePlots] = {"BootstrapPoints" -> Automatic, "PlotOptions" -> {},
	"PointOptions" -> {}, "Names" -> All, "SiteCounts" -> <||>};

(*"SiteCounts" has to be given here as well as to MBARSurfaces, or the overlaid points are drawn
	with NN = 1 against a surface built with the real site count and sit an order of magnitude
	below it -- a mismatch that looks like a reweighting error and is not one.*)
SurfacePlots[_, _, $Failed, ___] := $Failed;
SurfacePlots[$Failed, ___] := $Failed;

SurfacePlots[data_Association, grid_, surfaces_Association, OptionsPattern[]] :=
Module[{boot = OptionValue["BootstrapPoints"], popts = OptionValue["PlotOptions"],
		ptopts = OptionValue["PointOptions"], which = OptionValue["Names"], sel},
	If[boot === Automatic,
		boot = BootstrapValues[data, "SiteCounts" -> OptionValue["SiteCounts"]]];
	sel = If[which === All, Keys[surfaces], Intersection[Keys[surfaces], which]];
	AssociationMap[
		Function[nm,
			Show[
				surfaceGraphic[grid, surfaces[nm], nm, popts],
				If[boot === None || !KeyExistsQ[boot, nm],
					Graphics3D[{}],
					ListPointPlot3D[boot[nm], PlotRange -> All,
						PlotStyle -> Directive[PointSize[0.02], Red], Sequence @@ ptopts]],
				PlotRange -> All]],
		sel]
];

(* ::Section:: *)
(*Histories*)

SwapAcceptance[res_List] :=
Module[{h = histOf[res]},
	If[h === $Failed, Message[CTLSurfaceData::badres]; Return[$Failed]];
	h[[4]][[All, "swapAccept"]]/h[[4]][[All, "swapTry"]]
];

End[];
EndPackage[];
