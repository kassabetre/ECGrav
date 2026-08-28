#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*TemperingMixing -- how well the replicas mixed*)

(*
	Load after Needs["ECGrav`"].  Analyses the histories a GraphParallelTempering (or
	GraphCTLSchedule) run returns, and answers the question HomogeneousHamiltonian.md section 7
	step 4 says is the one that matters: not what the swap acceptance was, but whether
	configurations actually travelled from the hot end to the cold end and back.

	WHAT "history" ACTUALLY HOLDS -- this is the thing to get right before reading anything here.

	On an accepted swap between slots k and l the driver exchanges the two STATES and leaves each
	slot's external field where it is (MCSims.wl:7999-8006), then pushes each slot's newest history
	entry onto the other (8011-8020).  So a slot is a fixed field value, and what moves is the
	configuration.  Each configuration is tagged by the field of the slot it started in, and

		histories[[k, "history"]]  =  the ordered list of configuration TAGS that occupied slot k

	which is the transpose of what it looks like.  It is NOT a walk through field space: on this
	run only 32% of consecutive entries in slot 1's buffer are even adjacent to slot 1 in the swap
	graph, where a trajectory would be 100% by construction.

	RECOVERING THE TRAJECTORIES.  Each accepted swap appends one entry to each of the two slots
	involved, and the tags always form a permutation of the slots.  That is enough to replay the
	whole run: an event at slots k and l is admissible only when each slot's next unread entry is
	the other's current occupant, and any two admissible events are necessarily slot-disjoint --
	if they shared a slot they would both claim its next entry, which names one tag.  Disjoint
	events commute, so the replay is order-independent and the per-configuration trajectories are
	unique.  Verified on a 2402-swap run: identical under first / last / random tie-breaking, and
	every buffer read to its end.

	The swap graph is not needed -- the tag constraint picks the events on its own, and the graph
	comes back out as a by-product (41 edges recovered, exactly the schedule's).
*)

BeginPackage["TemperingMixing`", {"ECGrav`"}];

MixingData::usage =
"MixingData[res] replays the run behind a GraphParallelTempering or GraphCTLSchedule return and \
bundles the trajectories, the swap statistics and the rung ordering.  Rungs are numbered by \
increasing beta, so rung 1 is the hottest and rung n the coldest.";

RoundTrips::usage =
"RoundTrips[mix] counts full traversals per configuration.  A round trip is coldest rung -> \
hottest rung -> coldest rung; the reverse direction is reported too.";

ReplicaFlow::usage =
"ReplicaFlow[mix] gives the Katzgraber replica-flow fraction f at each rung: of the visits to \
that rung, the share made by configurations currently heading toward the hot end.  It should \
fall smoothly from 1 at the cold end to 0 at the hot end; a flat stretch is a bottleneck.";

OccupancyKind::usage =
"OccupancyKind[mix] is \"sweeps\" when the run's chart carried the per-sweep permutation \
(column 7) and \"arrivals\" when the trajectories had to be replayed from the histories.";

OccupancyMatrix::usage =
"OccupancyMatrix[mix] gives m[[rung, configuration]]: sweeps spent there when the chart carries \
the per-sweep permutation, arrivals when the trajectories were replayed.  OccupancyKind says which.";

MixingReport::usage =
"MixingReport[mix] summarises the run: round trips, per-rung acceptance, flow, occupancy \
uniformity and the correlation-time flags.";

MixingPlots::usage =
"MixingPlots[mix] returns an association of diagnostic plots: \"RoundTrips\", \"Acceptance\", \
\"Flow\", \"Occupancy\", \"Trajectory\".";

TrajectoryPlot::usage =
"TrajectoryPlot[mix, j] draws configuration j's rung against swap number.";

Begin["`Private`"];

histOf[res_List] := Which[
	Length[res] === 2 && ListQ[res[[2]]] && Length[res[[2]]] === 4, res[[2]],
	Length[res] === 4 && AssociationQ[res[[2]]], res,
	True, $Failed];

MixingData::badres = "Expected the 4-element return of GraphParallelTempering or GraphCTLSchedule.";
MixingData::incomplete =
"Replay stalled with `1` of `2` accepted swaps consumed.  The histories are inconsistent with a \
pure exchange of states -- check that the run was not resumed or edited.";

(*Consecutive duplicates removed: the WALK, one entry per move. Round trips and arrival counts
	mean the same thing in both modes because both are computed from this.*)
walkOf[perSweep_List] := First /@ Split[perSweep];

(*Chart mode. Column 7 names the configuration occupying each slot at each sweep, so the whole
	permutation history is there and nothing has to be reconstructed. Swap[] exchanges a matching
	of non-intersecting edges, so a configuration moves at most once per sweep and this is the
	complete trajectory at sweep resolution -- and, unlike the replay, it carries DWELL.*)
perSweepFromChart[chart_Association, slots_List, n_Integer, col_Integer] :=
Module[{idx = AssociationThread[slots -> Range[n]], tagAt, nSweeps, traj},
	tagAt = Map[idx, Values[chart[[All, All, col]]], {2}];
	nSweeps = Length[First[tagAt]];
	traj = ConstantArray[0, {n, nSweeps}];
	Do[traj[[tagAt[[k, t]], t]] = k, {k, n}, {t, nSweeps}];
	traj
];

(*Two shapes reach here. GraphParallelTempering keys its chart straight to the measurement rows;
	GraphCTLSchedule and GraphMultiHistogram wrap them as <|"minusBetaF"->…, "data"->rows|>. Reading
	Length[First[Values[chart]]] without normalising gives the number of KEYS in that wrapper -- 3
	or 4 -- as the sweep count, which is wrong and silently plausible.*)
chartRows[chart_Association] :=
	If[AssociationQ[First[Values[chart]]], chart[[All, "data"]], chart];

(*Which column, if any, carries the configuration tag.  Testing the row LENGTH is not enough and
	is actively wrong: the beta-schedule chart splices its observables into the row, so an
	eight-observable run has eleven columns and a "length >= 7" test reads observable four as a
	tag, producing Missing[KeyAbsent, …] where a slot index was wanted.

	Test the invariant instead -- a tag column holds a key of the chart in every row, and at every
	sweep the column read across slots is a PERMUTATION of the keys, since a swap exchanges
	configurations and never creates or destroys one.  Candidates are tried from the last column
	backwards, so this keeps working if the layout gains columns.*)
tagColumnOf[rows_Association, slots_List] :=
Module[{vals = Values[rows], len, keys = Alternatives @@ slots, ok},
	len = Length[First[First[vals]]];
	If[!ListQ[First[First[vals]]] || len < 4, Return[None]];
	ok[c_] := AllTrue[vals, MatchQ[#[[1, c]], keys] &] &&
		AllTrue[Transpose[vals[[All, All, c]]], Sort[#] === Sort[slots] &];
	SelectFirst[Range[len, 4, -1], ok, None]
];

(*A single-parameter run keys its chart on a bare Real, a multi-field one on a vector. The rung
	ordering is by beta either way, which is the first component or the number itself.*)
betaOf[k_] := If[ListQ[k], First[k], k];
nt0Of[k_] := If[ListQ[k] && Length[k] > 1, k[[2]]/k[[1]], 0.];

repField[reps_Association, key_String] :=
	If[AllTrue[Values[reps], KeyExistsQ[#, key] &],
		Values[reps[[All, key]]], Missing["NotAvailable"]];

MixingData[res_List] :=
Module[{h4, chart, hist, slots, n, idx, tags, ptr, atSlot, traj, ready, e, k, l, tk, tl,
		pairs, fired, order, reps, accepted, perSweep = None, source, tagCol},

	h4 = histOf[res];
	If[h4 === $Failed, Message[MixingData::badres]; Return[$Failed]];
	chart = chartRows[h4[[2]]]; reps = h4[[3]]; hist = h4[[4]];
	slots = Keys[chart];
	n = Length[slots];
	idx = AssociationThread[slots -> Range[n]];
	tags = Table[idx /@ hist[[k0, "history"]], {k0, n}];

	(*Replay.  Every unordered pair is a candidate; the tag constraint is what selects, and the
		pairs that fire come back as the swap graph.*)
	pairs = UndirectedEdge @@@ Subsets[Range[n], {2}];
	ptr = ConstantArray[1, n];
	atSlot = Range[n];
	traj = Table[{k0}, {k0, n}];
	fired = {};
	While[True,
		ready = Select[pairs, Function[ed, With[{a = ed[[1]], b = ed[[2]]},
			ptr[[a]] < Length[tags[[a]]] && ptr[[b]] < Length[tags[[b]]] &&
			tags[[a, ptr[[a]] + 1]] === atSlot[[b]] &&
			tags[[b, ptr[[b]] + 1]] === atSlot[[a]]]]];
		If[ready === {}, Break[]];
		e = First[ready]; k = e[[1]]; l = e[[2]]; tk = atSlot[[k]]; tl = atSlot[[l]];
		atSlot[[k]] = tl; atSlot[[l]] = tk;
		ptr[[k]]++; ptr[[l]]++;
		AppendTo[traj[[tl]], k]; AppendTo[traj[[tk]], l];
		AppendTo[fired, Sort[{k, l}]]];

	accepted = Total[Values[hist[[All, "swapAccept"]]]]/2;
	If[ptr =!= (Length /@ tags),
		Message[MixingData::incomplete, Length[fired], accepted]];

	(*Rungs by increasing beta: rung 1 hottest, rung n coldest.  In homogeneous c.O form beta is
		c_0, which is what the ladder is a ladder IN.*)
	order = Ordering[betaOf /@ slots];

	(*Prefer the chart. The replay is exact, but it recovers only the ORDER of moves; column 7
		additionally says how long each configuration stayed, which is what occupancy and the
		replica flow are properly weighted by.*)
	(*The chart's first row is written AFTER that sweep's Swap[], so the initial permutation is
		not in it -- configuration j starts at slot j, which is what the histories seeding encodes
		(MCSims.wl:7954). Prepend it to the WALK only; residence and the flow stay on the recorded
		sweeps, where every sweep is a whole sweep. Split then collapses the duplicate for any
		configuration that did not move in sweep 1. With this the chart walk and the replay walk
		agree element for element.*)
	tagCol = tagColumnOf[chart, slots];
	If[tagCol =!= None,
		perSweep = perSweepFromChart[chart, slots, n, tagCol];
		traj = Table[walkOf[Prepend[perSweep[[j]], j]], {j, n}];
		source = "chart",
		source = "replay"];

	<|"slots" -> slots, "n" -> n, "order" -> order, "source" -> source, "tagColumn" -> tagCol,
	  "perSweep" -> perSweep,
	  "betas" -> (betaOf /@ slots),
	  "nt0" -> (nt0Of /@ slots),
	  "trajectories" -> traj,
	  "swapAccept" -> Values[hist[[All, "swapAccept"]]],
	  "swapTry" -> Values[hist[[All, "swapTry"]]],
	  "swapGraph" -> UndirectedEdge @@@ DeleteDuplicates[fired],
	  "acceptedSwaps" -> Length[fired],
	  "sweeps" -> Length[First[Values[chart]]],   (*chart is normalised to rows above*)
	  (*Not every driver records every one of these -- the beta-schedule path omits corrTValues
		and eqlT entirely -- so they are looked up rather than parted.*)
	  "corrT" -> repField[reps, "corrT"],
	  "corrTValues" -> repField[reps, "corrTValues"],
	  "corrTMeasured" -> repField[reps, "corrTMeasured"],
	  "eqlT" -> repField[reps, "eqlT"],
	  "complete" -> (ptr === (Length /@ tags))|>
];

(*Hottest and coldest rung as SLOT indices.*)
ends[mix_] := {First[mix["order"]], Last[mix["order"]]};

countTrips[walk_List, from_, to_] :=
Module[{lab = None, c = 0},
	Do[Which[s === to, If[lab === "from", c++]; lab = "to",
			 s === from, lab = "from"], {s, walk}];
	c];

RoundTrips[mix_Association] :=
Module[{hot, cold},
	{hot, cold} = ends[mix];
	(*Every mapped value is parenthesised: "key" -> f[#] & /@ list parses as
		("key" -> f[#]) & /@ list, a LIST OF RULES sharing one key, which an association then
		collapses to the last element -- silently, as a scalar where a list was meant.*)
	<|"perConfiguration" -> (countTrips[#, cold, hot] & /@ mix["trajectories"]),
	  "perConfigurationReverse" -> (countTrips[#, hot, cold] & /@ mix["trajectories"]),
	  "total" -> Total[countTrips[#, cold, hot] & /@ mix["trajectories"]],
	  "withAtLeastOne" -> Count[countTrips[#, cold, hot] & /@ mix["trajectories"], _?Positive],
	  "perConfigurationPerSweep" ->
		N[Total[countTrips[#, cold, hot] & /@ mix["trajectories"]]/(mix["n"]*mix["sweeps"])],
	  "meanRoundTripSweeps" ->
		With[{tot = Total[countTrips[#, cold, hot] & /@ mix["trajectories"]]},
			If[tot === 0, Indeterminate, N[mix["n"]*mix["sweeps"]/tot]]],
	  "visitsHot" -> (Count[#, hot] & /@ mix["trajectories"]),
	  "visitsCold" -> (Count[#, cold] & /@ mix["trajectories"])|>
];

(*Katzgraber replica flow.  Label a configuration "up" from the moment it touches the coldest
	rung -- it can only go toward hot next -- and "down" from the moment it touches the hottest.
	f is the share of up-labelled visits at each rung.  A schedule with no bottleneck gives a
	smooth monotone descent; a plateau says configurations are piling up and turning back.*)
(*Weighted by SWEEPS when the chart carries the per-sweep permutation, which is the definition --
	f is a histogram of walkers over rungs, and a walker that lingers should count for the time it
	lingered. Falling back to arrivals, as a replayed run must, over-weights the high-degree rungs
	that are entered and left most often.*)
ReplicaFlow[mix_Association] :=
Module[{hot, cold, up, down, lab, series},
	{hot, cold} = ends[mix];
	series = If[mix["perSweep"] === None, mix["trajectories"], mix["perSweep"]];
	up = ConstantArray[0, mix["n"]];
	down = ConstantArray[0, mix["n"]];
	Do[
		lab = None;
		Do[
			Which[s === cold, lab = "up", s === hot, lab = "down"];
			Which[lab === "up", up[[s]]++, lab === "down", down[[s]]++],
			{s, walk}],
		{walk, series}];
	MapThread[If[#1 + #2 === 0, Indeterminate, N[#1/(#1 + #2)]] &, {up, down}]
];

(*Sweeps spent per rung when the chart carries them, arrivals otherwise. The two differ by more
	than a scale: a rung with six swap-graph neighbours is entered and left more often than one
	with four, so arrivals carry a structural bias that residence does not.*)
OccupancyMatrix[mix_Association] :=
	Transpose@Table[
		Module[{c = Counts[If[mix["perSweep"] === None,
					mix["trajectories"][[j]], mix["perSweep"][[j]]]]},
			Table[Lookup[c, k, 0], {k, mix["n"]}]],
		{j, mix["n"]}];

OccupancyKind[mix_Association] := If[mix["perSweep"] === None, "arrivals", "sweeps"];

MixingReport[mix_Association] :=
Module[{rt = RoundTrips[mix], f = ReplicaFlow[mix], m = OccupancyMatrix[mix], o = mix["order"],
		rate, occFrac, comps},
	rate = mix["swapAccept"]/mix["swapTry"];
	(*Uniformity of arrivals across rungs, per configuration: 1 means perfectly even.*)
	occFrac = Table[With[{col = m[[All, j]]},
		If[Total[col] === 0, 0., N[Min[col]/Max[col]]]], {j, mix["n"]}];
	(*A disconnected swap graph makes round trips impossible BY CONSTRUCTION -- configurations
		cannot reach rungs outside their own component -- so a zero count means the ladder was
		severed, not that the sampling was slow. Report it next to the count rather than leaving it
		to be inferred.*)
	comps = ConnectedComponents[Graph[Range[mix["n"]], mix["swapGraph"]]];
	<|"source" -> mix["source"],
	  "occupancyIn" -> OccupancyKind[mix],
	  "swapGraphConnected" -> (Length[comps] === 1),
	  "swapGraphComponents" -> Sort[Length /@ comps, Greater],
	  "roundTripsTotal" -> rt["total"],
	  "configurationsWithARoundTrip" -> StringForm["`1` of `2`", rt["withAtLeastOne"], mix["n"]],
	  "roundTripsPerConfiguration" -> N[rt["total"]/mix["n"]],
	  "meanRoundTripSweeps" -> rt["meanRoundTripSweeps"],
	  "sweeps" -> mix["sweeps"],
	  "acceptedSwaps" -> mix["acceptedSwaps"],
	  "acceptanceRange" -> MinMax[rate],
	  "acceptanceSpread" -> N[Max[rate]/Min[rate]],
	  "rungsVisitedPerConfiguration" -> (Length[DeleteDuplicates[#]] & /@ mix["trajectories"]),
	  "occupancyEvenness" -> occFrac,
	  "flowByRung" -> f[[o]],
	  (*f runs 0 at the hottest rung to 1 at the coldest by construction, so "monotone" means
		NON-DECREASING in mix["order"], which is hot to cold. Testing the reverse asks whether it
		decreases and answers False for every healthy ladder.*)
	  "flowIsMonotone" -> OrderedQ[f[[o]]],
	  "replayComplete" -> mix["complete"],
	  "corrT" -> mix["corrT"],
	  "corrTMeasured" -> mix["corrTMeasured"]|>
];

(* ::Section:: *)
(*Plots*)

MixingPlots[mix_Association] :=
Module[{o = mix["order"], rt = RoundTrips[mix], f = ReplicaFlow[mix], m = OccupancyMatrix[mix],
		bs, rate},
	bs = mix["betas"][[o]];
	rate = (mix["swapAccept"]/mix["swapTry"])[[o]];
	<|
	"RoundTrips" -> BarChart[rt["perConfiguration"],
		ChartLabels -> Range[mix["n"]], PlotLabel -> "round trips per configuration",
		AxesLabel -> {"configuration", None}, PlotRange -> All],
	"Acceptance" -> ListLinePlot[Transpose[{bs, rate}], PlotMarkers -> Automatic,
		PlotLabel -> "swap acceptance per rung (CTL aims to flatten this)",
		AxesLabel -> {"\[Beta]", "accepted / tried"}, PlotRange -> {0, All}],
	"Flow" -> ListLinePlot[Transpose[{bs, f[[o]]}], PlotMarkers -> Automatic,
		PlotLabel -> "replica flow f (1 at the cold end, 0 at the hot end)",
		AxesLabel -> {"\[Beta]", "f"}, PlotRange -> {0, 1}],
	"Occupancy" -> MatrixPlot[m[[o]], FrameLabel -> {"configuration", "rung (hot \[Rule] cold)"},
		PlotLabel -> "arrivals per rung per configuration", ColorFunction -> "TemperatureMap"],
	"Trajectory" -> TrajectoryPlot[mix, 1]
	|>
];

TrajectoryPlot[mix_Association, j_Integer] :=
Module[{o = mix["order"], rank, chartMode = mix["perSweep"] =!= None, series},
	rank = AssociationThread[o -> Range[mix["n"]]];
	series = If[chartMode, mix["perSweep"][[j]], mix["trajectories"][[j]]];
	ListLinePlot[rank /@ series,
		PlotLabel -> StringForm["configuration `1`: rung against `2`", j,
			If[chartMode, "sweep", "accepted swap"]],
		AxesLabel -> {If[chartMode, "sweep", "accepted swap"], "rung (1 = hottest)"},
		PlotRange -> {1, mix["n"]}, InterpolationOrder -> 0]
];

End[];
EndPackage[];
