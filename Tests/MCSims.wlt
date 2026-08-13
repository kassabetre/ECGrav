(* ::Package:: *)

(* Regression tests for the MCSims-domain public functions.

   Two kinds of tests:
     - Deterministic: Hamiltonians and data-analysis helpers are pure functions
       of their inputs, so exact values are asserted (floats via a tolerance
       SameTest).
     - Smoke tests: the Monte Carlo drivers are stochastic. We assert that they
       run and return a well-formed result, and that the diagnostic-plot
       functions still emit their plots (guards the plot-restoration fix).
   Run with the Wolfram TestReport tool. *)

With[{local = If[StringQ[$InputFileName] && $InputFileName =!= "",
        FileNameJoin[{DirectoryName[$InputFileName], "TestPrelude.wl"}], ""]},
    Get[If[local =!= "" && FileExistsQ[local], local,
        "/Users/012759760/Desktop/Research/ECGravMathematicaPackage/Paclet/ECGrav/Tests/TestPrelude.wl"]]];

(* Inert Hamiltonian / delta wrappers, matching how the notebook drives the MC
   functions (which expect an inert head[params] argument).

   These MUST be plain multi-argument definitions, not curried ones. The drivers bind
   hamiltonian_[hparams___] and delH_[delHparams___], which requires h[-1.,0.] and dH[-1.,0.]
   to stay SYMBOLIC so that hamiltonian=h and hparams=(-1.,0.); the driver then calls
   delH[graph, delHparams, row, col] (MCSims.wl:1839). Written curried --
   dH[jj_,ll_] := delHIsing[#1,jj,ll,#2,#3]& -- the partial application evaluates to a
   Function, delH binds to Function itself, and that four-argument call hits a three-slot
   pure function: Function::argb every sweep, with the energy delta coming back as an
   unevaluated expression rather than a number.

   Note the patterns accept either shape, so nothing is rejected at the boundary; the damage
   only surfaces later. Same trap as passing a Hamiltonian symbol that does not exist. *)
h[am_List, jj_Real, ll_Real] := ECGrav`HIsing[am, jj, ll];
dH[am_List, jj_Real, ll_Real, i_Integer, j_Integer] := ECGrav`delHIsing[am, jj, ll, i, j];

(* ---------- Hamiltonians (deterministic) ---------- *)

VerificationTest[ECGrav`HIsing[K4, -1.0, 0.0], -12., SameTest -> floatEq, TestID -> "HIsing-K4"];
VerificationTest[ECGrav`HIsing[K4, 1.0, 0.5], 15., SameTest -> floatEq, TestID -> "HIsing-K4-field"];
VerificationTest[ECGrav`delHIsing[C6, -1.0, 0.0, 1, 2], 2., SameTest -> floatEq, TestID -> "delHIsing-C6"];
VerificationTest[
    ECGrav`HWeightedFaceCounts[K4, -1.0, 0.0, 0.0, 0.0, 0.0], -4.,
    SameTest -> floatEq, TestID -> "HWeightedFaceCounts-K4-vertices"
];
VerificationTest[ECGrav`HLaplacian[C6, 1.0], 0., SameTest -> floatEq, TestID -> "HLaplacian-cycle-is-0"];
VerificationTest[ECGrav`H1dCombManifold[C6, 1.0], 0., SameTest -> floatEq, TestID -> "H1dCombManifold-cycle-is-0"];

(* delHIsing must equal the difference of HIsing before/after toggling one edge *)
VerificationTest[
    With[{i = 1, j = 3, toggled = C6},
        ECGrav`delHIsing[C6, -1.0, 0.0, i, j]
        - (ECGrav`HIsing[ReplacePart[C6, {{i, j} -> 1 - C6[[i, j]], {j, i} -> 1 - C6[[j, i]]}], -1.0, 0.0]
           - ECGrav`HIsing[C6, -1.0, 0.0])],
    0.,
    SameTest -> floatEq, TestID -> "delHIsing-matches-HIsing-difference"
];

(* ---------- Aggregating-data helpers (deterministic) ---------- *)

VerificationTest[ECGrav`LogSumExp[{0., 0., 0.}], Log[3.], SameTest -> floatEq, TestID -> "LogSumExp-log3"];
VerificationTest[ECGrav`LogSumExp[{1., 2., 3.}], 3.4076059644443806, SameTest -> floatEq, TestID -> "LogSumExp-123"];
VerificationTest[ECGrav`CorrelationTime[0, {1., 2., 3., 4.}], 1.25, SameTest -> floatEq, TestID -> "CorrelationTime-lag0-variance"];
VerificationTest[ECGrav`SpecificHeat[{1., 2., 3., 4.}, 5, 1.5], 0.75, SameTest -> floatEq, TestID -> "SpecificHeat"];
VerificationTest[ECGrav`Susceptibility[{1., 2., 3., 4.}, 10, 1.5], 25.0, SameTest -> floatEq, TestID -> "Susceptibility"];
VerificationTest[
    ECGrav`InternalEnergy[0.5, <|0.5 -> -2.3|>, <|0.5 -> {1., 2., 3.}|>], 2.0,
    SameTest -> floatEq, TestID -> "InternalEnergy-single-beta"
];

(* ---------- Free-energy / reweighting helpers (deterministic oracles) ---------- *)

(* Empirical autocorrelation is 1 at lag 0 by construction (normalized); lag-1 hand-computed. *)
VerificationTest[ECGrav`EmpCorrelationTime[0, {1., 2., 3., 4.}], 1., SameTest -> floatEq, TestID -> "EmpCorrelationTime-lag0-is-1"];
VerificationTest[ECGrav`EmpCorrelationTime[1, {1., 3., 2., 4.}], -0.5, SameTest -> floatEq, TestID -> "EmpCorrelationTime-lag1"];

(* dCv/dT per site, hand-computed from (-beta^3/NN)(2(<E^2>-<E>^2)(1+beta<E>) - beta(<E^3>-<E><E^2>)). *)
VerificationTest[ECGrav`DSpecificHeat[{1., 2., 3., 4.}, 5, 1.5], -1.6875, SameTest -> floatEq, TestID -> "DSpecificHeat"];

(* Multi-histogram reweighting evaluated at the single measured beta collapses to closed
   forms (the free-energy factor cancels): -betaF equals the supplied value; the extrapolated
   observable equals its sample mean; Cv/T equals beta^2 (<E^2> - <E>^2); and the
   self-consistent free energy is gauge-fixed to 0 for a single temperature. *)
VerificationTest[
    ECGrav`NegativeBetaTimesFreeEnergy[0.5, <|0.5 -> -2.3|>, <|0.5 -> {1., 2., 3.}|>],
    -2.3, SameTest -> floatEq, TestID -> "NegativeBetaTimesFreeEnergy-single-beta"
];
VerificationTest[
    ECGrav`CvOverT[0.5, <|0.5 -> -2.3|>, <|0.5 -> {1., 2., 3.}|>],
    0.5^2 (Mean[{1., 4., 9.}] - 2.^2), SameTest -> floatEq, TestID -> "CvOverT-single-beta"
];
VerificationTest[
    ECGrav`ExtrapolatedExpectationValue[0.5, <|0.5 -> -2.3|>, <|0.5 -> {1., 2., 3.}|>, <|0.5 -> {10., 20., 30.}|>],
    20., SameTest -> floatEq, TestID -> "ExtrapolatedExpectationValue-single-beta"
];
VerificationTest[
    First[Values[Quiet[
        ECGrav`ComputeMinusBetaTimesFreeEnergy[<|0.5 -> {1., 2., 3.}|>], ParallelTable::nopar]]],
    0., SameTest -> floatEq, TestID -> "ComputeMinusBetaTimesFreeEnergy-single-beta"
];

(* ---------- Exact enumeration (deterministic) ---------- *)

(* K4 and C6 are the two lowest levels; energies -12 and -6.
   (The $KernelCount == 0 regression is guarded by a separate test at the end of
   this file, which must run last because it closes the parallel kernels.) *)
VerificationTest[
    Sort[Round[Keys[ECGrav`LowEnergyStates[{K4, C6}, ECGrav`HIsing[#, -1.0, 0.0] &, 2]]]],
    {-12, -6},
    TestID -> "LowEnergyStates-energies"
];

(* High-beta expectation of AvgDeg is dominated by K4 (avg degree 3). Passing a list of
   observable *functions* as the 3rd argument must dispatch cleanly: no expected-messages
   argument is given, so this also asserts that the stray Part::partd (emitted by the old
   numeric-obsValues pattern test while probing the overload) is now silenced. *)
VerificationTest[
    ECGrav`ExactExpectationValue[{K4, C6}, {{1.0, 1.0}}, {ECGrav`AvgDeg},
        ECGrav`HIsing[#, -1.0, 0.0] &, 1.0][[1, 3, 1]],
    2.9975273768433652,
    SameTest -> looseEq, TestID -> "ExactExpectationValue-AvgDeg-highbeta"
];

(* ---------- Monte Carlo drivers (stochastic smoke tests) ---------- *)

VerificationTest[
    Block[{Print = Null &}, Module[{r},
        SeedRandom[101];
        r = ECGrav`GraphMetropolis[C6, 0.02, ECGrav`HIsing[#, -1.0, 0.0] &, {Total[#, 2]/2.0 &}, 20];
        {Head[r], Length[r], Length[r[[2]]]}]],
    {List, 2, 20},
    TestID -> "GraphMetropolis-smoke"
];

VerificationTest[
    Block[{Print = Null &}, Module[{r},
        SeedRandom[102];
        r = ECGrav`GraphSweepReplica[C6, 0.2, h[-1.0, 0.0], 20, 0.2, 0];
        {Head[r] === List, KeyExistsQ[r[[1]], "minEnergy"]}]],
    {True, True},
    TestID -> "GraphSweepReplica-smoke"
];

(* The minimum a sweep reports must include the state the sweep STARTS from. The per-step
   comparison sits inside If[accept==1,...], so before this was fixed only a state entered by
   an ACCEPTED move could ever become the minimum: a sweep that started on a low
   configuration and accepted nothing reported minEToBeat itself -- an energy no visited
   state ever had -- with an empty minEstates to go with it.
   That is the mechanism behind the GraphCTLSchedule / GraphMultiHistogram complaint. A
   replica swap moves a configuration into a replica without any accepted move, and under a
   different external field its energy is a number no sweep has weighed, so the receiving
   replica never counted it.
   Nothing here is pinned to a seed. K4 is the exact ground state of h[-1.,0.] (energy -12.),
   and every move out of it costs +4, so at beta 10 the acceptance weight is Exp[-40] ~ 4*10^-18
   and the chain provably sits at -12. for all 120 steps. minEToBeat is tested both ABOVE the
   seed energy (0.) and BETWEEN it and zero (-5.), because the old code returned the threshold
   in both cases. Both overloads carried the same hole. *)
VerificationTest[
    Block[{Print = Null &}, Module[{withDelH, noDelH},
        withDelH = Table[
            ECGrav`GraphSweepReplica[K4, 10.0, h[-1.0, 0.0], dH[-1.0, 0.0], 20, thr, 0],
            {thr, {0.0, -5.0}}];
        noDelH = ECGrav`GraphSweepReplica[K4, 10.0, h[-1.0, 0.0], 20, 0.0, 0];
        {#[[1, "minEnergy"]] & /@ withDelH,
         #[[1, "minEstates"]] === {K4} & /@ withDelH,
         noDelH[[1, "minEnergy"]],
         noDelH[[1, "minEstates"]] === {K4}}]],
    {{-12., -12.}, {True, True}, -12., True},
    TestID -> "GraphSweepReplica-counts-its-seed"
];

(* The two numbers a sweep returns must be the hamiltonian's own values, not accumulations.
   The running energy is built incrementally (energy += delE) and drifts a few ulps over a
   sweep. On machine reals == and < share one RELATIVE tolerance, so that drift is normally
   invisible -- but a relative tolerance has nothing to scale at zero. h[0.1,0.1] is repulsive
   in both couplings, so its ground state is the empty graph at exactly 0., and the sweep used
   to report minima like -4.163336342344337*^-16: a minimum BELOW the least energy the model
   can produce, which then compares unequal to 0. and drops the states tying it.
   Order[a,b]==0 rather than ==, ===, or floatEq, because all three are tolerant on machine
   reals and would pass on exactly the values this is meant to reject. Order is bitwise.
   Several seeds because which one lands on the zero-energy state is not fixed; the assertions
   hold for every run either way. *)
VerificationTest[
    Block[{Print = Null &}, Module[{seedG, runs, hz},
        hz = ECGrav`HIsing[#, 0.1, 0.1] &;
        seedG = Normal[AdjacencyMatrix[CycleGraph[7]]];
        runs = Table[SeedRandom[s];
            ECGrav`GraphSweepReplica[seedG, 4.0, h[0.1, 0.1], dH[0.1, 0.1], 40, 1000.0, 0],
            {s, 6}];
        {(* no reported minimum lies below the model's own floor of zero *)
         AllTrue[runs, #[[1, "minEnergy"]] >= 0. &],
         (* each reported minimum is bit-for-bit the hamiltonian of the states reported with it *)
         AllTrue[runs, With[{r = #},
             AllTrue[r[[1, "minEstates"]], Order[hz[#], r[[1, "minEnergy"]]] == 0 &]] &],
         (* and the end state's energy is the hamiltonian's value, not the accumulated one *)
         AllTrue[runs, Order[hz[#[[2, "graph"]]], #[[2, "energy"]]] == 0 &],
         (* the run really does reach the zero-energy ground state, so the case above is live *)
         AnyTrue[runs, #[[1, "minEnergy"]] == 0. &]}]],
    {True, True, True, True},
    TestID -> "GraphSweepReplica-canonical-energies"
];

(* The same hole reached the exported drivers, which is where it was reported.
   GraphComputeCorrelationTime starts its running minimum AT minEToBeat with no states and
   relies entirely on GraphSweepReplica to improve it, so a chain frozen below the threshold
   returned the threshold and an empty minEstates. Same provably-frozen K4 construction as
   above, so this is deterministic rather than a seed pin.
   GraphEquilibriate is not covered here because it never had the hole: it initialises its
   minimum from hamiltonian[seedGraph] directly rather than from the threshold. *)
VerificationTest[
    Block[{Print = Null &}, Module[{ct, ctNoDelH},
        ct = ECGrav`GraphComputeCorrelationTime[K4, 10.0, h[-1.0, 0.0], dH[-1.0, 0.0], 21, 0.0, 1, 0];
        ctNoDelH = ECGrav`GraphComputeCorrelationTime[K4, 10.0, h[-1.0, 0.0], 21, 0.0, 1, 0];
        {ct[[1, "minEnergy"]], ct[[1, "minEstates"]] === {K4},
         ctNoDelH[[1, "minEnergy"]], ctNoDelH[[1, "minEstates"]] === {K4}}]],
    {-12., True, -12., True},
    TestID -> "GraphComputeCorrelationTime-counts-frozen-seed"
];

(* GraphEquilibriate must run AND emit its energy-vs-time plot. *)
VerificationTest[
    Module[{n}, SeedRandom[103];
        n = emittedGraphicsCount[ECGrav`GraphEquilibriate[C6, 0.2, h[-1.0, 0.0], 0]];
        n >= 1],
    True,
    TestID -> "GraphEquilibriate-emits-plot"
];

(* GraphComputeCorrelationTime must run AND emit its autocorrelation plot. *)
VerificationTest[
    Module[{n}, SeedRandom[104];
        n = emittedGraphicsCount[
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 21, 2.0, 1, 0]];
        n >= 1],
    True,
    TestID -> "GraphComputeCorrelationTime-emits-plot"
];

(* GraphMultiHistogram must emit its plots and return internally consistent output.
   The driver is not seed-reproducible -- two runs under the same SeedRandom give different
   numbers -- so there is no value to pin. Instead assert what must hold for ANY correct run:
   every replica's recorded energy equals the hamiltonian recomputed from its recorded graph,
   the recorded ground-state energy is the energy of the recorded ground states, and the
   reweighting produced finite free energies. Both halves come from one run because the call
   costs upwards of ten seconds. *)
VerificationTest[
    Module[{n, res, obs, HH}, SeedRandom[105];
        HH = ECGrav`HIsing[#, -1.0, 0.0] &;
        obs = {HH, Total[#, 2]/(10 (10 - 1)/2) &};
        (* The high-beta replicas of this run freeze onto the complete graph, which is the
           ground state of h[-1.0, 0.0], so both watched scalars sit still for the whole
           comparison window and GraphEquilibriate correctly reports that it had nothing to
           test. That is the diagnostic working, not a failure, so it is quieted by name --
           and only by name, so any other message still fails the test. *)
        {n, res} = emittedGraphicsAndResult[
            Quiet[ECGrav`GraphMultiHistogram[C6, 0.1, 1.2, h[-1.0, 0.0], obs, 40, 0],
                  ECGrav`GraphEquilibriate::nosignal]];
        {n >= 1,
         Length[res] === 3,
         groundStatesConsistentQ[res[[1]], HH, 6],
         (* groundStatesConsistentQ passes vacuously on an empty minEstates, which is exactly
            what the ground-state search used to return when no ACCEPTED move beat the
            threshold, so demand states and demand that the reported minimum is at or below
            every energy the run wrote into its own chart (column 3). That last inequality is
            the reported symptom stated directly: a returned minimum above an energy the run
            recorded is a minimum that is not the minimum. *)
         res[[1]]["minEstates"] =!= {},
         res[[1]]["minEnergy"] <= Min[Join @@ (#["data"][[All, 3]] & /@ Values[res[[2]]])] + 10.^-9,
         mcReplicasConsistentQ[res[[3]], HH, 6],
         AllTrue[Values[res[[3]]], IntegerQ[#["eqlT"]] && #["eqlT"] > 0 && #["corrT"] >= 2 &],
         Sort[Keys[res[[3]]]] === Sort[Keys[res[[2]]]],
         AllTrue[Values[res[[2]]], MatchQ[#["minusBetaF"], _Real] &]}],
    {True, True, True, True, True, True, True, True, True},
    TestID -> "GraphMultiHistogram-invariants"
];

(* GraphParallelTempering must run cleanly on a tiny job (regression for the beta-history
   Part::take). NN=1 with two nearby betas drives a replica's swapAccept up to NN, which
   used to underflow the fixed-length history buffer: the trim slice -(swapAccept+1);;-1
   asked for positions -2..-1 of a length-1 list. SeedRandom[1] deterministically reaches
   swapAccept==NN on both replicas, exercising that boundary. With no expected-messages
   argument this also guards the Mod[_,0] progress-print glitch (printCase==0 when NN<5).
   The inner drivers emit diagnostic plots via Print, so Print is blocked. *)
VerificationTest[
    Block[{Print = Null &}, Module[{r},
        SeedRandom[1];
        r = ECGrav`GraphParallelTempering[C6, {0.1, 0.2}, h[-1.0, 0.0], dH[-1.0, 0.0],
            {Total[#, 2]/2. &}, 0, 1, 0, -100.0];
        {Head[r], Length[r]}]],
    {List, 4},
    TestID -> "GraphParallelTempering-smoke-no-Part-take"
];

(* GraphParallelTempering must return internally consistent output. Same reasoning as the
   GraphMultiHistogram invariants above: nothing can be pinned to a value, so assert the
   consistency any correct run has.
   minEToBeat is 0.0 rather than the -100.0 used by the smoke test above, because nothing
   ever beats -100.0 for this system: minEstates comes back EMPTY and the ground-state check
   passes vacuously. At 0.0 the threshold is reachable, so the check has something to test.
   Also asserts the swap bookkeeping, which is pure counting and cannot depend on the RNG. *)
VerificationTest[
    Block[{Print = Null &}, Module[{r, HH}, SeedRandom[1];
        HH = ECGrav`HIsing[#, -1.0, 0.0] &;
        r = ECGrav`GraphParallelTempering[C6, {0.1, 0.2}, h[-1.0, 0.0], dH[-1.0, 0.0],
            {Total[#, 2]/2. &}, 0, 1, 0, 0.0];
        {Length[r] === 4,
         r[[1]]["minEstates"] =!= {},
         groundStatesConsistentQ[r[[1]], HH, 6],
         (* As in the multihistogram invariants: the reported minimum must be at or below
            every energy the run recorded in its chart (column 3), keyed here by beta. *)
         r[[1]]["minEnergy"] <= Min[Join @@ (#[[All, 3]] & /@ Values[r[[2]]])] + 10.^-9,
         mcReplicasConsistentQ[r[[3]], HH, 6],
         Sort[#["beta"] & /@ Values[r[[3]]]] === {0.1, 0.2},
         AllTrue[Values[r[[4]]], 0 <= #["swapAccept"] <= #["swapTry"] &],
         AllTrue[Values[r[[4]]], Length[#["history"]] > 0 &]}]],
    {True, True, True, True, True, True, True, True},
    TestID -> "GraphParallelTempering-invariants"
];

(* ---------- Ground-state search ---------- *)

(* GradDescent is deterministic steepest descent: it returns a valid (6x6, symmetric, 0/1)
   adjacency matrix and never increases the energy -- the defining invariant of descent. *)
VerificationTest[
    Block[{Print = Null &}, Module[{r},
        r = ECGrav`GradDescent[C6, dH[-1.0, 0.0], 10];
        {Dimensions[r], r === Transpose[r], SubsetQ[{0, 1}, Union@Flatten[r]],
         ECGrav`HIsing[r, -1.0, 0.0] <= ECGrav`HIsing[C6, -1.0, 0.0]}]],
    {{6, 6}, True, True, True},
    TestID -> "GradDescent-descent-invariant"
];

(* SGradDescent (softmax stochastic descent) returns <|minE, minEstates, LastState|>;
   smoke-test both the with-delH and no-delH overloads. *)
VerificationTest[
    Block[{Print = Null &}, Module[{r}, SeedRandom[201];
        r = ECGrav`SGradDescent[C6, h[-1.0, 0.0], dH[-1.0, 0.0], 0.5, 5];
        {Head[r], KeyExistsQ[r, "minE"], NumberQ[r["minE"]]}]],
    {Association, True, True},
    TestID -> "SGradDescent-with-delH-smoke"
];
VerificationTest[
    Block[{Print = Null &}, Module[{r}, SeedRandom[202];
        r = ECGrav`SGradDescent[C6, h[-1.0, 0.0], 0.5, 5];
        {Head[r], KeyExistsQ[r, "minE"]}]],
    {Association, True},
    TestID -> "SGradDescent-no-delH-smoke"
];

(* SimulatedAnnealing (inert h[..]/dH[..] form) anneals betai -> betaf and returns an
   association including the minimum energy visited. *)
VerificationTest[
    Block[{Print = Null &}, Module[{r}, SeedRandom[203];
        r = ECGrav`SimulatedAnnealing[C6, h[-1.0, 0.0], dH[-1.0, 0.0], 0.1, 2.0, 5, 3];
        {Head[r], KeyExistsQ[r, "minE"], NumberQ[r["minE"]]}]],
    {Association, True, True},
    TestID -> "SimulatedAnnealing-smoke"
];

(* ---------- Bug #3 regression: LowEnergyStates with no parallel kernels ----------
   MUST BE LAST: CloseKernels[] forces $KernelCount == 0 so the divide-by-zero
   guard (Max[$KernelCount, 1]) is exercised. Without the fix this failed with
   Power::infy and returned $Failed. The ParallelTable::nopar warning (parallel
   call falling back to sequential with no kernels) is benign and quieted. *)
VerificationTest[
    Module[{r}, CloseKernels[];
        r = Quiet[
            Sort[Round[Keys[ECGrav`LowEnergyStates[{K4, C6}, ECGrav`HIsing[#, -1.0, 0.0] &, 2]]]],
            ParallelTable::nopar];
        r],
    {-12, -6},
    TestID -> "LowEnergyStates-no-kernels-divide-guard"
];


(* The two user settings have to reach subkernels, because a subkernel runs its own copy of
   the package and would otherwise read its own defaults. SyncParallelSettings used to do
   that with DistributeDefinitions, which cannot: it silently skips every symbol whose context
   is a loaded package, and ECGrav` is one as soon as the subkernels have the package. It
   returned {} and shipped nothing while the helper recorded success, so every parallel driver
   went on using the defaults and a user raising the budget saw no effect.

   Asserted through OwnValues on the subkernel and never through ParallelEvaluate[sym]: that
   inlines the master's value at send time and reports the master's setting whether or not
   anything arrived, which is exactly what kept the broken version looking correct. The
   DistributeDefinitions leg is the negative control -- it is what the assertion above would
   look like if the old mechanism were still in place.

   Subkernels are loaded from source first, because ParallelNeeds["ECGrav`"] on its own
   resolves to the *installed* paclet and would test whatever was last built rather than this
   tree; the ParallelNeeds after it is a no-op for loading but registers the context, and
   registration is what makes DistributeDefinitions skip the symbols. Loading the package on
   subkernels without registering it does not reproduce the failure -- which is why this is the
   configuration to test, since it is the one the drivers are used in. The pool is closed and
   relaunched afterwards so later files get the prelude's clean state. *)
VerificationTest[
    Module[{pushed, shipped, src = FileNameJoin[{$ECGravRoot, "Kernel", "ECGrav.wl"}]},
        CloseKernels[]; Quiet[LaunchKernels[2]];
        With[{path = src}, ParallelEvaluate[Get[path]]];
        ParallelNeeds["ECGrav`"];
        pushed = Block[{ECGrav`$ECGravMaxEquilibriationSweeps = 777,
                        ECGrav`$ECGravMaxCorrelationSweeps = 333},
            ECGrav`Private`SyncParallelSettings[];
            ParallelEvaluate[{Last /@ OwnValues[ECGrav`$ECGravMaxEquilibriationSweeps],
                              Last /@ OwnValues[ECGrav`$ECGravMaxCorrelationSweeps]}]];
        shipped = DistributeDefinitions[ECGrav`$ECGravMaxCorrelationSweeps];
        CloseKernels[]; Quiet[LaunchKernels[]];
        {DeleteDuplicates[pushed], shipped}],
    {{{{777}, {333}}}, {}},
    TestID -> "SyncParallelSettings-reaches-subkernels"
];
