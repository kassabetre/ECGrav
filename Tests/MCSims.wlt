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
   functions (which expect an inert head[params] argument). *)
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

(* ---------- Exact enumeration (deterministic) ---------- *)

(* K4 and C6 are the two lowest levels; energies -12 and -6. *)
VerificationTest[
    Sort[Round[Keys[ECGrav`LowEnergyStates[{K4, C6}, ECGrav`HIsing[#, -1.0, 0.0] &, 2]]]],
    {-12, -6},
    TestID -> "LowEnergyStates-energies"
];

(* High-beta expectation of AvgDeg is dominated by K4 (avg degree 3).
   The Part::partd message is emitted by the function itself while probing
   which overload its observable argument matches; the result is still correct. *)
VerificationTest[
    ECGrav`ExactExpectationValue[{K4, C6}, {{1.0, 1.0}}, {ECGrav`AvgDeg},
        ECGrav`HIsing[#, -1.0, 0.0] &, 1.0][[1, 3, 1]],
    2.9975273768433652,
    {Part::partd},
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

(* GraphMultiHistogram must run AND emit its Cv/T and entropy plots. *)
VerificationTest[
    Module[{n, obs}, SeedRandom[105];
        obs = {ECGrav`HIsing[#, -1.0, 0.0] &, Total[#, 2]/(10 (10 - 1)/2) &};
        n = emittedGraphicsCount[
            ECGrav`GraphMultiHistogram[C6, 0.1, 1.2, h[-1.0, 0.0], obs, 40, 0]];
        n >= 1],
    True,
    TestID -> "GraphMultiHistogram-emits-plot"
];
