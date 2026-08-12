(* ::Package:: *)

(* Shared setup for the ECGrav .wlt test suite.

   Each .wlt file Get[]s this prelude first. It loads the package straight
   from source (Kernel/ECGrav.wl) rather than from an installed paclet, so the
   tests always exercise the current working tree, and defines the fixtures and
   float-comparison helpers used across the suite. *)

(* Paclet root. Under TestReport, $InputFileName does not point at this file,
   so prefer the known absolute root and only use a $InputFileName-derived path
   if it actually resolves to this prelude. Tests/ is a sibling of Kernel/. *)
$ECGravRoot = "/Users/012759760/Desktop/Research/ECGravMathematicaPackage/Paclet/ECGrav";
If[StringQ[$InputFileName] && $InputFileName =!= "" &&
        FileExistsQ[FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]], "Kernel", "ECGrav.wl"}]],
    $ECGravRoot = ParentDirectory[DirectoryName[$InputFileName]]];

(* CI / relocation override: an explicit ECGRAV_ROOT environment variable, when it points
   at a real checkout, wins over both the hardcoded default and the $InputFileName guess.
   This lets the suite run from a repo checked out to any path (e.g. GitHub Actions), where
   the hardcoded absolute path does not exist and $InputFileName does not resolve under
   TestReport. Tests/ci-run.wls sets this before invoking TestReport. *)
With[{envRoot = Environment["ECGRAV_ROOT"]},
    If[StringQ[envRoot] && envRoot =!= "" &&
            FileExistsQ[FileNameJoin[{envRoot, "Kernel", "ECGrav.wl"}]],
        $ECGravRoot = envRoot]];

Get[FileNameJoin[{$ECGravRoot, "Kernel", "ECGrav.wl"}]];

(* Much of MCSims runs on ParallelTable, and LowEnergyStates divides by
   $KernelCount, so it hard-fails (0^(-1)) if no parallel kernels exist yet.
   The notebook workflow always calls LaunchKernels[] before using the package;
   mirror that here so the suite exercises the package as it is actually used. *)
If[$KernelCount == 0, Quiet[LaunchKernels[]]];

(* ----- Float comparison helpers (for use as VerificationTest SameTest) ----- *)

(* scalar / same-shape numeric arrays, absolute tolerance.
   The difference is wrapped in a list before flattening so that scalar
   comparisons (where N[#1]-N[#2] is atomic) work as well as array ones. *)
approxEq[tol_] := (Max[Abs[Flatten[{N[#1] - N[#2]}]]] < tol &);
floatEq = approxEq[10.^-6];
looseEq = approxEq[10.^-3];

(* ----- Fixtures: pure simplicial complexes (facet lists) ----- *)

fst1 = {{1, 2, 3}, {4, 5, 6}, {2, 4, 7}};
fst2 = {{1, 2, 3}, {3, 5, 6}, {3, 5, 7}, {4, 9, 10}};

(* Triangulated surfaces with known topology (independent oracles) *)
tetrahedron = {{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}};
octahedron  = {{1, 2, 3}, {1, 2, 4}, {1, 3, 5}, {1, 4, 5},
               {6, 2, 3}, {6, 2, 4}, {6, 3, 5}, {6, 4, 5}};
mobiusStrip1 = {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {4, 5, 1}, {5, 1, 2}};
mobiusStrip2 = {{1, 2, 3}, {2, 3, 5}, {3, 4, 5}, {4, 5, 6}, {2, 4, 6}, {1, 2, 6}};
kleinbottle = {{1, 2, 5}, {2, 5, 6}, {2, 3, 6}, {3, 6, 7}, {1, 3, 7}, {1, 5, 7},
               {4, 5, 6}, {4, 6, 8}, {6, 7, 8}, {7, 8, 9}, {5, 7, 9}, {4, 5, 9},
               {1, 4, 8}, {1, 3, 8}, {3, 8, 9}, {2, 3, 9}, {2, 4, 9}, {1, 2, 4}};
projectivePlane = {{1, 2, 3}, {1, 2, 5}, {1, 3, 4}, {1, 4, 6}, {1, 5, 6},
                   {2, 3, 6}, {2, 4, 5}, {2, 4, 6}, {3, 4, 5}, {3, 5, 6}};
torus = {{13, 15, 16}, {7, 13, 15}, {14, 15, 16}, {10, 14, 16}, {3, 10, 14},
         {10, 12, 16}, {5, 10, 12}, {11, 12, 16}, {11, 13, 16}, {2, 11, 13},
         {7, 9, 15}, {5, 9, 12}, {8, 14, 15}, {8, 9, 15}, {5, 8, 9}, {2, 5, 8},
         {6, 11, 12}, {6, 9, 12}, {6, 7, 9}, {3, 6, 7}, {4, 8, 14}, {4, 6, 11},
         {3, 4, 14}, {3, 4, 6}, {2, 4, 11}, {2, 4, 8}, {1, 5, 10}, {1, 7, 13},
         {1, 3, 10}, {1, 3, 7}, {1, 2, 13}, {1, 2, 5}};

(* Seed complexes for the MCMC driver tests. mcmcSeed24 uses the whole p*M label pool, so
   the free-vertex chain starts at the maximum vertex count; mcmcSeed34 has p = 3, which the
   fixed-vertex tests need because for p = 2 the facets are the edges and the edge-count
   observable could never fluctuate. *)
mcmcSeed24 = Partition[Range[8], 2];                          (* p = 2, M = 4, on [8] *)
mcmcSeed34 = {{1, 2, 3}, {1, 4, 5}, {2, 4, 6}, {3, 5, 7}};    (* p = 3, M = 4, on [7] *)

(* The seven surfaces in the order used by the notebook's orientability test *)
surfaces = {mobiusStrip1, mobiusStrip2, kleinbottle, projectivePlane,
            torus, tetrahedron, octahedron};

(* Clique-graph views of a couple of the complexes *)
octaG = ECGrav`GraphFromCliques[octahedron];
tetraG = ECGrav`GraphFromCliques[tetrahedron];
torusG = ECGrav`GraphFromCliques[torus];  (* clique complex is a clean torus triangulation *)

(* ----- Fixtures: graph adjacency matrices for the MC / Hamiltonian tests ----- *)

K4 = Normal[AdjacencyMatrix[CompleteGraph[4]]];
C6 = Normal[AdjacencyMatrix[CycleGraph[6]]];

(* Count graphics objects a (stochastic) function Print[]s to the front end.
   Used to assert that the diagnostic-plot functions still emit their plots. *)
SetAttributes[emittedGraphicsCount, HoldAll];
emittedGraphicsCount[expr_] := Module[{captured},
    captured = Reap[Block[{Print = Sow[{##}, "grph"] &}, expr], "grph"][[2]];
    Count[Flatten[captured, Infinity], _Graphics | _Legended | _GraphicsBox, Infinity]
];

(* As above, but hands back the driver's return value too, so a test can assert the plots
   AND the result from a single (expensive) run. *)
SetAttributes[emittedGraphicsAndResult, HoldAll];
emittedGraphicsAndResult[expr_] := Module[{captured, res},
    captured = Reap[Block[{Print = Sow[{##}, "grph"] &}, res = expr], "grph"][[2]];
    {Count[Flatten[captured, Infinity], _Graphics | _Legended | _GraphicsBox, Infinity], res}
];

(* Collect the arguments of every `name` message an expression issues, as a list of
   argument lists (empty if it never fires). Two reasons this is not done with Check or
   $MessageList: it recovers the filled-in slots, not just the fact that something fired,
   and the handler runs before Quiet and before the General::stop repeat limit, so the
   caller can wrap the expression in Quiet -- keeping the test log clean -- and still see
   every message, including the fourth and later copies of the same one. *)
SetAttributes[messageArguments, HoldRest];
messageArguments[name_String, expr_] := Module[{collected = {}},
    Internal`HandlerBlock[{"Message", Function[m,
        Replace[m, Hold[Message[MessageName[s_, tag_], args___], _] :>
            If[SymbolName[s] <> "::" <> tag === name, AppendTo[collected, {args}]]]]},
        expr];
    collected
];

(* The correlation-time integral that RandomPureSimplicialComplexMCMCCorrelationTime runs
   over each measured observable: sum the autocorrelation, normalized by its lag-0 value,
   over lags 0, 1, 2, ... and stop at the first non-positive term, which is itself included
   in the sum; if the autocorrelation never turns over, stop one lag short of lastLag. Both
   are deliberate conventions carried over from the pre-lazy tabulated form, so an oracle
   that reproduces the driver's answer has to reproduce them too.

   Ceiling of the sum, floored at 2, is what the driver reports. Returns that alongside
   which of the two exits was taken, so a test can show it exercised both.

   Which convention applies is only visible in the answer on a SHORT run. The terminating
   non-positive term is by definition the first one at or below zero, so it is small, and
   over a few hundred sweeps one lag at the end of the range is worth about 10^-3 -- far
   below Ceiling. Over 20 sweeps the same lag is a tenth of the integral, which is why the
   short-run correlation-time test exists. *)
mcmcCorrT[series_, lastLag_] :=
    Module[{norm, corrSum = 0.0, corr, t = 0, turnedOver = False},
        norm = ECGrav`CorrelationTime[0, series];
        If[norm == 0.0, norm = 1.0];
        While[!turnedOver && t <= lastLag,
            corr = ECGrav`CorrelationTime[t, series]/norm;
            Which[
                corr <= 0, corrSum += corr; turnedOver = True,
                t + 1 <= lastLag, corrSum += corr];
            t++];
        {Max[Ceiling[corrSum], 2], turnedOver}
    ];

(* Every measurement row a RandomPureSimplicialComplexMCMC run Sows carries both the complex
   it measured and the numbers it measured on it:
   {sweep index, energy, vertex or edge count, one value per operator, complex}. Recomputing
   those numbers from the recorded complex must reproduce the recorded ones, which ties the
   reported numbers to the reported states -- a row that kept a stale complex, or columns that
   slipped out of step with the operator list, fails here where a shape check would not.

   The recomputation goes back through the package's own sweep at NN = 0: no steps run, so
   what comes back is the weight and count it derived from the complex it was handed. Extra
   arguments are passed straight through, so callers write ...[meas, ops, 1] for the
   free-vertex chain and ...[meas, ops, True, 1] for the fixed-vertex one. *)
mcmcMeasurementsConsistentQ[measurements_, operators_, sweepArgs___] :=
    AllTrue[measurements, Function[row,
        Module[{complex = Last[row], probe},
            probe = ECGrav`RandomPureSimplicialComplexMCMCSweep[complex, sweepArgs, 0];
            ECGrav`PureComplexQ[complex] &&
            floatEq[row[[2]], probe["energy"]] &&
            row[[3]] === If[KeyExistsQ[probe, "vertexCount"], probe["vertexCount"],
                                                              probe["edgeCount"]] &&
            Take[row, {4, -2}] === Through[operators[complex]]]]];

(* ----- Invariants for the stochastic MC drivers -----
   GraphMultiHistogram and GraphParallelTempering are not seed-reproducible: two runs under
   the same SeedRandom give different numbers, so nothing can be pinned to an expected value.
   What CAN be asserted is internal consistency of whatever a run produced -- properties that
   hold for any correct run and fail on malformed output. *)

(* A valid n-vertex simple-graph adjacency matrix. *)
adjacencyMatrixQ[m_, n_Integer] := MatrixQ[m, IntegerQ] && Dimensions[m] === {n, n} &&
    m === Transpose[m] && SubsetQ[{0, 1}, Union @ Flatten[m]] &&
    Diagonal[m] === ConstantArray[0, n];

(* Every replica must carry a valid graph whose energy, recomputed from scratch, matches the
   energy the driver recorded. This ties the reported numbers to the reported states, so a
   replica that was never actually evolved -- or that came back malformed -- fails here. *)
mcReplicasConsistentQ[replicas_, ham_, n_Integer] :=
    AllTrue[Values[replicas],
        adjacencyMatrixQ[#["state"]["graph"], n] &&
        floatEq[#["state"]["energy"], ham[#["state"]["graph"]]] &];

(* The recorded ground-state energy must be the energy of the recorded ground states. Note
   minEstates is legitimately empty when nothing beat the minEToBeat threshold, in which case
   minEnergy is that threshold and there is nothing to cross-check. *)
groundStatesConsistentQ[gs_, ham_, n_Integer] :=
    AllTrue[gs["minEstates"], adjacencyMatrixQ[#, n] &] &&
    (gs["minEstates"] === {} || floatEq[gs["minEnergy"], Min[ham /@ gs["minEstates"]]]);
