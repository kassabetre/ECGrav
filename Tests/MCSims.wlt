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

(* The thermodynamic metric that places CTL replicas is the Fisher information, which for an
   exponential family in these natural parameters is the COVARIANCE of the conjugate observables.
   It used to be assembled from the second moments alone -- E[O_i O_j], under a comment claiming
   Cov -- with every entry passed through Max[0.0,...]. Both halves were wrong, and both are
   invisible at one external field: the missing subtraction is then a smooth positive inflation,
   and the clip never fires because a second moment E[O^2] cannot be negative.

   The fixture is chosen so a single case exercises both. The two observables have nonzero means,
   so the subtraction matters; they are anti-correlated, so the covariance is negative; and their
   cross second moment is *also* negative, so the old clip zeroed it outright. Under the old
   expression this returns {{5.5, 0.}, {0., 3.}} against a true {{5.25, -3.25}, {-3.25, 2.75}} --
   the diagonal inflated by E[O]^2 and the off-diagonal destroyed, which is exactly the cross term
   that makes one multi-field schedule better than one schedule per axis.

   The oracle is the population covariance computed directly from the samples, so it shares no
   code with the function under test. Rows 3 and 4 are the vacuity guards: the result must differ
   from the second-moment matrix, and the off-diagonal must come back negative. *)
VerificationTest[
    Module[{a = {-2., -1., 1., 4.}, b = {3., 1., -1., -1.}, m1, m2, pop, oracle, f},
        f = ECGrav`Private`CTLMetricFromMoments;
        pop[u_, v_] := Mean[(u - Mean[u]) (v - Mean[v])];
        m1 = {Mean[a], Mean[b]};
        m2 = {{Mean[a^2], Mean[a b]}, {Mean[a b], Mean[b^2]}};
        oracle = {{pop[a, a], pop[a, b]}, {pop[a, b], pop[b, b]}};
        {Max[Abs[Flatten[f[1.0, m1, m2] - oracle]]] < 10.^-12,   (* is the covariance *)
         Max[Abs[Flatten[f[2.0, m1, m2] - 4 oracle]]] < 10.^-12, (* scales as bt^2 *)
         f[1.0, m1, m2] =!= m2,                                  (* not the second moment *)
         Negative[f[1.0, m1, m2][[1, 2]]],                       (* off-diagonal keeps its sign *)
         f[1.0, {1., 1.}, {{0.999, -0.5}, {-0.5, 2.}}],          (* diagonal floored, off-diag not *)
         f[1.0, {Mean[a]}, {{Mean[a^2]}}]}],                     (* one field, unchanged *)
    {True, True, True, True, {{0., -1.5}, {-1.5, 1.}}, {{5.25}}},
    TestID -> "CTLMetricFromMoments-is-the-covariance"
];

(* The CTL interpolates its metric over the bounding BOX of the bootstrap field table, and a box
   is not the region the bootstrap run covered -- the corners are reached by extrapolation, and a
   metric invented there can place replicas where nothing was sampled. The effective sample size
   is what the density is masked with to prevent that, so it has to be right in two regimes: equal
   to the sample count where every sample contributes, and collapsing towards 1 where one sample
   carries the weight.

   The oracle is the Kish ratio written directly, without logs, sharing no code with the
   implementation. It is used here as a correctness check only: the direct form does not fail
   under extrapolation, since WL escapes to arbitrary precision rather than overflowing, but that
   escape is why the shipped version works in log space -- 34x slower one field-width out on a
   60-sample fixture, and unfinished after ten minutes on a 720-sample one. *)
VerificationTest[
    Module[{mbF, meas, mbF1, meas1, oracle, ess},
        ess = ECGrav`Private`MBAREffectiveSampleSize;
        oracle[bt_, h_, f_, m_] := Module[{ks = Keys[f], nv = Length /@ m, u},
            u = Flatten@Table[Table[
                1/Sum[nv[[Key[j]]] Exp[-f[[Key[j]]]] Exp[bt (h - j) . m[[Key[i], s]]], {j, ks}],
                {s, 1, nv[[Key[i]]]}], {i, ks}];
            Total[u]^2/Total[u^2]];
        mbF = <|{0.} -> 0., {1.} -> 0.|>;
        meas = <|{0.} -> {{1.}, {2.}, {3.}}, {1.} -> {{1.}, {2.}, {3.}}|>;
        (* one field with the target on it: every denominator is identical, so every sample
           carries the same weight and the effective count is the actual count *)
        mbF1 = <|{0.} -> 0.3|>;
        meas1 = <|{0.} -> {{1.}, {5.}, {9.}, {-2.}}|>;
        {AllTrue[{0., 0.5, 60., 800.},
            Abs[ess[1.0, {#}, mbF, meas] - oracle[1.0, {#}, mbF, meas]] < 10.^-9 &],
         Abs[ess[1.0, {0.}, mbF1, meas1] - 4] < 10.^-9,     (* equal weights -> N *)
         1 <= ess[1.0, {0.5}, mbF, meas] <= 6,              (* bounded by 1 and N *)
         ess[1.0, {0.}, mbF, meas] > 5,                     (* well supported *)
         ess[1.0, {800.}, mbF, meas] < 2.5}],               (* extrapolated: collapses *)
    {True, True, True, True, True},
    TestID -> "MBAREffectiveSampleSize-matches-Kish-and-collapses-off-support"
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

(* The deep-reject regime, where the acceptance test used to stop being representable.

   From the complete graph with J = -1 every proposal REMOVES an edge and so is uphill by 8,
   making beta*delE = 3200 at beta 400 -- far past the ~709 where Exp[-beta*delE] ceases to be
   a normalized double. The old form underflowed to 0. and rejected, which was the right answer
   arrived at by luck; the log form compares Log[u] < -3200 and rejects on the arithmetic.
   Either way the chain must sit still and return exact numbers, which is what this pins.

   Both labelings are asserted because they take different paths: the unlabeled one multiplies
   by a ratio of automorphism group orders, and that product underflows in its own right even
   where Exp alone would not have.

   This does not prove the message is gone -- that was checked separately against the previous
   code (55 underflow messages versus 0, with identical energies, magnetizations and edge
   counts over runs at beta 7.5 and 20). What it guards is that the rewrite did not disturb
   the deep-reject path, which is where a sign slip or a wrong comparison threshold would
   show up as a chain that wanders instead of freezing. *)
VerificationTest[
    Module[{K6 = Normal[AdjacencyMatrix[CompleteGraph[6]]], runs},
        runs = Table[SeedRandom[77];
            Quiet@Block[{Print = Null &},
                ECGrav`GraphSweepReplica[K6, 400.0, h[-1.0, 0.0], dH[-1.0, 0.0], 40, 0.0, u]],
          {u, {0, 1}}];
        {(* nothing moved: every uphill proposal was rejected *)
         AllTrue[runs, #[[2, "graph"]] === K6 &],
         (* and the reported numbers are exact reals, not underflow debris *)
         AllTrue[runs, NumberQ[#[[2, "energy"]]] && #[[2, "energy"]] == -60. &],
         AllTrue[runs, #[[1, "minEnergy"]] == -60. &],
         (* the regime really is past the representable range, so the case is live *)
         400.0*ECGrav`delHIsing[Normal[AdjacencyMatrix[CompleteGraph[6]]], -1.0, 0.0, 1, 2] > 709}],
    {True, True, True, True},
    TestID -> "GraphSweepReplica-deep-reject-at-extreme-beta"
];

(* beta enters the sweep in exactly one place -- logRatio = Log[selectionProb] - delE*beta -- and
   nothing else in the driver reads it. That is what makes a HOMOGENEOUS reparameterisation
   possible: rewrite the hamiltonian so every coupling is explicit, H = Sum_i c_i O_i, run at
   bt = 1 with c scaled by beta, and the chain is not merely similar but IDENTICAL, because the
   code only ever forms the product beta*delH.

   This matters because the multi-field tempering path (GraphCTLSchedule, GraphMultiHistogram,
   GraphParallelTempering) tempers the field vector c at FIXED bt, so beta is not one of its
   coordinates. Promoting the beta-conjugate part of the hamiltonian to a conjugate observable
   with its own coupling makes it one, and that trick is only sound while this identity holds. If
   a second use of beta is ever added to the acceptance, this test is what notices.

   HIsing is linear and homogeneous in (J,L) with no constant term, so scaling beta into the
   couplings is exactly the reparameterisation: beta*HIsing[am,J,L] == HIsing[am,beta J,beta L],
   checked here to 4.4e-16. The last column is the vacuity guard -- a frozen chain would return
   the seed from both runs and compare equal for no reason. *)
VerificationTest[
    Module[{res},
        res = Table[
            Module[{orig, hom},
                SeedRandom[9091];
                orig = ECGrav`GraphSweepReplica[C6, bt, h[-1.0, 0.3], dH[-1.0, 0.3], 25, 0.0, 0];
                SeedRandom[9091];
                hom = ECGrav`GraphSweepReplica[C6, 1.0, h[bt*-1.0, bt*0.3],
                        dH[bt*-1.0, bt*0.3], 25, 0.0, 0];
                {orig[[2, Key["graph"]]] === hom[[2, Key["graph"]]],
                 Abs[bt*orig[[2, Key["energy"]]] - hom[[2, Key["energy"]]]] < 10.^-9,
                 orig[[2, Key["graph"]]] =!= C6}],
            {bt, {0.4, 0.7, 2.0, 5.0}}];
        {AllTrue[res, #[[1]] &],                         (* identical chain *)
         AllTrue[res, #[[2]] &],                         (* recorded energy scales by beta *)
         AllTrue[res, #[[3]] &],                         (* and the chain actually moved *)
         Abs[0.7*ECGrav`HIsing[C6, -1.0, 0.3]
             - ECGrav`HIsing[C6, 0.7*-1.0, 0.7*0.3]] < 10.^-12}],
    {True, True, True, True},
    TestID -> "GraphSweepReplica-homogeneous-reparameterisation-is-identical"
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

(* ---------- GraphComputeCorrelationTime: the integration rule ----------

   Differential test of the reported correlation times against the independent implementation
   of the documented integration rule in the prelude (mcmcCorrT), which is the same rule both
   sides of the package run. This is what pins the lazy lag walk: it computes only the lags up
   to the turnover, and must still produce exactly what tabulating every lag produced.

   The driver measures its own observables internally, so the test recovers them: one supplied
   operator records the graph it is handed, which is the post-sweep graph the driver measures
   in the same iteration, and the energy and magnetization are deterministic functions of it
   (mag being Total[Flatten[g]]/(v(v-1)), written exactly as the driver writes it).

   The second operator is a sweep counter -- a monotone ramp whose autocorrelation never turns
   over, which forces the lag walk to run to the end of its range instead of stopping after a
   handful of terms. Without it every value sits on the floor of 2 and the comparison is
   nearly vacuous: turnedOver below records that BOTH branches of the integral were taken --
   terminating on a non-positive term, and running out of lags one short of lastLag.

   Vacuity guard: the same oracle on shuffled series, whose autocorrelation is destroyed, must
   give a different answer. *)
VerificationTest[
    Module[{recorded = {}, sweepNo = 0, r, numsweeps, rows, oracle, shuffled, nV = 6},
        SeedRandom[601];
        r = Quiet@Block[{Print = Null &},
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0], 31, 0.0,
                {Function[g, AppendTo[recorded, g]; 1.0*Total[Flatten[g]]],
                 Function[g, 1.0*(++sweepNo)]}, 0]];
        numsweeps = Length[recorded];
        rows = {ECGrav`HIsing[#, -1.0, 0.0] & /@ recorded,
                Total[Flatten[#]]*1.0/(nV (nV - 1)) & /@ recorded,
                1.0*Total[Flatten[#]] & /@ recorded,
                1.0*Range[numsweeps]};
        oracle = mcmcCorrT[#, numsweeps - 10] & /@ rows;
        SeedRandom[9];
        shuffled = mcmcCorrT[RandomSample[#], numsweeps - 10][[1]] & /@ rows;
        {numsweeps, Last[recorded] === r[[2, "state", "graph"]],
         r[[2, "corrTValues"]],
         r[[2, "corrTValues"]] === oracle[[All, 1]],
         r[[2, "corrT"]] === Max[oracle[[All, 1]]],
         oracle[[All, 2]],                        (* turnover reached on each series? *)
         shuffled =!= oracle[[All, 1]]}],
    {155, True, {3, 2, 2, 53}, True, True, {True, True, True, False}, True},
    TestID -> "GraphComputeCorrelationTime-integration-rule"
];

(* The same rule on the overload that takes no delH. It is a separate body from the one above
   -- and until this was written it was the only one of the four with no stuck-observable
   handling at all -- so it gets its own oracle rather than being assumed equivalent. *)
VerificationTest[
    Module[{recorded = {}, sweepNo = 0, r, numsweeps, rows, oracle, nV = 6},
        SeedRandom[602];
        r = Quiet@Block[{Print = Null &},
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 31, 0.0,
                {Function[g, AppendTo[recorded, g]; 1.0*Total[Flatten[g]]],
                 Function[g, 1.0*(++sweepNo)]}, 0]];
        numsweeps = Length[recorded];
        rows = {ECGrav`HIsing[#, -1.0, 0.0] & /@ recorded,
                Total[Flatten[#]]*1.0/(nV (nV - 1)) & /@ recorded,
                1.0*Total[Flatten[#]] & /@ recorded,
                1.0*Range[numsweeps]};
        oracle = mcmcCorrT[#, numsweeps - 10] & /@ rows;
        {numsweeps, r[[2, "corrTValues"]],
         r[[2, "corrTValues"]] === oracle[[All, 1]],
         r[[2, "corrT"]] === Max[oracle[[All, 1]]],
         oracle[[All, 2]]}],
    {155, {2, 2, 2, 53}, True, True, {True, True, True, False}},
    TestID -> "GraphComputeCorrelationTime-integration-rule-no-delH"
];

(* $ECGravCorrelationRunMultiplier is the other half of the pair that sets the measurement run,
   and the half that binds whenever the multiple is below the cap -- which at the default cap of
   1000 is every eqlT under 200. Until it became a setting the multiple was a literal 5, so a
   user who raised the cap to lengthen a short run got no change whatever, and had no knob that
   worked. The recording operator counts the sweeps directly, as on the complex side.

   The fractional case is here because a multiplier, unlike a sweep cap, invites one: 2.5*31 is
   77.5, and the run length has to come out a whole number or numsweeps-10 is not an Integer and
   LazyCorrelationTime does not match its own argument pattern -- it would return unevaluated and
   be indexed as if it were an association. Rounded up, so no setting can ask for fewer sweeps
   than it names. *)
VerificationTest[
    Module[{sweepsFor},
        sweepsFor[mult_, eqlT_, cap_] := Module[{rec = {}},
            Block[{ECGrav`$ECGravCorrelationRunMultiplier = mult,
                   ECGrav`$ECGravMaxCorrelationSweeps = cap},
                Quiet@Block[{Print = Null &}, SeedRandom[603];
                    ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0],
                        eqlT, 0.0,
                        {Function[g, AppendTo[rec, g]; 1.0*Total[Flatten[g]]]}, 0]]];
            Length[rec]];
        {sweepsFor[5, 31, 1000],    (* the former literal, unmoved *)
         sweepsFor[10, 31, 1000],   (* raising it lengthens the run ... *)
         sweepsFor[2, 31, 1000],    (* ... and lowering it shortens it *)
         sweepsFor[2.5, 31, 1000],  (* fractional: Ceiling[77.5], not 77 or unevaluated *)
         sweepsFor[10, 31, 100]}],  (* the cap still wins where it is the smaller *)
    {155, 310, 62, 78, 100},
    TestID -> "GraphComputeCorrelationTime-run-length-multiplier"
];

(* A frozen observable must be excluded, named, and reported with ITS OWN stuck value.
   Both halves were broken. The overload with delH indexed the value two rows above the one it
   named, so for the operator at row 4 it printed the magnetization's value; the overload
   without delH did no filtering at all, leaving a Null in the autocorrelation table which
   indexed its way into corrT and returned it as an unevaluated Sum rather than a number.
   Hence the IntegerQ assertions: a shape check alone passed on the broken output.

   The constant operator sits at row 4 (energy, mag, fluctuating operator, frozen operator),
   which is exactly where the misindexed value differed from the right one. *)
VerificationTest[
    Module[{args3, args4, r3, r4, ops},
        ops = {Function[g, 1.0*Total[Flatten[g]]], Function[g, 7.5]};
        SeedRandom[603];
        args3 = messageArguments["GraphComputeCorrelationTime::stuck",
            r3 = Quiet@Block[{Print = Null &},
                ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0],
                    21, 0.0, ops, 0]]];
        SeedRandom[603];
        args4 = messageArguments["GraphComputeCorrelationTime::stuck",
            r4 = Quiet@Block[{Print = Null &},
                ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 21, 0.0, ops, 0]]];
        {args3, args4,
         VectorQ[r3[[2, "corrTValues"]], IntegerQ], IntegerQ[r3[[2, "corrT"]]],
         VectorQ[r4[[2, "corrTValues"]], IntegerQ], IntegerQ[r4[[2, "corrT"]]],
         Length[r3[[2, "corrTValues"]]], Length[r4[[2, "corrTValues"]]]}],
    {{{Function[g, 7.5], 7.5}}, {{Function[g, 7.5], 7.5}},
     True, True, True, True, 4, 4},
    TestID -> "GraphComputeCorrelationTime-frozen-observable"
];

(* When nothing fluctuates, both operator-list overloads report ::alldefault and leave every
   correlation time at the floor. K4 at beta 10 is the provably-frozen construction used by
   GraphComputeCorrelationTime-counts-frozen-seed above, so this is deterministic. *)
VerificationTest[
    Module[{ad3, ad4, r3, r4, ops = {Function[g, 1.0*Total[Flatten[g]]]}},
        ad3 = messageArguments["GraphComputeCorrelationTime::alldefault",
            r3 = Quiet@Block[{Print = Null &},
                ECGrav`GraphComputeCorrelationTime[K4, 10.0, h[-1.0, 0.0], dH[-1.0, 0.0],
                    21, 0.0, ops, 0]]];
        ad4 = messageArguments["GraphComputeCorrelationTime::alldefault",
            r4 = Quiet@Block[{Print = Null &},
                ECGrav`GraphComputeCorrelationTime[K4, 10.0, h[-1.0, 0.0], 21, 0.0, ops, 0]]];
        {Length[ad3], Length[ad4], ad3 === ad4,
         ad3[[1, 1]], ad3[[1, 2]],
         r3[[2, "corrTValues"]], r4[[2, "corrTValues"]],
         r3[[2, "corrT"]], r4[[2, "corrT"]]}],
    {1, 1, True, {{-12.}, {1.}, {12.}}, {2, 2, 2}, {2, 2, 2}, {2, 2, 2}, 2, 2},
    TestID -> "GraphComputeCorrelationTime-all-frozen"
];

(* A run too short to have a lag range at all must still return a NUMBER.
   eqlT 1 gives numsweeps = 5, so lastLag = numsweeps - 10 = -5 and there is no lag to
   integrate over. The tabulated form took FirstPosition's default -- a bare integer rather
   than a position list -- indexed it, and handed back corrT as an unevaluated expression;
   every one of the four overloads did this, and a caller then multiplied its sweep count by
   it. IntegerQ is the whole point of the test: the broken values were still Max[...] shaped
   and compared equal to nothing. *)
VerificationTest[
    Module[{ops = {Function[g, 1.0*Total[Flatten[g]]]}, rs},
        SeedRandom[604];
        rs = Quiet@Block[{Print = Null &},
            {ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0], 1, 0.0, 1, 0],
             ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 1, 0.0, 1, 0],
             ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0], 1, 0.0, ops, 0],
             ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 1, 0.0, ops, 0]}];
        {#[[2, "corrT"]] & /@ rs,
         AllTrue[rs, IntegerQ[#[[2, "corrT"]]] &],
         VectorQ[rs[[3, 2, "corrTValues"]], IntegerQ],
         VectorQ[rs[[4, 2, "corrTValues"]], IntegerQ]}],
    {{2, 2, 2, 2}, True, True, True},
    TestID -> "GraphComputeCorrelationTime-short-run-returns-number"
];

(* corrTMeasured separates the two ways corrT can come back as 2: computed from an
   autocorrelation integral that really did land on the floor, versus never computed at all.
   The number alone cannot tell them apart -- the same problem the "converged" key was added to
   the equilibriators to fix, and the reason a stuck run is not simply given a bigger corrT:
   corrT is the sweep count the drivers pass to GraphSweepReplica, so inflating it would
   multiply the cost of every production run that froze, and at high beta freezing usually
   means the ground state was found rather than that anything went wrong.

   Four cases across all four overloads. A normal run measures it. K4 at beta 10 is the
   provably-frozen construction used above: nothing fluctuates, ::alldefault fires, and the
   flag is False. A single constant operator alongside a live one must NOT flip it, since
   corrT is the max over the observables that do fluctuate. And eqlT 1 gives numsweeps 5, so
   there is no lag range at all and no term enters the integral -- the case no message reports,
   where the flag is the only signal. *)
VerificationTest[
    Module[{normal, allFrozen, oneStuck, tooShort, ops = {Function[g, 1.0*Total[Flatten[g]]]}},
        SeedRandom[606];
        normal = Quiet@Block[{Print = Null &},
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0], 21, 0.0, 1, 0]];
        allFrozen = Quiet@Block[{Print = Null &},
            ECGrav`GraphComputeCorrelationTime[K4, 10.0, h[-1.0, 0.0], dH[-1.0, 0.0], 21, 0.0, ops, 0]];
        SeedRandom[607];
        oneStuck = Quiet@Block[{Print = Null &},
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 21, 0.0,
                Append[ops, Function[g, 7.5]], 0]];
        SeedRandom[608];
        tooShort = Quiet@Block[{Print = Null &},
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 1, 0.0, 1, 0]];
        {normal[[2, "corrTMeasured"]],
         allFrozen[[2, "corrTMeasured"]],
         oneStuck[[2, "corrTMeasured"]],
         tooShort[[2, "corrTMeasured"]],
         (* all four report corrT as an ordinary integer, so the flag carries information the
            number does not -- and the two False cases both report exactly 2 *)
         AllTrue[{normal, allFrozen, oneStuck, tooShort}, IntegerQ[#[[2, "corrT"]]] &],
         {allFrozen[[2, "corrT"]], tooShort[[2, "corrT"]]}}],
    {True, False, True, False, True, {2, 2}},
    TestID -> "GraphComputeCorrelationTime-corrTMeasured"
];

(* The same run through ALL FOUR overloads, which is what the test above only appears to do:
   each of its four cases exercises a different overload, so every overload gets exactly one
   verdict, and the operator-list-with-delH form drew the frozen case -- the one where the
   right answer is False anyway. That overload set corrTMeasured nowhere at all: it initialised
   the key to False and never assigned it, so it agreed with the test for the wrong reason
   while reporting "never measured" for every external-field tempering run in the package, next
   to corrTValues holding integrals in the hundreds.

   The shape of the bug is the 1.8.1 one, and so is the shape of the fix: assert the behaviour
   at every CALL SITE, not once on a representative. Two rows, because a flag that is hardwired
   to either constant passes a one-row test -- eqlT 21 gives a normal run and all four must say
   True; eqlT 1 gives numsweeps 5, no lag range, and all four must say False. *)
VerificationTest[
    Module[{ops = {Function[g, 1.0*Total[Flatten[g]]]}, four},
        four[eq_] := (SeedRandom[611]; Quiet@Block[{Print = Null &},
            {ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0], eq, 0.0, 1, 0],
             ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], eq, 0.0, 1, 0],
             ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0], eq, 0.0, ops, 0],
             ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], eq, 0.0, ops, 0]}]);
        {#[[2, "corrTMeasured"]] & /@ four[21],
         #[[2, "corrTMeasured"]] & /@ four[1],
         #[[2, "corrT"]] & /@ four[1]}],
    {{True, True, True, True}, {False, False, False, False}, {2, 2, 2, 2}},
    TestID -> "GraphComputeCorrelationTime-corrTMeasured-every-overload"
];

(* Both operator-list overloads emit the autocorrelation plot. The one without delH did not
   draw one at all before, which is why it is asserted separately from the plot test above. *)
VerificationTest[
    Module[{ops = {Function[g, 1.0*Total[Flatten[g]]]}}, SeedRandom[605];
        {emittedGraphicsCount[
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], dH[-1.0, 0.0], 21, 2.0, ops, 0]] >= 1,
         emittedGraphicsCount[
            ECGrav`GraphComputeCorrelationTime[C6, 0.2, h[-1.0, 0.0], 21, 2.0, ops, 0]] >= 1}],
    {True, True},
    TestID -> "GraphComputeCorrelationTime-operator-overloads-emit-plot"
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

(* The external-field chart carries a SEVENTH column: the tag of the configuration occupying
   each slot at that sweep, read off histories[[i, "history", -1]] immediately after Swap[].

   The invariant that makes it worth anything is that at every sweep the column, read across
   slots, is a PERMUTATION of the schedule's field vectors -- swaps exchange configurations, so
   none is ever lost or duplicated. That is asserted here sweep by sweep rather than by spot
   check, because a stale or mis-timed read would still produce plausible-looking field vectors.

   Swap[] exchanges a matching of non-intersecting edges, so a configuration moves at most once
   per sweep and the column is a complete trajectory at sweep resolution; that is what lets
   TemperingMixing read dwell times instead of replaying the histories. *)
VerificationTest[
    Block[{Print = Null &},
    Module[{r, chart, labels, edges, o0, o1, hHom, col7},
        o0 = 1.0*Total[#, 2]/2 &;
        o1 = 1.0*Tr[MatrixPower[#, 3]]/6 &;
        hHom[am_List, c0_Real, c1_Real] := c0*o0[am] + c1*o1[am];
        labels = <|1 -> {0.3, 0.05}, 2 -> {0.55, 0.05}, 3 -> {0.8, 0.05}|>;
        edges = {UndirectedEdge[1, 2], UndirectedEdge[2, 3]};
        SeedRandom[612];
        r = ECGrav`GraphParallelTempering[C6, 1.0, hHom[], {edges, labels},
                {o0, o1}, {Total[#, 2]/2. &}, 12, 0, -100.0];
        chart = r[[2]];
        col7 = Transpose[Values[chart[[All, All, 7]]]];
        {Length[First[First[Values[chart]]]],
         Sort[Keys[chart]] === Sort[Values[labels]],
         AllTrue[col7, Sort[#] === Sort[Values[labels]] &],
         Length[col7]}]],
    {7, True, True, 12},
    TestID -> "GraphParallelTempering-chart-column7-is-a-permutation"
];

(* The beta-schedule chart is nested too, and for the same reason: {sweep, beta, energy,
   {obs...}, tag}. It used to splice the observables into the row with Flatten, which put every
   column after the third at a position depending on how many observables the caller passed --
   so no fixed-position column could follow them, and this row could not carry a tag at all.

   Energy stays at column 3 in both layouts, which is the only column the package itself reads
   out of a beta chart (two SpecificHeat call sites and
   ComputeMinusBetaTimesFreeEnergy[chart[[All, All, 3]]]), so the change is internal-consumer
   clean; that is asserted here rather than assumed. *)
VerificationTest[
    Block[{Print = Null &},
    Module[{r, chart, row},
        SeedRandom[613];
        r = ECGrav`GraphParallelTempering[C6, {0.1, 0.3, 0.6, 1.0}, h[-1.0, 0.0], dH[-1.0, 0.0],
                {Total[#, 2]/2. &, 1.0*ECGrav`EulerChi[#] &}, 0, 20, 0, 0.0];
        chart = r[[2]];
        row = chart[[1, 1]];
        {Length[row],
         ListQ[row[[4]]], Length[row[[4]]],          (* the observables, nested not spliced *)
         MemberQ[Keys[chart], row[[5]]],             (* the tag *)
         AllTrue[Transpose[Values[chart[[All, All, 5]]]], Sort[#] === Sort[Keys[chart]] &],
         VectorQ[chart[[1, All, 3]], NumericQ],      (* energy still at column 3 *)
         AssociationQ[ECGrav`ComputeMinusBetaTimesFreeEnergy[chart[[All, All, 3]]]]}]],
    {5, True, 2, True, True, True, True},
    TestID -> "GraphParallelTempering-beta-chart-is-nested-with-a-tag"
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

(* The softmax over edges must survive a table where every move is uphill.

   SGradDescent picks the next edge with probability proportional to Exp[-beta*deltaE]. From
   the complete graph with J = -1 every toggle REMOVES an edge and costs +8, so at beta 200
   every exponent is -1600 and the unshifted form underflowed the whole table to 0.; Total came
   out 0., the normalisation made every weight Indeterminate, and RandomChoice then failed
   outright. Measured on the previous code this run raised ten distinct message types --
   General::munfl, Power::infy, RandomChoice::wghtv, Part::partw and more -- and returned
   through a cascade of errors rather than a choice.

   Subtracting the maximum before exponentiating leaves the distribution untouched (the factor
   cancels) while making the largest weight exactly 1, so Total is never 0. Note the sign: the
   exponent is -beta*deltaE, so it is the DOWNHILL edges -- the ones this is meant to favour --
   that made it large and positive and pushed Exp into arbitrary precision at the other end.

   Asserted as message-free rather than merely non-crashing, because the old failure announced
   itself loudly and a shape check alone passed straight through it. The benign case is already
   covered by the two smoke tests above, which pin that well-conditioned runs are unchanged. *)
VerificationTest[
    Module[{K6 = Normal[AdjacencyMatrix[CompleteGraph[6]]], msgs, r},
        SeedRandom[302];
        msgs = messageArguments["General::munfl",
            r = Block[{Print = Null &},
                ECGrav`SGradDescent[K6, h[-1.0, 0.0], dH[-1.0, 0.0], 200.0, 10]]];
        {Length[msgs],
         Head[r], NumberQ[r["minE"]], r["minE"] == -60.,
         (* it really did stay in the all-uphill regime the whole way *)
         r["LastState"] === K6,
         (* and the case is live: every exponent is far past the representable range *)
         200.0*ECGrav`delHIsing[K6, -1.0, 0.0, 1, 2] > 709}],
    {0, Association, True, True, True, True},
    TestID -> "SGradDescent-softmax-survives-all-uphill"
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

(* Regression: the energy-callback probe used the wrong ARITY in every external-field driver.

   In an external-field run the callbacks' parameters are a row of the field table, not hparams:
   the driver does Apply[hamiltonian,externalFieldTable[[i]]] and hparams is empty, so the
   callback is ultimately called as h[graph,field]. The probe nonetheless evaluated h[graph],
   which for a hamiltonian written the way those drivers require -- h[am_List,f_Real], passed as
   h[] -- comes back unevaluated. Every external-field GraphMultiHistogram, GraphCTLSchedule and
   external-field GraphParallelTempering overload therefore returned $Failed with ECGrav::badham
   before doing any work, rejecting a correct hamiltonian. Beta-tempering overloads were never
   affected: there hparams really is the parameter list.

   Pinned at the probe rather than through a driver because the drivers equilibrate before they
   return anything -- minutes -- and the probe is where the contract lives. Assertions 3 and 4
   are the bug itself and keep this from passing vacuously; 5 confirms the fix did not defang the
   check. 6 records a known limit rather than an aspiration: a slot-based pure Function ignores
   extra arguments, so it satisfies THIS probe and is caught one level in, by the inner driver's
   own probe on the applied form. *)
VerificationTest[
    Module[{hF, dHF, am = Normal[AdjacencyMatrix[CompleteGraph[4]]]},
        hF[a_List, f_Real] := ECGrav`HIsing[a, -1.0, 0.0] + f*Total[a, 2];
        dHF[a_List, f_Real, i_Integer, j_Integer] := 1.0;
        Quiet[{
            ECGrav`Private`HamiltonianUsableQ[am, hF, -0.5],
            ECGrav`Private`DelHUsableQ[am, dHF, -0.5],
            ECGrav`Private`HamiltonianUsableQ[am, hF],
            ECGrav`Private`DelHUsableQ[am, dHF],
            ECGrav`Private`HamiltonianUsableQ[am, undefinedHamHead, -0.5],
            ECGrav`Private`HamiltonianUsableQ[am, ECGrav`HIsing[#, -1.0, 0.0] &, -0.5]}]],
    {True, True, False, False, False, True},
    TestID -> "EnergyCallbackProbe-external-field-arity"
];

(* Regression: the SAME arity bug survived 1.8.1 in the two schedule-driven
   GraphParallelTempering overloads (MCSims.wl:7880 with delH, 8178 without).

   Those drivers take the callback parameters from a replica's field LABEL --
   Apply[hamiltonian, replicaSwapEdgesLabels[[2,Key[i]]]] -- which is a third calling
   convention alongside hparams (beta tempering) and the field table (GraphMultiHistogram,
   GraphCTLSchedule). 1.8.1 fixed the other two and left these probing hparams, which is
   empty when the hamiltonian is passed as h[]. Feeding a schedule straight from
   GraphCTLSchedule -- GraphParallelTempering[seed, bt, h[], dH[], sched[[2;;3]], ...] --
   therefore returned $Failed with ECGrav::badham before doing any work.

   Pinned at the CALL SITE, not at the probe, which is exactly what
   EnergyCallbackProbe-external-field-arity above does NOT cover: that test proves the probe
   behaves correctly given the right arguments, and this one proves the driver hands it the
   right arguments. Stubbing the two probes lets the gate fail deliberately, so the driver
   returns before equilibrating -- seconds instead of minutes. The assertion is that the
   field value -0.5 reaches the probe at all; with the bug the probe sees only {graph, head}
   and both entries come back False, so this cannot pass vacuously. *)
VerificationTest[
    Module[{hF, dHF, seed = Normal[AdjacencyMatrix[CompleteGraph[5]]], sched, withDelH, noDelH},
        hF[a_List, f_Real] := ECGrav`HIsing[a, -1.0, 0.0] + f*Total[a, 2];
        dHF[a_List, f_Real, i_Integer, j_Integer] := 1.0;
        sched = {{1 \[UndirectedEdge] 2}, <|1 -> {-0.5}, 2 -> {0.5}|>};

        (* ham -> True so delH is reached; delH -> False so the gate fails and the driver
           returns immediately *)
        withDelH = Reap[Block[{
              ECGrav`Private`HamiltonianUsableQ = (Sow[{##}]; True) &,
              ECGrav`Private`DelHUsableQ = (Sow[{##}]; False) &},
            Quiet[ECGrav`GraphParallelTempering[seed, 1.0, hF[], dHF[], sched,
              {ECGrav`EulerChi[#] &}, {ECGrav`EulerChi[#] &}, 1, 0, 1000.0]]]][[2]];

        (* the no-delH overload probes only the hamiltonian, so that stub must fail the gate *)
        noDelH = Reap[Block[{
              ECGrav`Private`HamiltonianUsableQ = (Sow[{##}]; False) &},
            Quiet[ECGrav`GraphParallelTempering[seed, 1.0, hF[], sched,
              {ECGrav`EulerChi[#] &}, {ECGrav`EulerChi[#] &}, 1, 0, 1000.0]]]][[2]];

        {MemberQ[Flatten[withDelH], -0.5],   (* hamiltonian probe got the label *)
         Count[Flatten[withDelH], -0.5] >= 2, (* delH probe got it too          *)
         MemberQ[Flatten[noDelH], -0.5]}],
    {True, True, True},
    TestID -> "EnergyCallbackProbe-schedule-driven-uses-replica-label"
];

(* ComputeMinusBetaTimesFreeEnergy: the MBAR/WHAM self-consistency solve, rewritten onto
   MBARSolveFreeEnergies. The old form recomputed the inner sum for every target state, giving
   O(K^3 N) interpreted work per iteration; the term carrying the target is beta_k*E_s, which has
   no j in it, so it factors out of the inner LogSumExp and the outer loop collapses. On a real
   11-beta, 999-sample run this took 218 s on eleven kernels and now takes 0.38 s on one.

   Leg 1 pins the fast form against the LITERAL DEFINITION -- the nested sum written out
   independently here, iterated a fixed number of times rather than to a tolerance -- so it tests
   the algebra rather than a refactor of itself. Legs 2 and 3 are the same data expressed as
   scalar keys and as one-element vector keys, which must give identical answers: that is what
   would catch a wrong A in one overload but not the other -- but only away from bt = 1, where
   dropping the factor is the identity and an earlier draft of this test passed with that exact
   break. Legs 3 and 4 are the vacuity guards on bt itself. Leg 5 pins the mean-centred gauge.

   The convergence test changed with this rewrite and deliberately so: the old
   Sqrt[Sum[(1 - F_k/Fnew_k)^2]] divides by a quantity whose zero point is pure gauge, so it read
   Indeterminate in the F_1 = 0 gauge on a converged fixture, and it stopped about 12x short of
   the accuracy it appeared to promise. *)
VerificationTest[
    Module[{betas = {0.3, 0.9, 1.7}, dat, ref, plain, got2, got3, got2bt, got3bt, ctr},
        SeedRandom[314];
        dat = Association@Table[b -> RandomVariate[NormalDistribution[-2.0*b, 1.5], 40], {b, betas}];
        (* the definition, written out, iterated a fixed 400 times *)
        plain = Module[{K = 3, logn, F, nF, ntab},
            ntab = Length /@ Values[dat]; logn = Log[1.0*ntab];
            F = ConstantArray[1.0, K];
            Do[nF = Table[
                   ECGrav`LogSumExp[Flatten[Table[Table[
                     -ECGrav`LogSumExp[Table[
                        logn[[j]] - F[[j]] + (betas[[k]] - betas[[j]])*dat[[Key[betas[[i]]], s]],
                        {j, K}]],
                     {s, ntab[[i]]}], {i, K}]]], {k, K}];
               F = nF - Mean[nF], {400}];
            F];
        ref = Values[ECGrav`ComputeMinusBetaTimesFreeEnergy[dat]];
        (* Same 1-D data as scalar-keyed and as vector-keyed external fields, at bt != 1.
           bt = 1 would make legs 2 and 4 vacuous: dropping the bt factor from one overload is
           the identity there, and an earlier draft of this test passed with exactly that break. *)
        SeedRandom[271];
        Module[{hs = {-0.7, 0.1, 0.9}, o, vec},
            o = Association@Table[h -> RandomVariate[NormalDistribution[-3.0*h, 2.0], 40], {h, hs}];
            vec = Association@Table[{h} -> Transpose[{o[h]}], {h, hs}];
            got2 = Values[ECGrav`ComputeMinusBetaTimesFreeEnergy[o, 1.7]];
            got3 = Values[ECGrav`ComputeMinusBetaTimesFreeEnergy[vec, 1.7]];
            got2bt = Values[ECGrav`ComputeMinusBetaTimesFreeEnergy[o, 1.0]];
            got3bt = Values[ECGrav`ComputeMinusBetaTimesFreeEnergy[vec, 1.0]]];
        ctr = Values[ECGrav`ComputeMinusBetaTimesFreeEnergy[dat]];
        {Max[Abs[(ref - Mean[ref]) - (plain - Mean[plain])]] < 10.^-6,   (* 1: matches the definition   *)
         Max[Abs[got2 - got3]] < 10.^-12,                                (* 2: the overloads agree      *)
         Max[Abs[got2 - got2bt]] > 10.^-3,                               (* 3: scalar overload uses bt  *)
         Max[Abs[got3 - got3bt]] > 10.^-3,                               (* 4: vector overload uses bt  *)
         Abs[Mean[ctr]] < 10.^-12}],                                     (* 5: mean-centred gauge       *)
    {True, True, True, True, True},
    TestID -> "ComputeMinusBetaTimesFreeEnergy-matches-definition"
];

(* The MBAR solve used to return silently when it hit its iteration cap without converging, so a
   caller got free energies worse than the tolerance implied with no indication -- and every
   quantity reweighted from them inherited that quietly. It now warns.

   Leg 1 is the warning firing on a deliberately starved budget. Leg 2 is the negative control and
   the thing that keeps this honest: the same data with a real budget must be SILENT, so the test
   fails if the condition is ever inverted or the threshold made vacuous. Leg 3 pins that the free
   energies are still returned rather than swallowed -- this is a warning, not a failure. *)
VerificationTest[
    Module[{betas = {0.3, 0.9, 1.7}, dat, A, Ev, logn, res, warned, quiet},
        SeedRandom[9];
        dat = Association@Table[b -> RandomVariate[NormalDistribution[-2.0*b, 1.5], 40], {b, betas}];
        Ev = N[Join @@ Values[dat]]; logn = Log[N[Length /@ Values[dat]]];
        A = Outer[Times, Ev, betas];
        warned = Quiet[
            Check[res = ECGrav`Private`MBARSolveFreeEnergies[A, logn, 10.^-8, 5]; False,
                  True, ECGrav`ComputeMinusBetaTimesFreeEnergy::noconv],
            ECGrav`ComputeMinusBetaTimesFreeEnergy::noconv];
        quiet = Quiet[
            Check[ECGrav`Private`MBARSolveFreeEnergies[A, logn, 10.^-8, 30000]; True,
                  False, ECGrav`ComputeMinusBetaTimesFreeEnergy::noconv],
            ECGrav`ComputeMinusBetaTimesFreeEnergy::noconv];
        {warned,                                                (* 1: warns when starved      *)
         quiet,                                                 (* 2: silent when it converges *)
         Length[res] === 3 && VectorQ[First[res], NumericQ]}],   (* 3: still returns the F      *)
    {True, True, True},
    TestID -> "MBARSolveFreeEnergies-warns-on-nonconvergence"
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
                        ECGrav`$ECGravMaxCorrelationSweeps = 333,
                        ECGrav`$ECGravCorrelationRunMultiplier = 9},
            ECGrav`Private`SyncParallelSettings[];
            ParallelEvaluate[{Last /@ OwnValues[ECGrav`$ECGravMaxEquilibriationSweeps],
                              Last /@ OwnValues[ECGrav`$ECGravMaxCorrelationSweeps],
                              Last /@ OwnValues[ECGrav`$ECGravCorrelationRunMultiplier]}]];
        shipped = DistributeDefinitions[ECGrav`$ECGravMaxCorrelationSweeps];
        CloseKernels[]; Quiet[LaunchKernels[]];
        {DeleteDuplicates[pushed], shipped}],
    {{{{777}, {333}, {9}}}, {}},
    TestID -> "SyncParallelSettings-reaches-subkernels"
];
