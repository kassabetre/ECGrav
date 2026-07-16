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

(* The seven surfaces in the order used by the notebook's orientability test *)
surfaces = {mobiusStrip1, mobiusStrip2, kleinbottle, projectivePlane,
            torus, tetrahedron, octahedron};

(* Clique-graph views of a couple of the complexes *)
octaG = ECGrav`GraphFromCliques[octahedron];
tetraG = ECGrav`GraphFromCliques[tetrahedron];

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
