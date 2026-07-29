# ECGrav — Emergent Combinatorial Gravity

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Version](https://img.shields.io/badge/version-1.3.1-blue.svg)
![Wolfram Language](https://img.shields.io/badge/Wolfram%20Language-14.2%2B-red.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21684432.svg)](https://doi.org/10.5281/zenodo.21684432)

A Wolfram Language paclet for studying models in which geometry — and ultimately
gravity — **emerges from the statistical mechanics of combinatorial structures**
(graphs and simplicial complexes).

ECGrav provides two complementary toolkits:

- **Combinatorics of complexes and graphs** — build and analyze pure simplicial
  complexes and graphs: dimensions, Euler characteristic, homology/holes,
  links and spheres, manifold and sphere predicates, automorphism groups, and
  isomorphism classes.
- **Monte Carlo statistical mechanics** — sample ensembles of graphs/complexes
  under energy functionals (Hamiltonians) that reward geometric, manifold-like
  configurations, and extract thermodynamic observables via Metropolis MC,
  parallel tempering, multi-histogram reweighting, and exact enumeration.

> **Status:** 1.3.0. Actively developed research code. See [ROADMAP.md](ROADMAP.md)
> for what's planned and [Tests/README.md](Tests/README.md) for known issues.

For the mathematics of the model — the energy functionals, observables, and
statistical mechanics, plus background on the emergent-gravity motivation — see
**[Theory.md](Theory.md)**.

---

## Requirements

- Wolfram Language / Mathematica **14.2 or later**
- Multiple parallel kernels are recommended — many Monte Carlo routines use
  `ParallelTable`, and some (e.g. `LowEnergyStates`, parallel tempering) expect
  `LaunchKernels[]` to have been called.

## Installation

### Option 1 — load the source directly (no build)

```wolfram
PacletDirectoryLoad["/path/to/ECGrav"];   (* the cloned repo directory *)
Needs["ECGrav`"];
```

### Option 2 — install a released paclet

Download `ECGrav-<version>.paclet` from the
[GitHub Releases](https://github.com/kassabetre/ECGrav/releases) page, then:

```wolfram
PacletInstall["/path/to/ECGrav-1.3.0.paclet"];
Needs["ECGrav`"];
```

### Option 3 — build the paclet from source

Build the distributable paclet with the bundled script (it regenerates the docs and
verifies the archive), then install what it produced:

```bash
wolframscript -file build.wls
```

```wolfram
PacletInstall["/path/to/ECGrav/build/ECGrav-1.3.0.paclet", ForceVersionInstall -> True];
Needs["ECGrav`"];
```

For the Monte Carlo routines, also start parallel kernels:

```wolfram
LaunchKernels[];
ParallelNeeds["ECGrav`"];
```

## Quick start

### Combinatorics on a simplicial complex

```wolfram
Needs["ECGrav`"];

(* An octahedron: 8 triangular facets on 6 vertices — a triangulated 2-sphere *)
octahedron = {{1,2,3},{1,2,4},{1,3,5},{1,4,5},{6,2,3},{6,2,4},{6,3,5},{6,4,5}};

EulerChi[octahedron]              (* 2        — Euler characteristic of S^2      *)
FVector[octahedron]               (* {6,12,8} — 6 vertices, 12 edges, 8 faces    *)
CountHoles[octahedron]            (* {1,0,1}  — Betti numbers b0, b1, b2          *)
CombinatorialSphereQ[octahedron]  (* True                                        *)
```

### A Metropolis Monte Carlo run

```wolfram
Needs["ECGrav`"];

seed = Normal[AdjacencyMatrix[CycleGraph[10]]];   (* seed graph as an adjacency matrix *)
H    = HIsing[#, -1.0, 0.0] &;                     (* Ising energy; J = -1 rewards edges *)

(* Run 50 Metropolis sweeps at inverse temperature 0.5, measuring the average degree. *)
{finalGraph, series} = GraphMetropolis[seed, 0.5, H, {AvgDeg}, 50];

Dimensions[series]                 (* {50, 3} — one {sweep, energy, observable} row per sweep *)
ListLinePlot[series[[All, 3]]]     (* average degree climbing as the graph densifies *)
```

`GraphMetropolis` returns the last graph visited and the `{sweep, energy, observable}`
table along the trajectory — this mirrors the [Getting Started tutorial](Documentation/English/Tutorials/GettingStartedWithECGrav.nb).
Higher-level drivers (`GraphEquilibriate`,
`GraphComputeCorrelationTime`, `GraphParallelTempering`, `GraphMultiHistogram`)
build on this and emit diagnostic plots as they run.

## What's in the package

ECGrav exports 113 public functions across two subpackages. A selection:

### Simplicial complexes & graphs

| Area | Functions (examples) |
| --- | --- |
| Representations & conversions | `FacetIncidenceMatrix`, `FacetAdjacencyMatrix`, `GraphFromCliques`, `CliquesFromFacetIncidence`, `ComplexFromFacetLabeledVertexList` |
| Local structure | `Sph`, `Bll`, `Lnk`, `Str`, `Deg`, `HyperDeg`, `FacetDeg` |
| Dimension & topology | `AvgKDim`, `SpectralDim`, `HausdorffDim`, `EulerChi`, `CountHoles`, `FVector`, `GVolume` |
| Manifold predicates | `CombinatorialManifoldQ`, `ClosedCombinatorialManifoldQ`, `CombinatorialSphereQ`, `DSphereQ`, `ClosedDGraphQ`, `DGraphQ`, `OrientableCombinatorialManifoldQ` |
| Symmetry & enumeration | `SimplicialComplexAutomorphismGroupOrder`, `IsomorphicSimplicialComplexQ`, `ChooseNonIsomorphicSimplicialComplexes`, `NumPureComplexes`, `RankComb`, `UnrankComb` |
| Random complexes | `RandomVertexLabeledPureSimplicialComplex`, `RandomFacetLabeledPureSimplicialComplex`, `RandomUnlabeledPseudoManifold` |

### Monte Carlo & statistical mechanics

| Area | Functions (examples) |
| --- | --- |
| Hamiltonians | `HIsing`, `HWeightedFaceCounts`, `HEdgeDeg`, `HLaplacian`, `H1dCombManifold`, `H2dCombManifold`, `H2dPseudoCombManifold` |
| Ground-state search | `LowEnergyStates`, `GradDescent`, `SGradDescent`, `SimulatedAnnealing` |
| Sampling | `GraphMetropolis`, `GraphSweepReplica`, `GraphEquilibriate`, `GraphParallelTempering` |
| Reweighting & schedules | `GraphMultiHistogram`, `GraphCTLSchedule`, `GraphCEITempSchedule`, `ComputeMinusBetaTimesFreeEnergy`, `ExtrapolatedExpectationValue` |
| Exact & analysis | `ExactExpectationValue`, `CorrelationTime`, `SpecificHeat`, `Susceptibility`, `ErrorBootstrap`, `LogSumExp` |

Use the built-in help for details on any symbol:

```wolfram
?ECGrav`*        (* list every public symbol *)
?GraphMetropolis (* usage message for one function *)
```

Reference documentation pages are available for every public function in the Wolfram
Documentation Center, along with an **ECGrav** guide (API by theme) and a
**Getting Started with ECGrav** tutorial.

## Tests

The `Tests/` directory holds a `.wlt` regression suite (87 tests) that loads the
package **from source**:

```wolfram
Needs["MUnit`"];
TestReport["/path/to/ECGrav/Tests"]
```

See [Tests/README.md](Tests/README.md) for the suite's design, coverage, and a
list of known issues.

## Versioning

This repository follows semantic-ish versioning. The `main` branch is the
current release (1.3.1). Earlier releases are preserved as git tags:

- `v1.3.1` — current: documentation & tooling — model specification (`Theory.md`),
  citation metadata + Zenodo DOI, and a reproducible build script
- `v1.3.0` — Phase 3 quality & CI: bug fixes, expanded test suite (87 tests),
  GitHub Actions CI, and the collapse to a single `ECGrav` public context
- `v1.2.0` — cleanup, regression test suite, bug fixes, MIT license
- `v1.1.0`, `v1.0.0` — the earlier codebase (preserved for reference)

## License

Released under the [MIT License](LICENSE).

## Authors

- **Kassahun Betre** — author
- **Khai Luong** — contributor

If you use ECGrav in academic work, please cite this repository (a companion
paper is planned).
