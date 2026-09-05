# ECGrav — Emergent Combinatorial Gravity

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Version](https://img.shields.io/badge/version-1.14.0-blue.svg)
![Wolfram Language](https://img.shields.io/badge/Wolfram%20Language-15.0%2B-red.svg)
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

> **Status:** 1.14.0. Actively developed research code. See [ROADMAP.md](ROADMAP.md)
> for what's planned and [Tests/README.md](Tests/README.md) for known issues.

For the mathematics of the model — the energy functionals, observables, and
statistical mechanics, plus background on the emergent-gravity motivation — see
**[Theory.md](Theory.md)**. The exact samplers and counters over pure complexes have
their own specifications: **[FacetLabeledCount.md](FacetLabeledCount.md)**,
**[FacetLabeledSampler.md](FacetLabeledSampler.md)**,
**[UnlabeledCount.md](UnlabeledCount.md)**, and
**[UnlabeledSampler.md](UnlabeledSampler.md)**.

Running parallel tempering over inverse temperature **and** an external field at the same time
needs no new machinery — the shipped multi-field path does it once the hamiltonian is written in
homogeneous form, with beta as one of the tempered couplings.
**[HomogeneousHamiltonian.md](HomogeneousHamiltonian.md)** is the recipe, its justification, and
the units caveat that comes with it.

The schedule builder's cost is dominated by the MBAR reweighting it interpolates the thermodynamic
metric from, and that cost grows as the square of the bootstrap table — a few seconds at one
external field, hours at two. **[MBARWeights.md](MBARWeights.md)** specifies the rewrite that
reduces every MBAR quantity to a ratio of one shared weight vector, with the algebra, the
reference implementation, and the measurements.

---

## Requirements

- Wolfram Language / Mathematica **15.0 or later**
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
PacletInstall["/path/to/ECGrav-1.14.0.paclet"];
Needs["ECGrav`"];
```

### Option 3 — build the paclet from source

Build the distributable paclet with the bundled script (it regenerates the docs and
verifies the archive), then install what it produced:

```bash
wolframscript -file build.wls
```

```wolfram
PacletInstall["/path/to/ECGrav/build/ECGrav-1.14.0.paclet", ForceVersionInstall -> True];
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

ECGrav exports 116 public functions across two subpackages. A selection:

### Simplicial complexes & graphs

| Area | Functions (examples) |
| --- | --- |
| Representations & conversions | `FacetIncidenceMatrix`, `FacetAdjacencyMatrix`, `GraphFromCliques`, `CliquesFromFacetIncidence`, `ComplexFromFacetLabeledVertexList` |
| Local structure | `Sph`, `Bll`, `Lnk`, `Str`, `Deg`, `HyperDeg`, `FacetDeg` |
| Dimension & topology | `AvgKDim`, `SpectralDim`, `HausdorffDim`, `EulerChi`, `CountHoles`, `FVector`, `GVolume` |
| Manifold predicates | `CombinatorialManifoldQ`, `ClosedCombinatorialManifoldQ`, `CombinatorialSphereQ`, `DSphereQ`, `ClosedDGraphQ`, `DGraphQ`, `OrientableCombinatorialManifoldQ` |
| Symmetry & enumeration | `SimplicialComplexAutomorphismGroupOrder`, `IsomorphicSimplicialComplexQ`, `ChooseNonIsomorphicSimplicialComplexes`, `NumVertexLabeledPureComplexes`, `NumFacetLabeledPureComplexes`, `NumUnlabeledPureComplexes`, `RankComb`, `UnrankComb` |
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

The `Tests/` directory holds a `.wlt` regression suite (188 tests) that loads the
package **from source**:

```wolfram
Needs["MUnit`"];
TestReport["/path/to/ECGrav/Tests"]
```

See [Tests/README.md](Tests/README.md) for the suite's design, coverage, and a
list of known issues.

## Versioning

This repository follows semantic-ish versioning. The `main` branch is the
current release (1.14.0). Earlier releases are preserved as git tags:

- `v1.14.0` — current: two observables that were silently wrong when misused are guarded at
  the source. `CountHoles[c, k]` now returns the k-th Betti number `b_k` rather than `b_{k-1}`,
  so `CountHoles[c, 1]` is the loop count and no longer the component count — **re-run any saved
  analysis that calls the two-argument form**; and `FractionInLargestComponent` no longer accepts
  a facet list, which it used to read as an adjacency matrix whenever it happened to be square.
  Also 31 new tests, `EnsembleWeights.md`, and `md2pdf.py`
- `v1.13.0` — a run records which configuration is at which slot at every sweep,
  so round trips through the tempering ladder can be counted — with `TemperingMixing.wls`
  to count them and `TemperingPlots.wls` to read continuous surfaces off the same run.
  Beta-schedule chart rows are no longer flattened, so their observables move from
  `chart[[…, 4 ;;]]` to `chart[[…, 4]]`; energy stays at column 3
- `v1.12.1` — the MBAR free-energy solve warns when it does not converge —
  it used to return silently on hitting its iteration cap, leaving free energies worse
  than the tolerance implied with nothing to say so
- `v1.12.0` — the MBAR free-energy solve is rebuilt on one shared matrix and
  its convergence test replaced — the self-consistency iteration was `O(K^3 N)` per pass
  and its stopping rule divided by a quantity whose zero point is pure gauge, stopping
  about twelve times short of the accuracy it appeared to promise
- `v1.11.1` — `GraphParallelTempering` rejected correct hamiltonians when handed a
  schedule from `GraphCTLSchedule` — its two schedule-driven overloads probed the energy
  callbacks with `hparams`, which is empty, where the driver itself calls them with a
  replica's field label
- `v1.11.0` — the constant-thermodynamic-length metric is rebuilt from one
  shared MBAR weight vector — the free energy cancels out of every expectation value, all
  interpolation grids share one weight vector, and the sum over bootstrap field points is
  hoisted out of the grid entirely, taking a two-component schedule build from over seven
  hours to a fraction of a second with the numerics unchanged to 2e-14
- `v1.10.1` — four defects in the constant-thermodynamic-length schedule
  builder, all of them invisible with a single external field — a replica's field-vector
  components could be permuted, the thermodynamic metric was built from second moments
  rather than covariances and had its off-diagonals clipped to zero, and replicas could be
  placed on a metric extrapolated into regions the bootstrap run never sampled
- `v1.10.0` — the correlation-time run length is fully settable —
  `$ECGravCorrelationRunMultiplier` replaces the multiplier that was hard-coded to 5,
  which was the half of `Min[5*eqlT, $ECGravMaxCorrelationSweeps]` that binds on every
  short run, so raising the cap alone did nothing
- `v1.9.1` — the empty complex had no Euler characteristic — `EulerChi[{}]`
  returned `$Failed` rather than 0, which any code taking the link of an edge hits as
  soon as the edge's endpoints share no neighbour
- `v1.9.0` — the graph correlation time caught up with its complex-space
  twin — the lag table is walked and stopped at the turnover rather than tabulated
  in full, which was quadratic in the run length; three defects fixed along the way,
  including a `corrT` that came back unevaluated on short runs; and every
  correlation-time result now reports `corrTMeasured`, separating a measured
  interval from the floor it falls back to
- `v1.8.1` — the energy-callback check rejected correct hamiltonians in
  every external-field driver — it probed the callbacks at the arity they have when
  the parameters come from `hparams`, not the one they have when the driver takes
  them from the field table
- `v1.8.0` — the facet-labeled sampler's draw grouped by completion
  key — 2–103× faster with the same distribution and the same public interface,
  because the number of distinct candidate *weights* per facet is a small constant
  where the number of candidates is not
- `v1.7.1` — the ground-state search reported minima that were not
  minima — a sweep never weighed the state it started from, so the state a replica
  swap handed it was invisible, and reported energies are now the hamiltonian's own
  values rather than incremental accumulations
- `v1.7.0` — equilibriation and correlation time, audited and rebuilt — the
  convergence test now compares like with like, watches two observables, and says so
  when neither varied; requires Wolfram Language 15.0+
- `v1.6.0` — fully unlabeled pure complexes — a count of the isomorphism
  classes (`NumUnlabeledPureComplexes`) completing the trio of labellings, and the
  unlabeled random generator made exact instead of rejection-based, up to 14,000×
  faster
- `v1.5.0` — pure-complex counting and random generation — the counting
  functions renamed and joined by a facet-labeled count, three random generators
  fixed where they silently returned malformed output, and the facet-labeled
  sampler made exact instead of rejection-based
- `v1.4.0` — the Monte Carlo hot path computed from adjacency matrices
  instead of `Graph` objects — `HyperDeg`, `delH2dCombManifold` and `EulerChi` are
  8–22× faster with unchanged results, and `GraphEquilibriate` now reports progress
- `v1.3.1` — documentation & tooling — model specification (`Theory.md`),
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
