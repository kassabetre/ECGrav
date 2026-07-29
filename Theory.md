# The ECGrav model — background and specification

> **How to read this document.** Sections 1–4 are a *code-grounded specification*:
> they state precisely what the package computes, with each definition matching the
> implementation in `Kernel/Subpackages/`. Sections 5–6 are an *outline for the
> physics narrative* — the motivation and interpretation that make this "emergent
> combinatorial gravity." Those parts are marked **[Author]** and are intentionally
> left for Kassahun to fill in; the surrounding scaffolding is here to be filled, not
> to assert physics on the author's behalf.

Notation used throughout: $A$ is a symmetric $0/1$ adjacency matrix with zero
diagonal; $\deg(i)=\sum_j A_{ij}$; $c_q$ is the number of $q$-cliques in the graph
(so $c_1$ = vertices, $c_2$ = edges, $c_3$ = triangles, …). A pure simplicial
complex is represented either by its **facet list** (each facet a sorted vertex
list) or, when it is a flag/clique complex, by the graph whose maximal cliques are
the facets (`GraphFromCliques` / `CliquesFromFacetIncidence`).

---

## 1. The configuration space

A configuration is a finite graph (equivalently its clique complex) or a pure
simplicial complex on a labelled vertex set. The package moves between three
equivalent encodings:

- **Adjacency matrix** $A$ — the working representation for the Monte Carlo code.
- **Facet list** — a list of maximal simplices, for the combinatorial-topology code.
- **Clique complex** — the simplicial complex whose $q$-simplices are the
  $q$-cliques of a graph; this is how a graph acquires higher-dimensional topology.

Ensembles are built either over all graphs on $N_0$ vertices (edge-toggling Monte
Carlo) or over pure complexes of fixed purity and facet count (the
`RandomUniform…` and `…MCMC` generators).

---

## 2. Observables (what "geometry" is measured by)

These are the geometric and topological functionals evaluated on a configuration.

### Topology
- **`FVector`** — the $f$-vector $(c_1, c_2, \dots)$, the simplex counts by dimension.
- **`EulerChi`** — the Euler characteristic $\chi = \sum_{q\ge 1} (-1)^{q-1} c_q$
  (alternating sum of the $f$-vector).
- **`CountHoles`** — the Betti numbers $(b_0, b_1, \dots)$ of the clique complex
  (ranks of the homology groups); $b_0$ is the number of connected components,
  $b_1$ the number of independent cycles, and so on.
- **Manifold predicates** — `CombinatorialManifoldQ`, `ClosedCombinatorialManifoldQ`,
  `CombinatorialSphereQ`, and the graph analogues `DGraphQ` / `ClosedDGraphQ` /
  `DSphereQ`. A *closed combinatorial manifold* is connected, has every codimension-1
  face in exactly two facets, and has every vertex link recursively a closed
  combinatorial manifold; a *combinatorial sphere* additionally has the Euler
  characteristic of a sphere, $\chi = 1 + (-1)^{d}$.

### Emergent dimension
Two independent probes of the effective dimension of a configuration, both defined
by a scaling exponent as a function of a "radius" $s$:

- **`SpectralDim`** — from the return probability of an $s$-step walk. With
  $P(s) = \tfrac1N \sum_i (A^s)_{ii}\big/\sum_j (A^s)_{ij}$ the average return
  probability, the spectral dimension is
  $$ d_s(s) \;=\; -\,\frac{2\,\log P(s)}{\log s}, $$
  the discrete form of $P(s)\sim s^{-d_s/2}$.
- **`HausdorffDim`** — from volume growth. With $V(s) = \tfrac1N\sum_i \mathrm{vol}(s,i)$
  the average number of vertices within $s$ steps (or within metric radius $s$, using
  a supplied distance matrix; `GVolume`), the Hausdorff dimension is
  $$ d_H(s) \;=\; \frac{\log V(s)}{\log s}, $$
  the discrete form of $V(s)\sim s^{d_H}$.
- **`AvgKDim`** — an average *local* combinatorial dimension, derived from the clique
  structure around each vertex (see the source for the exact recursion).

### Local structure
`Sph` (the sphere/link neighbourhood of a vertex), `Lnk`, `Str`, `Deg`, and
`HyperDeg` (the number of top simplices containing a given face — e.g. the number of
triangles on an edge). These are the building blocks the manifold-targeting energies
below are written in terms of.

---

## 3. The energy functionals (Hamiltonians)

Each Hamiltonian assigns a real energy $H(A)$ to a configuration; the Monte Carlo
then samples $e^{-\beta H}$ (Section 4). The functionals fall into three groups.

### 3.1 Building block — weighted simplex counts
**`HWeightedFaceCounts`** is a linear combination of clique counts,
$$ H \;=\; \sum_{q=1}^{5} J_q\, c_q, $$
i.e. a chemical potential $J_q$ on each $q$-clique (vertices, edges, triangles,
4- and 5-cliques). The manifold-targeting energies are assembled from it.

### 3.2 Degree-based energies
- **`HIsing`**$(A;J,L)$ — with $\mathrm{wedges} = \sum_i \binom{\deg i}{2}$ the number
  of length-2 paths,
  $$ H \;=\; \frac{J}{2}\sum_i \deg(i)\big(\deg(i)-1\big) \;+\; L\,c_2 \;=\; J\,(\#\text{wedges}) + L\,(\#\text{edges}). $$
- **`HLaplacian`**$(A;J)$ — the Dirichlet energy of the degree field, with
  $L = \mathrm{diag}(\deg) - A$ the graph Laplacian and $d$ the degree vector,
  $$ H \;=\; J\,d^{\mathsf T} L\, d \;=\; J\!\!\sum_{(i,j)\in E}\!\!\big(\deg(i)-\deg(j)\big)^2, $$
  which rewards degree-homogeneous (locally regular) graphs.
- **`HEdgeDeg`**$(A;J,D_1,D_2)$ — a quartic, degree-regularizing energy tuned by
  target parameters $D_1,D_2$,
  $$ H \;=\; \frac{J}{2N}\Big(\operatorname{Tr}A^4 - \operatorname{Tr}\!\big[A^2\,(2D_1 K + 2D_2 I)\big]\Big),\qquad K = \mathbf{1}-I. $$
  (The source's comment lists an additive constant $D_1^2 N(N-1)$; it is dropped in
  code and is immaterial to sampling, since a constant offset does not change the
  ensemble.)

For each degree-based energy there is a matching **`delH…`** giving the exact energy
change under a single edge toggle, so the Metropolis step is $O(N)$ rather than a
full re-evaluation (Section 4).

### 3.3 Manifold-targeting energies — geometry as the ground state
These are the heart of the model: energies engineered so that **manifold-like
complexes are the minima**. Each is a sum of non-negative penalties for violating a
local manifold condition.

- **`H1dCombManifold`**$(A;J)$
  $$ H \;=\; \frac{J}{N_0}\Big(-2c_2 + 2c_3 + \tfrac12\textstyle\sum_i \deg(i)(\deg(i)-1) + N_0 + (\#\text{comp}-1)\Big). $$
  A connected cycle (every vertex degree 2, no triangles) sits at $H=0$; these
  **1-dimensional closed manifolds are the ground states**.

- **`H2dCombManifold`**$(A;J)$
  $$ H \;=\; J\Big(\underbrace{c_2 - 3c_3 + 4c_4}_{\text{simplex balance}}
     + \underbrace{\textstyle\sum_{e\in E}\binom{\mathrm{hdeg}(e)}{2}}_{\text{each edge in 2 triangles}}
     + \underbrace{N_0\,\lvert \#\text{comp}-1\rvert}_{\text{connected}}
     + \underbrace{\textstyle\sum_v \lvert \#\text{comp}(\mathrm{Sph}(v))-1\rvert}_{\text{each link a cycle}}\Big), $$
  where $\mathrm{hdeg}(e)$ is the number of triangles on edge $e$ and $\mathrm{Sph}(v)$
  is the link of vertex $v$. The four terms enforce the local defining conditions of
  a **closed 2-manifold**: each edge in exactly two triangles, and each vertex link a
  single cycle, on a connected complex — so triangulated closed surfaces are the
  ground states.

- **`H2dPseudoCombManifold`**$(A;J)$ — the analogue that targets 2-dimensional
  *pseudo*-manifolds (each ridge in two facets, but links not required to be spheres);
  it replaces the vertex-link term with a 2-path-connectivity term.

`delH2dCombManifold` gives the corresponding single-edge-toggle energy change.

---

## 4. The statistical mechanics

Given a Hamiltonian $H$ and inverse temperature $\beta$, the package samples the
Gibbs ensemble
$$ p(A) \;\propto\; e^{-\beta H(A)}, $$
and estimates observables as ensemble averages $\langle \mathcal{O}\rangle$.

- **Sampling** — `GraphMetropolis` proposes single-edge toggles accepted with the
  Metropolis rule using the incremental $\Delta H$ from the `delH…` functions;
  `GraphEquilibriate` and `GraphComputeCorrelationTime` handle burn-in and
  autocorrelation; `GraphParallelTempering` swaps replicas across a temperature
  ladder to cross barriers; `GraphMultiHistogram` reweights multi-temperature data.
- **Free energy & reweighting** — `ComputeMinusBetaTimesFreeEnergy`,
  `NegativeBetaTimesFreeEnergy`, `ExtrapolatedExpectationValue`, and the
  `…Schedule` helpers estimate $-\beta F$ and extrapolate observables in $\beta$.
- **Thermodynamics** — `InternalEnergy`, `SpecificHeat`, `Susceptibility`,
  `CvOverT`, `ErrorBootstrap` for the standard response functions and their errors.
- **Ground states** — `LowEnergyStates`, `GradDescent`, `SGradDescent`, and
  `SimulatedAnnealing` search for minima of $H$ (the "emergent geometries").
- **Exact checks** — `ExactExpectationValue` enumerates small systems to validate
  the Monte Carlo.

The overall picture: **fix a Hamiltonian whose minima are manifold-like, and the
low-temperature ensemble concentrates on configurations with emergent geometry**;
the observables of Section 2 then measure that geometry.

---

## 5. The physics — emergent combinatorial gravity **[Author]**

> This section is the scientific heart of the note and is left for the author. The
> code specifies *what* is computed; this section should say *why*. Suggested
> structure (prompts, not claims):
>
> - **Motivation.** Why should spacetime geometry be emergent rather than
>   fundamental, and why are graphs / simplicial complexes the right microscopic
>   degrees of freedom?
> - **The gravity dictionary.** What plays the role of the metric, curvature, the
>   Einstein–Hilbert action, and the gravitational path integral in this discrete
>   setting? How do the manifold-targeting Hamiltonians of Section 3.3 encode a
>   preference for geometry — and how (if at all) do their couplings map to physical
>   constants (e.g. a cosmological-constant-like term, Newton's constant)?
> - **Continuum limit / phases.** What is the expected phase structure as $\beta$ and
>   the couplings vary, and in which phase does a smooth $d$-dimensional geometry
>   emerge? What order parameters (Section 2's dimensions, $\chi$, Betti numbers)
>   distinguish the phases?
> - **Results.** Summarize what the simulations show.

## 6. Relation to other approaches & references **[Author]**

> Situate the model relative to the broader literature — the author should confirm
> and cite. The general landscape it sits near includes: causal and Euclidean
> dynamical triangulations; random graphs and random simplicial complexes;
> spectral-dimension and dimensional-reduction studies in quantum gravity; and
> other combinatorial / emergent-geometry programmes. Add the specific references
> (and the associated paper, once it exists) here.

---

*Sections 1–4 track the implementation; if the code changes, update them (and the
`build.wls` / test suite will catch behavioural drift). Sections 5–6 are the
author's to complete.*
