# ECGrav Roadmap

Working plan for turning ECGrav into a well-documented, openly shared Wolfram
Language paclet. Checkboxes track progress; phases are roughly ordered by
priority. This is a living document — update it as items land.

_Status snapshot: 2026-07-17._

**Decisions made**
- License: **MIT** ✅
- Citation: no associated paper yet — revisit when one exists.
- Remote: an existing GitHub repo holds an older version; plan is to preserve
  the old state (tag/branch) and make this cleaned-up history canonical.

---

## Phase 0 — Foundation ✅ (done)

- [x] Initialize git repository with a baseline commit
- [x] `.gitignore` — ignore `.DS_Store`, local Claude settings, and the expanded
      `build/ECGrav/` tree while keeping the distributable `.paclet` tracked
- [x] Remove debug `Print` scaffolding; convert warning `Print`s to `Message[]`
- [x] Route simulation progress messages through `PrintTemporary`
- [x] Restore the diagnostic-plot `Print`s that had been stripped
      (`GraphEquilibriate`, `GraphComputeCorrelationTime`, `GraphMultiHistogram`,
      `GraphCTLSchedule`, `ConstrainedProbConjugateField`, density plots)
- [x] Remove dead `LandauFreeEnergy` declaration (no implementation existed)
- [x] Fix `N` shadowing in `GraphMetropolis` (`N` → `vCount`)
- [x] Add `.wlt` regression suite under `Tests/` (63 tests, all passing)
- [x] Fix `Branchedness` `Null` leak on facet-list input
- [x] Fix `LowEnergyStates` divide-by-zero when no parallel kernels are launched
- [x] Add MIT `LICENSE`

## Phase 1 — GitHub foundations

- [x] Create/connect the GitHub remote and push (old versions preserved as tags `v1.0.0`, `v1.1.0`; released `v1.2.0`)
- [x] Top-level `README.md`: what Emergent Combinatorial Gravity is, install,
      a short quickstart, current status and limitations
- [x] `.gitattributes` marking `*.paclet` (and other binaries) as binary
- [x] `CHANGELOG.md` starting at 1.2.0
- [ ] `CITATION.cff` — deferred until there is a paper (software-only cite optional)

## Phase 2 — Documentation completeness

- [ ] Add usage messages for the 5 undocumented `*Conn` functions
- [ ] Author reference pages for the ~52 public functions that lack them
      (essentially the whole MCSims domain: Hamiltonians, MC drivers, tempering,
      free-energy tools)
- [ ] Add at least one **Guide** page mapping the API by theme
      (`Documentation/English/Guides/` is currently empty)
- [ ] Add a getting-started **Tutorial**: build a graph → run Metropolis →
      read observables → plot (`Documentation/English/Tutorials/` is empty)
- [ ] Make the README quickstart mirror the tutorial

## Phase 3 — Quality & CI

- [ ] Resolve or document known bugs (see below)
- [ ] Expand test coverage beyond the current ~40 / 116 symbols; gaps include
      parallel tempering, simulated annealing, gradient descent, the random-complex
      generators, and the free-energy extrapolation family
- [ ] Add CI (GitHub Actions) running the `.wlt` suite via Wolfram Engine on push/PR
- [ ] Address the ~119 `ReturnAmbiguous` warnings (bare `Return` inside `If`/`With`)
- [ ] Decide on the subpackage `BeginPackage` convention: every symbol is
      hand-qualified as `ECGrav\`Foo` because the subpackages open
      `BeginPackage["ECGrav\`MCSims\`"]` with no exports — fix or document it so a
      new function added without the qualifier doesn't silently vanish

## Phase 4 — Release & distribution

- [ ] Adopt semantic-versioning discipline and tag releases
- [ ] Attach the built `.paclet` to GitHub releases
- [ ] Decide distribution channel: Wolfram Paclet Repository submission vs.
      install-from-GitHub-release
- [ ] Turn the manual steps in `BuildingAndInstallingPaclets.nb` into a
      reproducible build script

## Phase 5 — Scientific documentation (optional, high value for research code)

- [ ] Background/theory note (the physics and mathematics of the model)
- [ ] Example-notebook gallery
- [ ] Links to associated paper(s) once available

---

## Known issues / bugs

Found while building the test suite; **not** yet fixed, and deliberately not
asserted as passing behaviour in the suite:

- **`CombinatorialSphereQ[torus]` returns `True`.** A torus (χ = 0) is not a
  sphere. Only the correct cases (tetrahedron, octahedron) are currently asserted.
- **`ExactExpectationValue` emits a stray `Part::partd`** while probing which
  overload its observable argument matches. The result is correct; the message is
  declared as expected in the test. Worth silencing at the source.

See `Tests/README.md` for the full list of test caveats and design notes.
