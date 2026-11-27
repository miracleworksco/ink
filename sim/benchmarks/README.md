# ext's benchmarks

Simulation (strength, stiffness, vibration, heat), sizing by it, and
topology optimization against hand calculations and reference data,
as the kernel's benchmarks/ are. Each file's claims are checked by the
model rule (ink-ext/tests/cli/models.rs: the quick ones by `cargo
test`), its key answers that are not claims by its test (`bNN_...` in
ink-ext/tests/cli/benchmarks.rs). `tools/` holds independent references (truss.py,
blade_beam.py, top88_mbb.py).

`holds/` holds claims on parts held by fixes (337 cases, 408 claims: the
statics and dynamics hunts', the holds reviews', a plate), each against its
exact answer: `python3 benchmarks/holds/run.py INK -j 3` counts right, wrong,
not settled and refused, against the totals its header records for main and
the holds branch (about 50 minutes at one case at a time).

`thin/` holds claims on material thinner than the cells (fins, webs,
tabs, ligaments, split cells), each with its true answer and how that is
known: `python3 benchmarks/thin/run.py INK` counts the false holds
and false fails a build of ink gives on them, each claim also asked the
other way (`max-stress >` for `<`) on the same solves. Every case, the
slow ones (`--slow`) too, must answer as it says, its declared `known:`
wrong answers apart (about 20 minutes: the area's slow set);
`--slow --only NAME` reruns one.

`memory/`: a solve's memory estimate (what it warns by, before it
allocates, where it may not fit) against the peak it reaches, on the models
there at several cell sizes and every study shipped here, in examples
and in any models given: `python3 benchmarks/memory/run.py INK
[MODEL...]`.

Topology on the MBB beam of Andreassen et al.'s 88-line code (60 × 20
cells, volume 0.5, filter 2.4): the layout settles in about 60 iterations
at 1.52 times the solid domain's compliance; the 88-line code, ported
(`tools/top88_mbb.py`), gives 1.51 with ink's projection and update, 1.60
with its own Heaviside filter (ink-ext/tests/topology.rs checks 3%).

| # | benchmark | checked |
|---|---|---|
| 03 | lofted, twisted NACA blade under a tip load | tip deflection under 5 mm (FEA) |
| 05 | L bracket loaded and held by named faces | each face reference selects material; safety factor > 2 |
| 08 | helical compression spring under 50 N | deflection within 10% of G d^4 / 8 D^3 n (end-corrected); p99 stress within 25% of Wahl |
| 12 | Pratt truss bridge under a midspan load | chord stresses within 10%, deflection within 20% of the pin-jointed calc |
| 13 | pressure vessel ring vs Lame | radial growth within 5%, von Mises within 10%; sector and true pressure |
| 15 | topology-optimized cantilever bracket | layout settles; volume target; stiffer than a same-volume web; watertight |
| 16 | pin-fin heat sink | temperature rise within 10% of fin theory; 10 W heat balance |
| 17 | transient heat-up of a copper block | T(600 s) within 2%, 80 C crossing within 3% of lumped capacitance |
| 21 | gyroid sheet lattice block | volume fraction within 5%; one solid, no cavities; E*/Es in [0.1, 0.35] |
| 23 | two-material bar | series stiffness, mass and COM |
| 24 | cantilever snap-fit hook | strain, seated clearance, retaining gap, insertion with controls |
| 25 | modal: cantilever plate and a bracket | bending 1 within 3% of Euler-Bernoulli; torsion found; frequency floors |
| 28 | stepped shaft with keyway and hub | bending, torsion and von Mises within 5% of M c / I and T r / J; shoulder peaks within 25% of Peterson; deflection within 10%; fits proven; tight hub disproven |
| 30 | PETG two-leaf flexure guide | stiffness, surface stress and first mode vs beam theory |
| 36 | pin, roller and bearing supports | deflections and reactions vs beam theory; mechanisms reported |
