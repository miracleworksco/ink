# Benchmarks

Engineering tasks with hand calculations. Each file's header gives the
task, the hand calculation and what must hold. Each is checked two ways:

- its claims (masses, volumes, clearances: the hand numbers a claim can
  state) by the model rule (CONTRIBUTING.md, Test); a control is a claim
  meant to fail, and "disproven" below means it was shown false;
- its other answers (a solve's probes, interference, claims at other
  params or moments, estimates, meshes) by its test, `bNN_...` in
  crates/ink-cli/tests/cli/benchmarks.rs.

`cargo test` runs the quick ones. The slow ones (`SLOW` in
crates/ink-cli/tests/cli/models.rs) are validation:

```sh
benchmarks/run.sh            # every slow benchmark: its test, then its claims
benchmarks/run.sh 10 22      # those benchmarks, slow or not
```

`inputs/` holds two STLs: 22's is `ink mesh` of tools/22-part-source.ink;
33's is tools/33-gen-knob.py at n 240 (b33 regenerates and compares).
`tools/` holds independent references: hand calculations in Python, input
generators, and 32's solve at sampled moments.

Numbers are not reused. 14 (column buckling) is gone. These need
packages this release doesn't ship, and live with them: 3 (twisted
blade), 5 (bracket intent), 8 (helical spring), 12 (truss bridge), 13 (pressure
vessel), 15 (topology bracket), 16 (pin-fin heat sink), 17 (transient
heat-up), 21 (gyroid lattice), 23 (multi-material bar), 24 (snap-fit), 25 (modal
plate), 27 (cooling channel), 28 (stepped shaft), 29 (I-beam against
stiffness), 30 (compliant guide) and 36 (pin supports).

| # | task | checks |
|---|---|---|
| 01 | tube swept along a 3D Bezier, flanged ends | bore and bolt holes clear (proven, res 0.1) |
| 02 | involute spur gear pair, 20 / 30 teeth | disjoint over a tooth pitch (proven); wrong-phase control disproven |
| 04 | 12-vertex non-convex polygon plate | volume within 3.5% of the shoelace value (an estimate) |
| 06 | mass, centre of mass and inertia of a plate | within 0.1% of the analytic polygon moments |
| 07 | bolted flange, ISO M10 threads, placed by frames | hole and thread clearances proven; no overlaps |
| 09 | 6204 bearing housing with press fits | fit bands proven; only the declared fits overlap; mass within 0.5% |
| 10 | slider-crank over a full turn | rod angle within 16.7 deg, slider in its guide, rails clear (proven over t); no overlaps |
| 11 | two-stage spur gearbox, ratio 6 | both meshes proven at 21 angles; controls disproven; no overlaps over t |
| 18 | screen with 10 000 holes; 10 001-instance assembly | both forms' volumes within 0.5%; watertight; no collisions |
| 19 | spline sweep, varying sweep and lofts | volumes vs area x length and frustum; loft sections; watertight |
| 20 | organic mouse shell, redistanced | watertight; shell wall within 0.1 of 2 mm, naive shell off |
| 22 | imported STL: offset, shell, redistance | volume, ray distances, bore radii, shell topology |
| 26 | drafted, ribbed enclosure with a lid | lip clearance proven, looser control disproven; wall thickness; screw holes by probes |
| 31 | sheet-metal channel bracket | ligament rules proven, controls disproven; volume vs formed and blank |
| 32 | four-bar crank-rocker over a turn | angle range, transmission, clearances, pins (proven over t); solve vs closed form |
| 33 | 58k-triangle imported knob | error bound; volume; claims measured on its samples; remesh; rays |
| 34 | ISO 286 fits over their tolerance bands | fits proven; controls disproven by sizes their tolerances allow |
| 35 | tolerance stack, worst case | worst case proven; tight design disproven in band, proven at the values |
| 37 | slider placed by a solve, over tolerances | clearance and solved ranges proven; tight stop disproven in band |
