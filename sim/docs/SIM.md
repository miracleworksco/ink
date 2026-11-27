## Studies (simulation)

A study asks how a part answers its loads: how far it moves, how highly
it is stressed and how near yield, at what load it buckles, how it
vibrates, how hot it runs. ink solves it on the shape's field itself, with
no mesh, and decides claims on it (`max-displacement`,
`min-safety-factor`, ...) as it decides any other.

### Quick start

A 10 mm square aluminium bar on two supports, pressed at its middle,
with its hand calculation in its comments:

<!-- example examples/beam.ink -->
```lisp
; A 10 mm square aluminium bar on two supports 200 mm apart, pressed down
; with 1000 N (about 100 kg) at the middle. Is it stiff enough, and strong
; enough? By hand (E 68900 MPa, yield 240 MPa, I = 10^4 / 12 = 833 mm^4):
;   deflection  P L^3 / 48 E I = 2.90 mm, and 0.02 of shear: 2.93 mm
;   stress      (P L / 4) c / I = 300 MPa: a safety factor of 240 / 300 = 0.80
; ink require examples/beam.ink: stiff-enough holds (2.94 mm), and
; strong-enough fails (0.79): the bar would bend for good.
; expect-status: strong-enough fails estimate
(import "std")

(body bar (shape (move (box 210 10 10) 100 0 0)))   ; x from -5 to 205, y and z from -5 to 5

; A support: the bar's bottom edge at x0, held within 1 mm of it, free to turn about it.
(fn edge (x0) (- (len2 (- x x0) (+ z 5)) 1))

(study press
  (shape (shape-of bar))
  (elastic (material aluminum-6061))
  (fix pin (edge 0) (pivot (at 0 0 -5) (axis 0 1 0)))
  (fix roller (edge 200) (dofs z) (pivot (at 200 0 -5) (axis 0 1 0)))   ; held up only: free to slide
  (load hand (force 0 0 -1000N) (region (- (abs (- x 100)) 1)))          ; spread over the middle 2 mm
  (res 0.5))                                                             ; cells of 0.5 mm

(require stiff-enough (max-displacement press < 4))       ; mm
(require strong-enough (min-safety-factor press > 1.5))
```

Run from sim/ (seconds on a laptop's CPU):

<!-- run -->
```sh
ink require examples/beam.ink
```

```json
{"counts": {"yes": 1, "no": 1},
 "requirements": [{"name": "stiff-enough", "verdict": "yes", "value": 2.94, "res": 0.5, "grid_error": 0.004},
             {"name": "strong-enough", "verdict": "no", "value": 0.79, "res": 0.5, "grid_error": 0.069}]}
```

- `holds`: the bar bends 2.94 mm at most (by hand, 2.93), under the 4
  allowed. ink solved it on cells of 0.5 mm, and again on cells of 1 and 2
  mm; from how far the answer moved between them its grid error is about
  0.004 mm, far less than the 1.06 mm margin, so the claim is decided.
- `fails`: its least safety factor is 0.79 (by hand, 0.80: yield 240 MPa
  over a peak stress of 303), under the 1.5 asked, and checked the same
  way: its grid error, about 0.069, and what a stress's read may miss
  (Study claims) are far less than its 0.71 margin. `at` is where: the
  bottom face at the middle.
- `unknown`, which neither is here: the grids disagree by more than the
  margin, or the cells cannot see a feature that decides it (a wall
  thinner than two cells, a gap they glue shut). `undecided` says what,
  and `hint` what would decide it, often a finer `(res r)`.
- `accuracy: null`: every number is an estimate from the cells, not a
  bound (Study claims). `ink solve examples/beam.ink` gives the whole
  solve: displacement, stress and where it peaks, the safety factor,
  reactions (500 N at each support), the buckling factor and the solver's
  statistics.

### Forms

```lisp
(study name
  (shape field)                                   ; the solid to analyse (fuse parts with min)
  (elastic (E 2100MPa) (nu 0.38) (density 1270kg/m3) (yield 45MPa))   ; MPa, -, a density, MPa; each may be a field; nu defaults to 0.3
  (fix wall region [(dofs x y z)])                ; held at zero; region: a field (held where < 0) or (face name ...)
  (fix seat r (dofs normal))                      ; frictionless: the part's surface in r held along its own normal (a bore: radially)
  (fix press r (dofs z) (move 0 0 -1))            ; held moved by the displacement given (mm; fields of x y z)
  (fix pin r (pivot (at 0 0 0) (axis 0 1 0)))     ; its mean motion held, but free to turn about the axis (no axis: a ball joint)
  (fix guide r (slide 1 0 0))                     ; its mean motion held, but free to slide along the direction
  (fix lug r (pivot (at 0 0 0) (axis 0 1 0)) (rigid))   ; held rigidly, but free to turn as one about the axis
  (load books (force 0 0 -49N) (region r))         ; total force (N) spread over the material in r
  (load inside (pressure 1.5MPa) (region r))         ; surface pressure (MPa, any field) on the surface in r
  (load drive (torque 0 0 5000N*mm) (region r))       ; torque (N mm) about the centroid of the material in r
  ; a force's or a torque's components may read what a motion's model works out, (nth (reaction PATH) 0): taken at
  ; the motion's worst moment when the study is asked (DYNAMICS.md, A study's load from the motion)
  (gravity 0 0 -9.81m/s2)                             ; m/s^2, needs density
  (res 1)                                         ; grid cell size in mm (default: size / 100)
  (support-margin 2))                             ; mm; stress peaks this close to rigid fixes are reported apart

(study hot
  (shape field)
  (heat (k 0.2W/m/K))                                  ; conductivity W/(m K), may be a field
  (temperature base region 20degC)                    ; held temperature: may vary over the region (x y z), and in a transient in t
  (convection air (face heatsink) (coefficient 25W/m2/K) (ambient 25degC))   ; W/(m^2 K), C
  (source chip (power 5W) (region r)))             ; total heat input (W) over the material in r
```

`ink solve FILE --study NAME` runs a study by the finite cell method: a
regular grid covers the part, and cells the surface cuts are integrated
by sampling the field, so lattices and material fields need no mesh. The
grid's planes lie at multiples of `res` from the origin, so a face or a
region's edge on one is exact. Two cells share a node only where their
material touches, so the sides of a gap at least half a cell wide move
apart. Cells cannot hold a feature thinner than two of them across, of
material or of air: a cell's field across thin material holds it as if
it filled the cell (a fin, a web or a ligament too stiff, its stress and
motion unread), and across thin air moves its two sides as one (a slit
or a clearance shut, a slit's hinge held by one cell). Each cut cell's
samples measure the feature they lie in, its chord through their phase
(the run of material, or of air, along the field's gradient), and a cell
holding one thinner than two of its cells is split into cells half as
big, sampled again, and those still holding one again, until each
feature is two cells across, to an eighth of a cell: the solve says how many
(`cells_refined`) and how fine (`refined_res`). A cell under 3% material (a
sliver) is left out once splitting is done (`dropped_sliver_volume_mm3`),
and one beside a split cell its size is split as it is: a sliver is one
at the size of the cells round it.

A piece disconnected from every fix or convection surface is included in
transient heat: its initial temperature and heat capacity determine its
response. A modal study with no fix keeps it too: the part is free, its
rigid motions no modes (below). Other studies drop such pieces, a modal
study with fixes among them. The solve reports the dropped
volume (`dropped_floating_volume_mm3`). Study claims and claims that
read its sampled fields answer `unknown` and name the missing volume.

A feature still thinner than two of the finest cells (or where the
split cells would outnumber the grid's) is reported ("a feature thinner than 2 of its cells is
left"), and no claim on the solve is decided; refine (`res`). Except within
two of the finest cells of where its phase is resolved: a taper's edge
(a wedge, a blade's trailing edge, the air at a notch's root) is thinner
than any cells next to it, and its error there is the finest cells'
own, which a claim's coarser grids measure (Study claims). A cell takes
one value of each material property (E or k, nu, the yield claims read,
the density where the solve reads one: modes, heat over time, gravity;
and 1/E or 1/k, what cells stacked in series sum, so a dip toward 0 or a
film far below its neighbours is judged as the cells average it),
so a material feature thinner than the cells (an insulating layer, a
stiffer strip) would be lost alike on every grid a claim checks. A cell
whose properties are not one value across it (by the fields' bounds over
its box) is split as a thin feature is, unless its variation is one the
grids measure (its range over boxes a quarter, a half and the whole of
the cell growing as a smooth property's does, at least linearly and at
most quadratically, read as a claim's three grids read a value) or it is
one interface between materials each two cells thick beside it, each
reaching the step (along an axis, to within a cell of the finest
cells; at the finest, at any angle). A cell holding two materials takes
their laminate's element: each one's share of the cell and the face's
normal, measured along lines through its samples; its stiffness the two
in series across the face and side by side along it; each material's
stress its own strain (the cell's, with the jump the face makes) times
its own stiffness. A flat face along the cells is so held where it lies
on every grid (a soft band 0.2 mm off planes 2 mm apart reads its stress
to 1e-5, where one stiffness per cell read it 5.3% low on all three). A
face oblique to the cells, or curved, is not: where it crosses a row of
nodes, the node takes the soft side's stretch into the stiff cells
beside it, an error that does not shrink with the cells (a soft band in
a pulled bar turned 1 degree read 3.5 MPa against 2.45), so a claim on
the solve is unknown, saying so; a face at a small angle, its cells
still in layers, is not refused. A cell holding a third material takes
its centre's values. A film can still pass beside other
variation inside its own finest cell; a stiff one so hidden reads its
stress and safety factor low by about its E over the E beside it, since
stress is read point by point, not averaged.
One left so at the finest cells is reported ("a material thinner than 2
of its cells is left") and no claim on it is decided, holding or failing
(the layer lost, or filling a cell its centre lies in). A material's
face on a cell's face splits nothing, and a layer whose value lies
between its neighbours' is not seen (the cell's range is theirs). A gap the
cells' samples cannot part at all is glued shut. Where the part does
not join the gap's sides round within the cells about it
(a slit, a saw cut, a clearance; not a corner, a groove under a cell deep
or a notch's root, whose error shrinks with the cells), the solve says
so ("glued shut": the gap's width, where, and the cells that part it): it
solved the part with the gap filled, and decides no claim on it,
holding or failing (Study claims). A gap's sides join only where bounds show the material solid
across it; pieces that meet flush (a lattice set on its skins, a union
of boxes) are, as everywhere in ink. Where neither that nor air is shown
(a hairline under any sample), the solve decides no claim either. A notch sharper than about 15 degrees
is glued at every cell size, so claims on it stay unknown. Some parts'
solves stop converging on fine cells (a bracket with bores held under
washers, below 1 mm cells): their claims are then unknown. Results:

- maximum displacement, or maximum and minimum temperature;
- peak von Mises stress, or peak heat flux, with where (`at`), and
  stress percentiles (of the cells' centres). The peak is read as a
  probe reads it, at the part's surface (a mesh of it at half the cells
  it lies in, a split cell's finest child's)
  and at the cells' centres, so stress in bending and at a concentration
  is read where it peaks. A mesh's facets cut a convex edge's corner off,
  and stress in bending falls into the part from the edge: a square bar
  turned 37 degrees read its peak 0.14 of a cell inside its edge, 4% low
  on every grid. That read stands; what the solve's own field reaches
  past it is measured by a climb from each cell's highest surface read,
  along the surface (its zero set) while the cell's stress rises and the
  point is the cell's or no cell's, within half a cell of it, and a
  claim's allowance over the read widens by it (Study claims). It (`max_von_mises_mpa`, what `max-stress` and
  the safety factor read) leaves out the support margin of each fix
  that holds rigidly (`excludes` names them; a symmetry plane has none),
  whose own peak is `max_von_mises_near_supports_mpa` (Supports);
- the material as solved (`material`): E, nu, yield. An elastic study
  gives its Poisson's ratio, or a material that has it: none is assumed. Ink checks E in 1e-4..1.2e6 MPa, yield in
  0.001..100000 MPa, and Poisson's ratio strictly between -1 and 0.5.
  Thermal conductivity is checked in 1e-5..1e5 W/(m K). These limits
  apply to the value in the material. Use `(E 1e6Pa)` for a 1 MPa
  silicone. Bare `1e6` means one million MPa and also fits the range;
- minimum safety factor (with `yield`): the least of the yield over the
  stress at every point the stress is read, each with its own yield, so
  each material's own sets it where it is, and where (`at`). A cell not
  of one material (its E, nu or yield not one value over where it is
  read: its box and its reads, up to half a cell past it where a face in
  a cell under half material is read from it; a strip, a gradient, an
  interface between bodies) gives no value: no point read in it is
  surely its material's (a soft, weak skin read the bulk's stress
  against its own yield). It is bounded instead: its least yield over the
  stress's bound over that region (its element's field at the region's
  corners, at the stiffest E / (1 + nu) there: a material stiffer than
  the cell's mean, side by side with it along the load, carries that
  much more; a cell of two materials, each one's own stress over its own
  yield; a least yield not above 0, as a blend's interval can reach,
  bounds the factor at 0). Away from the
  supports the bound is left out as `yield-between-reads` (`ink solve`:
  `yield_between_reads`), checked on the value's grids with its grid
  error; in a support margin it counts in the margin's value, as
  singular as the margin's stress. A bound is the least the factor can
  be, so a claim it holds holds, with no value at all (a part whose
  yield varies everywhere answers on its bound, `band`); a claim fails
  only on a value. Material that varies inside cells can so hold a claim
  but never fail one: a false claim on it stays unknown. Cells under half material, the outer
  part of any surface's cells as well as walls thinner than about a
  cell, are neither read nor bounded (the climb above may reach their
  surface from a cell beside them, widening a claim's allowance, not
  its read);
- applied loads and support reactions, each `total` a force in N (a
  vector, with `moment_nmm` about the origin in N mm) or, in a heat
  study, a power in W (a number), and an equilibrium check
  (`equilibrium_imbalance`, covering both forces and `moment_imbalance`
  in elastic studies). The elastic check takes moments about the center
  of the part's bounds and scales them by the part's size and the loads
  and reactions;
- the largest rotation (static elastic), solver statistics, warnings;
- `accuracy: null`: every number is an estimate, as a claim on a study
  is (a claim compares three grids: Study claims).

Flags: `--points` probes displacement and stress, or temperature and flux, at
points (a table: Geometric queries), each with its `estimate`: one read on
one grid, its grid error not measured, not a bound; in a cell part
material, a stress that scatters about the true one by a share finer cells
do not shrink (a ring's bore read 26% low). A point in no cell of the
solve (in air, a bore, or material dropped as floating) has its row null
and its reason under `outside`, `[row, why]`. `--render out.png` colours the part
by stress or temperature, deformed, in render's views, with the colour
scale under them (the quantity, its unit and range, and how much the
displacement is exaggerated).

- **Elastic studies are linear** (small displacements and rotations).
  Above 5 degrees of rotation (`max_rotation_deg`) the study and its
  claims warn ("large rotation"); ink has no large-rotation solver.
- **What a static elastic study does not model**, so its stresses and
  safety factors say nothing of: buckling; yield and plasticity (past
  `yield` the answer is still elastic, read the safety factor); bolt
  preload; fatigue; contact, friction and separation (bodies are bonded,
  supports pull as well as push); thermal stress (heat and elastic
  studies are not coupled); large deflection; body loads other than
  gravity (a spinning part's, an accelerating frame's: a load's force
  cannot vary with x, y, z, and is refused if it does). Every static
  elastic solve that converges finds its linear buckling factor (below).
  Where its band reaches 1, the loads may buckle it: the study warns
  "buckles", and its stress, safety-factor, displacement and compliance
  claims that hold are `unknown` (Study claims). A safety factor counts
  a buckle as its band's low end, so a claim on it is `unknown` where
  that is under the claim's bound (a column that buckles at 2.8 times
  its loads holds `> 2`, not `> 5`). Where the eigensolve resolves no
  factor, it warns "buckling not resolved", and the same claims are
  `unknown`. A declared mirror includes buckles across its cut. An explicit
  `(excluding buckling)` can ask for the linear response while leaving
  out a buckle in the reflected part. The modeled portion's own factor
  still counts on a reduced model, even under this exclusion.
- **Linear buckling**: `(min-buckling-factor study > k)` on a static
  elastic study claims the least factor on its loads that buckles it,
  the lowest lambda of K phi = -lambda K_G phi (K_G the geometric
  stiffness of the solved stress), found by LOBPCG on the solver's
  multigrid at one to three static solves' work. `(study-field s
  displacement (mode 1))` reads its buckling mode, 1 mm at its peak. On
  4 cells across, Euler's columns read within 1.5%; a plate pressed in
  its plane on 2 cells through it, 2.5% low; a tube of radius 20 times
  its wall, 6% under thin shell theory, which leaves out the wall's
  shear. It is linear: imperfections lower a real part's buckling load,
  a thin shell's to half or less, so claim a margin; a load that follows
  the surface (a fluid's pressure) is taken as fixed in direction. With
  `(mirror name x 0)`, the other half is the reflection of the material,
  physical supports and loads. The static solution holds normal displacement
  at the cut. Its stress is kept for every buckling search. A symmetric
  mode holds normal displacement there, and an antisymmetric mode holds
  the two tangential components. The elastic and prestress energies of
  the reflected problem are invariant under the same reflections, so
  the modes split into these independent sets, including repeated modes.
  The solver checks every combination on up to three perpendicular cuts,
  and takes the lowest ends of the factor bands. The returned mode is the
  one with the lowest computed factor. As elsewhere these factors are
  finite-cell estimates. `sectors_solved` counts the searches and
  `symmetric_factor` gives the original portion's factor when finite.
  `antisymmetric_axes` lists the odd reflections for the returned mode.
  A free rigid rotation with no prestress force on any degree of freedom
  leaves the factor unchanged. For a turn that couples to the stress,
  the solver minimizes the geometric energy over its amplitude, or reports
  an unresolved factor if positive restoring energy cannot be established.
  The factor is a band (`band`): its mode's
  Rayleigh quotient, never under the lowest factor, down to what the
  eigensolve's residual allows (it stops at 1%). That low end takes the
  mode found as the lowest and reads its residual through a multigrid
  stand-in for the stiffness's inverse, so it is checked by what assumes
  neither: the stiffness plus the loads' geometric stiffness at that
  factor must be positive definite, which conjugate gradients on a
  random right-hand side either meet a direction against (a proof that
  the loads at that factor buckle the part, whatever the preconditioner)
  or converge without, to a hundred-thousandth of it (every mode has a
  part in a random vector far over that). It is checked at the low end
  found, or at 100 times the loads where that is more: the low end is
  then 100, what was checked. Where it fails the low end is the largest
  of its halves that passes. This is a numerical check to a
  tolerance, not a count of the matrix's negative pivots; a claim
  decides on the band widened by its grid error, an estimate as that is. Where the eigensolve cannot resolve a
  factor to 1%, the factor is unknown, saying so; where it bounds the
  least factor above 1000 it stops there, and a factor that moves that
  far between grids is left unknown by its grid error. Where its search
  finds no direction the loads compress (a bar pulled), it found no
  mode, which is not none: its factor is "at least" the largest its
  search excluded (above 1000, about 20 steps, three times a static
  solve), with no top. A pulled bar's clamp compresses its corners,
  which buckle at about 16,000 times its load, past where it stopped:
  `> 3` holds on it, `> 1e6` is `unknown`, `< 10` fails. That bound,
  like every factor's low end, takes the mode the search found as the
  lowest, which is not proven. Where no cell's stress is compressive at
  all, nothing can buckle: the factor is none, with no search
  (`iterations` 0). Where the only compression is under a roundoff
  allowance, there is no search either, and `low` is finite: the least
  stress that buckles a member the grid holds, over that compression
  (a rule of thumb for struts, plates and torsion, a thousand or more).
  One stress per cell is read, as the search reads it: compression that
  varies inside a cell is seen by neither. A clamp
  thicker than a cell compresses its outer cells slightly, and then the
  search decides.
- **`(modes n)`** on an elastic study with a density finds its lowest n
  natural frequencies instead of a static response: each mode's frequency
  in Hz and its effective mass fractions (along and about x, y, z). Loads
  are ignored. The rigid motions a part's supports leave free, or all six
  of each piece where it has no fix (a free plate, a ring on strings),
  are each a mode at 0 Hz, no vibration: they are found as a static
  study finds them, each piece's apart, kept out of every vector the
  eigensolve takes (mass-orthogonal to them), and the study says which
  and how many (`rigid_modes`). A free steel bar 200 x 10 x 10 reads
  1285.5, 3487.3 and 6686.1 Hz bending on cells of 1 mm, Timoshenko's
  1285, 3486 and 6676; on two rollers, 575.6 Hz for 575.6. A study with
  no fix is of a free part, and `min-frequency` is its least elastic
  mode. A study with fixes that leaves a motion free (a roller
  forgotten, or meant, as a beam's slides on two rollers) has a mode at
  0 Hz: `min-frequency` counts it (`rigid-motions`, `value_if_counted`
  0), so `>` is unknown until the motion is held or the claim says
  `(excluding rigid-motions)`. A piece its fixes do not reach is
  dropped, and its claims are not settled, as a static study's: a part
  cut from its clamp by a gap is not read at its free modes. A
  symmetry plane keeps only the modes symmetric about it, and the lowest
  is often one it drops (a strip modelled as half its thickness loses its
  bending: 819 Hz for 166), so a `min-frequency` claim on a study with
  one that holds is `unknown`: solve the whole part to include its antisymmetric vibration modes.
- **Transient heat**, `(transient (duration 600s) (steps 60) (initial
  20degC))`, integrates by backward Euler from t = 0, sources' magnitudes
  reading t (the study's own time, not the model's: a claim on it asks
  for no t-span); it needs `density` and `capacity` (J/(kg K)). Each step
  takes all the heat its sources give in it: their power's integral over
  the step, bounded piece by piece over its times, so a pulse shorter than
  a step is not missed, and every run's steps end at each moment a
  source's power may jump (a pulse's end is read where its heat stops,
  not a step later; a source that may jump more often than its `steps`
  is refused). A held temperature may read t: each step holds it as it
  is at the step's end (its holds' terms made again there, about half
  again a constant hold's solve time), and the steps end where it may
  jump too (over the cells its fix holds), the one ending there, or at
  the duration, holding its value before the jump. A held temperature
  reading t in a study that is not a transient is refused. The same run on steps twice and four times as long
  checks its steps as coarser cells check its grid: `max-temperature` is
  the peak over time, where and when it is, with a band of its error in
  time (an estimate, not a bound), and none where the three counts do not
  converge (each difference not smaller than the last, or of the other
  sign, by more than the solves resolve): steps that all spread a short
  pulse farther than it reached read its peak low alike, growing as they
  shrink, and no error is read from them. A claim
  its band crosses, or one with no error read, is `unknown`, saying to
  take more `(steps n)`. Results add a `transient`
  block with the peak temperature's history, its place and time, and that
  `time_error`; the fields are the final state. Elastic studies are static
  or modal.
  A transient study keeps pieces with heat capacity but no temperature
  fix or convection. The solve reports each one's volume and location.
  In every physics, when the mesh keeps multiple pieces, ink checks the
  original shape for a gap between each pair. It traces cell faces around
  a piece and requires strictly empty space across a closed surface.
  This also works for a pin inside a hole. If the shape could connect
  the pieces through a missed web, whole-study claims are unknown.
  Refine the cells across the possible connection.
- **Pressure** is a traction of -p n over the surface in the region, p
  any field (a hydrostatic head `(* 0.00981 (- 100 z))`), integrated over
  a mesh of the surface at half the cells it lies in (where cells split
  round a thin feature, half their finest child), each facet's part
  inside the region (so a face counts whole under a region however thin,
  `(- x 0.01)`); convection and a frictionless hold are integrated the
  same way, convection's heat removed reported as its reaction (W). A
  pressure or convection whose region meets none of the surface is
  refused. Convection alone can hold a heat study. When the shape lies
  exactly on the grid faces and the shape intervals establish this,
  the study uses those faces for pressure on each comparison grid.
  Where it can establish a convex part's bounding planes, it clips each
  plane against the others so the pressure covers the face out to its
  edges, including along an oblique taper or at a symmetry cut. A load
  region bounded by planes clips those faces again. Uniform pressure
  then integrates exactly over the selected area.
- Each such integral is checked against a second count of its region's
  surface: lines through the cells, as far apart as the mesh's planes
  (half a cell, half a split cell's finest child) and sampled halfway
  between them, counting where they cross the surface in the region. It is
  a second sampling, not a proof: a sheet between two of the mesh's planes
  thinner than half their spacing can pass both (a feature that thin is
  split, or its claims are unknown, before this). Where a region's rim may
  cut a line's cross-section, the line is counted again on 16 sub-lines
  across the rim, so a strip a line or two wide (a line load, a knife
  edge) is measured to half a sub-line at each rim; one narrower than a
  sub-line (1/32 of a cell) is not measured, and never flagged. The two
  are compared in patches of 4 grid cells to a side and over the region,
  each against the lines' noise and rim rounding, so a feature's surface
  lost in a large region still shows. `ink solve` reports
  `surface_measures`: `covered_mm2`, `lines_mm2`, `uncertainty_mm2` and
  `shortfall_mm2` (sum |n_a| dA: a face's area where it faces an axis).
  For pressure, `excess_mm2` reports where the applied area exceeds the
  independent count beyond its uncertainty. Refine across the region's
  rim when either discrepancy prevents a pressure claim.
  An integral short of the lines' by more than their own uncertainty, in a
  patch or over the region, loaded, cooled, heated or held less of the
  part than written, by how much no claim can say (a load missed on a
  lever moves a tip far more than its share of the region): it is warned
  of, and the study's claims are unknown, failing as well as holding.
  For pressure the solver also bounds the variation of the specified
  pressure field over each covered facet, including facets whose sample
  reads zero. Each facet is checked against its own sampled force.
  A pressure that varies linearly over a whole facet can be integrated
  directly, including its moment. For other fields with at most 3%
  possible variation, the solver can subdivide that facet and bound the
  remaining intensity error to the solve tolerance (at most 1e-8
  relative). If the variation remains unresolved, or force misses the
  cells, static claims reading that load are unknown. The message names
  the load and location to refine. A larger force on another face cannot
  reduce the discrepancy. `ink solve` reports the largest remaining
  fraction as `pressure_error_fraction` (unbounded for a potentially
  loaded facet whose sample reads zero). The separate
  surface-area check still applies to the selection itself.

```lisp
(import "std")
(body bar (shape (move (box 100 10 10) 50 0 0)))
(study heat-pulse (shape (shape-of bar)) (heat (k 0.2W/m/K) (density 1270kg/m3) (capacity 1300J/kg/K))
  (temperature wall (- x 0.5) 20degC)
  (source torch (power (* 1W (select (- t 60s) 1 0))) (region (- 95 x)))
  (transient (duration 600s) (steps 60) (initial 20degC))     ; initial temperature, C
  (res 2.5))
```

Read results with care: they are numerical estimates, not proofs. Smooth
parts converge at about second order in the cell size; sharp inside
corners and crack-like edges converge at first order or slower, and
their peak stresses grow with refinement (the idealised geometry): use
fillets, and read the percentiles. Supports hold material rigidly, so
the stress where a held region ends on the free surface is singular: it
grows as `res` shrinks. Stresses within `support-margin` of such a fix
(2 mm, and at least 2 cells: `within_mm`) are reported apart, as
`max_von_mises_near_supports_mpa`, which bounds nothing;
`max_von_mises_mpa`, the percentiles, `max-stress` and
`min-safety-factor` are the rest of the part. A symmetry plane has no
margin. Its displacement condition acts on the plane itself, so the stress
there contributes to the peak. A uniform normal restraint covering an entire
face also has no edge within the face. Ink counts its stress when the
held layer is at most a tenth of a cell thick. A hold that ends inside a face it holds normally (a pad) is a punch's edge, and
keeps its margin; so does a hold of two or more axes. A clamped plate's
edge stress, its design stress in plate theory, lies in its clamp's margin:
read it from a model of the clamp (the ring that grips it). The peak says which
fix's margin is nearest and how far (`nearest_excluded`: the fix, and the
distance to the centre of the nearest cell left out); where that cell is a
neighbour of the peak's, the peak is at the margin's edge, where the fix's
singular stress may reach past it: a `note` says so, and `max-stress` and
`min-safety-factor` claims warn. Check it on finer cells (a singular peak
grows), widen `support-margin`, or model the contact. A claim on them
also says what it leaves out (`excludes`: the margin's peak, as the
metric reads it), and is `unknown` where counting that would change its
answer (Study claims).
The solve warns when most of a load lies within the margin (a load beside
a clamp: its stress is all in the excluded part), and when most of it
lies inside a region a fix holds (it goes straight into the reaction).
Where the support's own stress matters, model its contact (a pin, a pad)
rather than a rigid fix.

**Materials.** `(material name)` gives typical properties as data (`E`,
`nu`, `density`, `yield`, `k`, `capacity`) from ../kernel/std/materials.ink (`ink
std --section materials` lists them). Each has one name: another spelling
(`PETG`, `aluminium`) is refused at load with the name to write. A
material can also be your own data with the same keys. A body's `(attr
material ...)` is one of these names (or such data); free text goes in
`(attr material-note "...")`.

- `(elastic (material petg))` and `(heat (material steel))` fill in every
  property; explicit keys override them.
- `(by-body target key [default])` is a field holding each body's
  attribute `key` (the nearest body's), else `default`: for an assembly
  of materials, `(elastic (E (by-body root E)) (nu (by-body root nu))
  (yield (by-body root yield)))`: each body's own yield sets the safety
  factor where it is. Bodies are bonded where they touch, and a bonded
  interface's edge on the surface is singular (its stress grows as the
  cells shrink, as a sharp corner's).
  A body's number is bare (its attr's, else its material's, in ink's
  unit), so a key whose quantity is not the bare default's (heat's `k`)
  gives the whole field its unit, and a default bare: `(heat (k (* 1W/m/K
  (by-body root k 0.2))))`. Below, a heater plate on a printed block, the
  block's temperature claimed alone (Results as fields):

<!-- example examples/two-materials-heat.ink -->
```lisp
; A 2 W aluminium heater plate on a printed PETG block: each body conducts
; by its own material, and the block's temperature is claimed alone.
(import "std")
(body block (shape (move (box 40 40 8) 0 0 4)) (attr material "petg"))
(body plate (shape (move (box 32 32 8) 0 0 12)) (attr material "aluminum-6061"))
(assembly root (instance block block) (instance plate plate))

(study warm (shape (shape-of root))
  (heat (k (* 1W/m/K (by-body root k 0.2))))    ; each body's k, else 0.2, given its unit whole
  (source heater (power 2W) (region (shape-of root/plate)))
  (convection air (- 0.5 z) (coefficient 10W/m2/K) (ambient 25degC))   ; all but the bottom
  (res 1))

(require hottest-under-120 (max-temperature warm < 120degC))
; One body's temperature: the block lies where the field is under 100 (C).
(require block-under-100 (subset (shape-of root/block) (- (study-field warm temperature) 100)))
```

<!-- run -->
```sh
ink require examples/two-materials-heat.ink
```
```json
{"ok": true, "all_yes": true, "counts": {"yes": 2}}
```


### Study claims

The metrics a study claim reads:

```
(max-stress study < v)  (max-displacement study < v)  (min-safety-factor study > v)  (max-temperature study < v)
(compliance study < v)  (min-frequency study > v)  (min-buckling-factor study > v)
```

A claim never answers past what its metric leaves out. Each thing left
out is listed under `excludes` (`key`; `value_if_counted`, the metric's
value were it counted: the claim's own `value` where counting it moves
nothing, `null` for an unmeasured possible effect or an unbounded value; `what`; `at`): `support-margin`, the
stress within a rigid fix's margin, and `buckling`, a buckle (a linear
buckling factor whose band reaches 1, or none resolved; on stress, safety
factor, displacement and compliance); `yield-between-reads`, the bound of the cells not of one material
(on safety factor); `rigid-motions`, a supported modal study's free
rigid motions (on frequency, 0 Hz); `symmetry-modes`, unsearched
mirror modes (vibration, or an ambiguous legacy hold). Neither those modes
nor `modeled-buckling`, the reduced model's own stability, can be excluded. Where
counting one (with its grid error, from its value on the claim's grids, as the value's is found; the support margin's, singular, as read) would change the answer, the claim is `unknown`:
`undecided` says what was left out and where, `hint` how to decide it
(model the support's contact; check the buckling load). A claim may
leave one out on purpose, its answer then on the rest, with a warning
saying so:

```
(max-stress study < v (excluding support-margin))   (min-safety-factor study > v (excluding buckling))
```

A study claim is decided, holding or failing, only where its grid saw
the part. It is read on three grids, each twice as coarse as the last:
its own cells and cells twice and four times as coarse, or, where those cannot read
the part (no cell of them at least half material, a load's region
empty, a gap glued), cells half and a quarter as fine, so that the
coarsest still reads it; it is answered on the finest (`res`). Its grid
error is estimated as Roache's grid convergence index: 1.25 |x - x2| /
(2^p - 1), never under 1.25 |x - x2|, where the three converge at an observed order p up to the
metric's (2; a stress's and a safety factor's 1), and three times their
whole movement, 3 |x - x4|, where they converge faster than that (a
feature the coarsest grid cannot see). Where they oscillate, or are close
but not in order (moving apart, or a swing growing by more than three), the
error is three times the larger difference, unless the solve has found
something that could make the value grow without bound as the cells shrink:
a load or hold less than a cell across (100 N on a 0.1 mm patch read
1.022e-4, 3.262e-4 and 7.736e-4 mm on cells of 4, 2 and 1 mm, and 6.336e-3
on 0.125); a peak stress at a hold's edge, a sharp corner, where materials
meet, or on a curve finer than its cells; a buckling or vibration mode that
shrinks with the cells. Then none is read, and the claim says why; a reading
the solve does not classify is refused the same way. Never under |x - x2|.
A difference no larger than what the solves themselves resolve is no
evidence of a trend, asked before any other: three grids with one such
are judged as two, by three times the other. A claim that holds or fails by less than its
grid error, or that has none, is `unknown`, `undecided`
giving the three values, `hint` a finer `(res r)`; so is one no such
three grids read, saying why. A peak stress's read can also miss in a
way all three grids share, which their differences cannot show. Within
1.5 of the claim's own leaf cells at the peak (a split cell's own size,
not `res`), the surface is surveyed for its tightest concave radius (a
hole, a fillet, a notch), by the field's sign alone, so any field with
that zero set reads alike. A radius under 1.45 of the cells the peak is
read in (the references at 1.5 read within the bound below) is read
only in part: a bore of radius 1 on 1 mm cells read 0.59 of its peak,
its grids within 2%. So is a crease: where the surface's normal turns,
across a span of 0.28125 of the claim's cells, by more than 0.001 rad
and by more than three quarters of its turn across 0.5625 (a curve's
turn halves with the span): a sharp inside corner, or a round under
about 0.4 of a cell. A turn of 0.001 rad or less across the shorter
span is read as a curve. There a `max-stress` or `min-safety-factor` claim is
`unknown`, holding or failing, saying the cells that would resolve it
(at a crease, that its stress is singular, or the round too small for
them). So is one whose coarser grids cannot resolve the peak's curve: no
finer grids are taken for it, as they are for a thin wall. Elsewhere a
claim is decided only past its grid error and what the read may miss:
over it, 22% over the cells the peak's radius spans (at least 4%: 15%
at 1.45 cells, 4.4% at 5), and 10% under it; the share over is of the
most the solve's own field reaches climbing the surface from the read
(Results: a convex edge oblique to the cells), so the read stands and
only the allowance widens. A curve's cut cells read
its peak low by less as its radius spans more cells (a flexure's 0.75 mm
fillet read 19, 15, 6.3 and 1.9% low on 0.5, 0.25, 0.125 and 0.0625 mm
cells), and three grids share the miss. Fitted on 150 reads of 81 inputs
with known answers (closed forms, converged CalculiX): past their grid
error the worst lay 2.2% over at 5 cells per radius and 8.2% under (at
a bonded interface), each allowed about twice and 1.2 times, not an
envelope. Known to miss it: a bonded interface of high stiffness
contrast (a ring of E 200000 inside 2000 read its least safety factor
8.4 to 9.0 against 11.4); and a field that
is not like a distance near its surface (a bore written `(pow (- 1.5 r)
3)` read 2.67 against 2.95 at 12 cells per radius), since a cut cell's
material is read off the field's tangent planes and a thin feature
judged by its value over its gradient, which take the field for the
distance to its surface. `ink solve` gives the peak's surface
(`surface`: its radius, cells and whether they resolve it). Where the claim's grid splits cells (a
thin feature), each other grid splits its cells where that one would,
as deep, so their cells differ there as their sizes do; a coarser one
may leave features thin that the finer resolved, which only widens
|x - x2|. Grids whose finest cells are not each at least twice the
last's (split to the same cells) are no such three. Its `estimate`
gives its cells and grid error, `grid_error` the error as a number. This is a check,
not a proof: two grids blind to the same feature agree. It costs each
claim two more solves: an eighth and a sixty-fourth of its own, or, for
a part too thin for coarser cells, 8 and 64 times it.

Elastic studies give `max-stress`, `max-displacement` and `compliance`;
`min-safety-factor` needs a ductile `(yield ...)` (std's glass, cast iron
and carbide have named strengths instead); `min-frequency` needs `(modes
n)`, `min-buckling-factor` a static study (no `(modes n)`); heat studies give
`max-temperature`. A solve that broke down or did
not reach its tolerance has no results: its claims are `unknown` with the
reason, `ink solve` reports `solver.failure`, and it is not cached.
Before a solve allocates, it estimates its memory from its cells
(`grid.memory_estimate_mb`; rough: 0.43 to 1.35 times what solves take),
and where that is more than the memory and swap available (or
`INK_MEMORY`, `8G` or bytes) it says so on stderr and goes on: it may
run out of memory.

A study solves one design. Where its shape (or anything it reads) reads
a param with a tolerance, a plain `ink require` has no solve covering the
band, so its claims are `unknown`: `unavailable` names the params and
their bands, and `hint` the flags that check it at one size, `--param
bore-d=8.4:8.4` (an end of the band) or `--nominal` (every band dropped).
Which size is worst is not proven, so each answer is that size's. `ink
solve` works at the values.
Solved studies are cached as the kernel's answers are (Results). A study
of an unbounded shape is solved in the window (`ink check` warns).

### Named faces

A region may be a named face instead of a field:

```lisp
(fix wall (face upright (normal 0 -1 0)))                  ; that face held
(load books (force 0 0 -49N) (region (face board (normal 0 0 1) (depth 1))))
```

The face is the surface where the name decides the value, filtered by the
normal (within 30°) when one is given, so the reference survives parameter
changes. A load's region is the band of material under it, `depth` deep,
in mm (default 1.5 grid cells); a fix holds the face itself, its depth
only where its band finds the face.

### Results as fields

`(study-field study quantity)` runs a study while the file loads and
returns its result as a sampled field, usable anywhere a field is
(colour, another study's shape, a claim): `von-mises` (MPa),
`displacement` (mm), `displacement-x`, `-y`, `-z` from elastic studies,
`temperature` and `heat-flux` from heat studies. `(study-field s q (mode
k))` reads mode k (from 1) of a modal study, its shape 1 mm at its peak;
of a static elastic study, `(mode 1)` is its buckling mode.
The deformed part, to first order, is the shape warped back by the
displacement, `(warp part (- x (study-field s displacement-x)) (- y
(study-field s displacement-y)) (- z (study-field s displacement-z)))`,
so claims can check a loaded part (a snap-fit hook against its catch),
for displacements under about a cell. The study must come before its use
and not depend on the field it produces. `examples/field-driven.ink`
grades a lattice by stress. std's `color-map` draws a result as a
contour: `(attr color (color-map (study-field s von-mises) 0 120))`. A claim
on one body's temperature reads the field over that body's shape: `(subset
(shape-of root/block) (- (study-field warm temperature) 100))` holds where
the block is everywhere under 100 C (an estimate: the field is the
cells', its grid error not measured).

### Supports

A fix holds zero displacement on its region: on every axis, or on the
axes `(dofs ...)` lists. The region's material is held, and the rest is
free: their boundary carries the held value (Nitsche's method), its
traction the fix's force (its reaction), so a hold neither slips nor
stiffens the cells its edge cuts (a cantilever turned off the grid reads
CalculiX's deflection within 0.5% on 4 cells through it). A fix on a
named face holds the face. A region's edge is traced in each cell it
crosses, each of its flat faces there, and each piece of it reads its
traction partly in its own cell and partly in the free cell beside it,
the more the neighbour's the less free material its own cell holds
behind it: the reading does not jump as the edge crosses a cell's
plane. Free material thinner than 3% of a cell that no cell beside can
carry (a sliver between a hold and the part's face) is held with the
fix, and the study says how much. Where the edge falls in a cell still
moves the reading by the grid's own error, which the three grids
measure. A fix that holds one direction alone says what
it is: `(support)` for a physical restraint (a roller, a guide, a seat),
or the cut is declared, `(mirror name axis coordinate)`. A cut of
symmetry leaves out the buckling modes antisymmetric about it, so one
that says neither (and leaves no pivot or slide free, which only a
support has) leaves static claims not settled, saying so; a whole flat
face so held is solved both ways first. A study's shape is proven
bounded (one that is not is refused: its study would answer for the
part in the window alone). A claim on
`(study-field s von-mises)` is not settled: that field is each cell's
stress at its centre, for showing and shaping; `max-stress` and
`min-safety-factor` measure the peaks. `(dofs normal)` is a frictionless contact: the
part's surface in the region, held (as a face is) along that surface's own normal (a
bore radially, so a shaft in it is free to turn), so a region holds the
same way however it is written. It holds every surface it takes in: a
band round a bore that reaches the plate's faces holds their rims along
the bore's axis, and a layer on a face holds the rims of the faces beside
it; keep the region inside the edges of the surfaces it should leave
free. Where other fixes hold some of a facet's normal's axes, it holds
what they leave: u . n held along the rest of n, at what their values
leave of its own. Where they hold all of the normal at another value it
is refused as a contradiction, but at a corner: the facet lies on a face
another fix's region is a layer on (its gradient nearer the face's
normal than its plane), the contact's region crosses that face (its
gradient nearer the plane), and the contact has a face of its own, one
no other fix holds along its normal. There it leaves the facet to the
other fix whole, and says how much: a floor's strip at y < 0.01, or
2 mm deep, meeting a top a fix moves along z. None of it reads the
cells. A pad on a moved face, with no face of its own, or a lid, its
region a layer on that face, is refused.

- `(move dx dy dz)` holds the region moved instead (numbers or fields of
  x, y, z, mm); its reaction is the force the move takes. Compliance has
  no value in a study a fix moves. A field drives a submodel: a fine
  study moved over its cut by a coarse study's displacement fields (cut a
  few section depths from where the answer is read).
- `(pivot (at x y z) (axis dx dy dz))` leaves a turn about the axis free
  (a pin through a bore, a hinge, a beam resting on an edge: a pivot on
  the material within a small radius of the line); `(pivot (at x y z))`
  every turn about the point (a ball joint); `(slide dx dy dz)` a slide
  (a guide). With `(dofs ...)`, only the axes listed are held, relative
  to the free motion: `(fix roller r (dofs z) (pivot (at 200 0 -5) (axis 0
  1 0)))` is a roller on a line. A pivot or a slide holds its region's
  mean motion (as a distributing support, RBE3, does): of the rigid
  motions of the region's boundary with the part's free material (a face
  it holds, or its edge inside the part), those its free motions leave,
  each weighed over that boundary by the hold's own strength; its
  material is the part's, free to deform. A pinned column's end face so
  spreads as it is squeezed, as a pin lets it. `(rigid)` holds the region
  rigidly instead but for its free motions, as one body. Take the default
  for a support whose region is only where the load goes in (a pinned or
  rolled end, a bearing seat, a hinge line): it stiffens nothing that is
  not there. Say `(rigid)` where something stiff really holds the region's
  shape (a steel pin in a plastic bore, a clamped lug that turns), since the
  default lets the region deform as the part's free material would. The solve reports how far each region moved along its free
  motions (`free_motions`: radians, mm), of its mean motion for a
  distributing one.
- Declare `(mirror sym x 0)` on a shape cut at x = 0. Use x, y or z
  and a coordinate in mm. Every load, material and physical restraint in
  the model is reflected. Applied forces, reported reactions and compliance
  refer to the modeled portion. Gravity must lie in each mirror plane and pressure must
  act on exposed faces. Use the whole part for a loading that does not
  reflect. A roller on a whole part can be written
  `(fix floor (- z 0.01) (dofs z) (support))`.
  An older zero-displacement layer over a complete flat face can describe
  either condition at any layer depth. A second unmarked fix on a declared
  cut must be removed or marked as a physical support. A support meeting the cut may slide
  parallel to it or turn about its normal. These shared motions are
  held in a buckling mode antisymmetric across that cut. Other free
  motions there require the whole part. An interpretation
  contradicted by the loads (such as a floor roller under transverse
  gravity) is discarded. Where possible ink checks each interpretation's
  buckling load. Declaring the intent removes the extra solves and lets
  a direct claim use the intended factor.
- A fix holds all of its region, both ways: ink has no support that only
  pushes, so a shelf fixed on two strips pulls where brackets would lift.
  A fix without free motions holds the region rigidly however thin it
  is, so a strip held along one axis still resists a turn across it: a
  shaft held on lines through its axis takes torque and bending moment.
  Leave a bearing's turns free with `pivot`.
- Mechanisms: before solving, an elastic study finds the rigid motions
  its supports leave free. Loads that drive one have no static solution:
  the study fails, naming the motion and the moment or force driving it.
  Balanced loads leave it at rest, with a warning. `benchmarks/36-pin-supports.ink`
  shows each.
