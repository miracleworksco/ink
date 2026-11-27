## Topology optimization

`(topology ...)` turns a static elastic study's shape into a design
domain. It searches for a stiff layout (low compliance) under the
study's supports and loads within a volume budget. The search is local,
from a uniform start: it finds a layout no small change improves, not
the stiffest there is.

### Quick start

A bracket laid out in a plate:

<!-- example examples/layout.ink -->
```lisp
; The stiffest bracket that uses 40% of a PETG plate, held at a wall and
; loaded at the bottom of its free end.
; ink topology examples/layout.ink: the layout found and its compliance
; (load times deflection, N mm: lower is stiffer) at each step;
; ink mesh examples/layout.ink --target bracket --out bracket.stl to print it.
(import "std")

(def domain (move (box 60 2 24) 30 0 0))         ; x from 0 to 60, z from -12 to 12
(study hang (shape domain) (elastic (material petg))
  (fix wall (- x 1))
  (load tip (force 0 0 -100N) (region (max (- 58 x) (+ z 9))))
  (res 1))
(topology layout (study hang) (volume-fraction 0.4))

; The layout as a part: material where its density is over 1/2.
(body bracket (shape (max domain (study-field layout solid))))
```

Run from sim/ (seconds to a minute on a laptop's CPU):

<!-- run -->
```sh
ink topology examples/layout.ink
```

```json
{"results": {"iterations": 52, "settled": true, "compliance_nmm": 459.9,
             "volume_fraction": 0.400, "compliance_start_nmm": 3533.4}}
```

- `compliance_nmm`: the load times the deflection it causes, N mm; lower
  is stiffer. The layout's is 460, from 3533 for the plate at 40%
  density everywhere: about 7.7 times as stiff for the same material.
- `settled`: the layout stopped changing, after 52 updates
  (`iterations`, of at most 60; `history` has a row for each). It is a local search: a
  layout no small change improves, not the stiffest there is.
- `accuracy: null`: an estimate, and its compliance is the density
  model's, not the part's. To trust the part, mesh it and give it a
  study of its own, as examples/topology.ink does.
- `part_buckling_factor`: the part's linear buckling factor, low to high
  (the loads times it buckle it), from the part solved as a study of its
  own (the density model runs no buckling analysis: grey cells buckle
  where no part does). Where the loads may reach it (its low end at most
  1) the answer warns, as a study does; where it is not known, it says
  why.

The layout is the field `bracket`'s shape, so it meshes, renders and
measures like any other:

<!-- run -->
```sh
ink mesh examples/layout.ink --target bracket --out bracket.stl
```

### Forms

```lisp
(study load (shape domain) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 1))
  (load tip (force 0 0 -200N) (region r))
  (res 1))
(topology layout (study load) [(volume-fraction 0.35)] [(iterations 60)] [(filter 3)] [(keep field)])
(def part (max domain (study-field layout solid)))   ; the optimized part, as a field
```

- The study comes first, and is a static elastic one whose fixes hold
  still (no `(modes ...)`, `(transient ...)` or `(move ...)`). Its grid
  (`res`) is the layout's; `ink topology FILE [--name NAME] [--res MM]`
  runs it and answers its history.
- The method is SIMP with the three-field scheme:
  - a design density per cell;
  - a density filter (radius `filter` in mm, default 1.5 cells; at least
    a cell is used);
  - a smoothed Heaviside projection, sharpened as the iterations proceed
    (its steepness doubling from 1 to 16), that drives densities to solid
    or void;
  - stiffness density³, updated by the method of moving asymptotes
    (Svanberg), each step meeting the volume.
- The steepness doubles once no cell's density moves by 0.01 in an
  iteration, or after a tenth of `iterations` (default 60, a whole
  number; a filter that is not a positive length is refused too). The run
  settles once none moves by 0.01 at the steepest with the volume within
  its budget (to 1e-3 of the domain), or stops after `iterations`. The
  answer's `settled` says which, and `iterations` how many updates it
  made; a run that did not settle warns: give it more iterations (unless
  what stays solid fills the budget, which no number settles). The summary
  (`compliance_nmm`, `volume_fraction` (default 0.3), `last_change`: the
  largest move from the design before) is the returned design's, the
  history's last row (the history has a row per design solved: the first,
  then one per update, each with the projection's steepness).
- Material inside `keep` stays solid, and so does a cell whose material
  a fix holds all of (its hold is in the cells beside it, so its own
  material would hold nothing the optimizer sees). Both count toward the
  volume budget: if they reach the volume fraction, nothing is left to
  optimize (the answer warns). Elsewhere supports and loads act on
  whatever material is left under them, so the optimizer keeps material
  where they need it.
- The answer is an estimate: `accuracy` null, `estimate` its cells,
  `method` the search. Its compliance is the density model's (grey
  cells counted at density³), not the extracted part's.
- `(study-field layout density)` is the density field. `(study-field
  layout solid)` is 1/2 - density, which is negative where material was
  kept. The density cells are the controls of a cubic B-spline with
  continuous first and second derivatives. The field fades to void beyond
  the cells, so the solid closes without a seam at the edge of the grid.
- The result is an ordinary field. Fillet it, measure it, export it, and
  verify it with a normal study of the extracted part
  (`examples/topology.ink`, `benchmarks/15-topology-bracket.ink`):
  its claims are the part's, decided as any study's (Study claims), so a
  member the cells glue to another, or too thin for them, leaves them
  unknown.
