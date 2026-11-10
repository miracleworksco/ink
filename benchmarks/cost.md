# What claims cost, against `ink check`'s estimates

`ink check`'s `cost_estimate` (crates/ink-lang/src/cost.rs) classes a
claim's time where what its work costs is measured. This is that
measurement: claims of real projects (name plates, a key rack, a coin,
lettering, a dispenser, bowls, a baseplate) and of benchmarks/, asserted
cold (`INK_NO_CACHE=1`) at 8 threads on a 20-core machine shared with
other work (load 20 to 90, so times run long), against what `ink check`
estimated for them at 8 threads. A claim asked at several cases is
compared on its first.

**Thickness, 58 claims.** The work counted: surface samples (area over
res squared, on the h, 2h and 4h grids) times the instructions an
evaluation runs near the surface. Time = 0.3 s + work / (3e6 per thread a
second), fitted on these same claims (not held out):

| estimated \ measured | under 1 s | 1-10 s | 10-60 s |
|---|---|---|---|
| under 1 s | 26 | 4 | |
| 1-10 s | 1 | 20 | 4 |
| 10-60 s | | 3 | |

46 of 58 in their class, none off by more than one. None took minutes,
so the classes past tens of seconds are extrapolated.

**Values, 130 claims:** all estimated and measured under a second.

**Proofs, 176 claims** (subset, within, disjoint, overlaps, clearance,
touches) took 0.1 ms to 8.8 s (median 30 ms), with little relation to
the size of their fields (rank correlation 0.32): of the 11 over a
second, 6 have fields under 300 instructions; of the 28 with fields over
5,000, 24 took under a second. A box inside itself at 1 µm cells takes
0.5 s. Their cells are where interval bounds cannot yet tell, which a
look at the shapes does not predict, so a proof's estimate gives its
budget and its field's size, and no class. `empty` claims are proofs
too (here, overhangs: up to 17 s).

**What estimating costs `ink check`.** Cold, best of two each, on three
lettering projects: 3.1 s against 3.4 s, 1.38 s against 1.39 s, 6.2 s
against 6.5 s; within the spread of single runs there (±30%). A body's
area comes from the survey `ink check` already makes; other shapes'
are sketched beside it.
