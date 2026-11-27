# Standard material inputs

`materials.ink` is plain IR data. Use `(material NAME)` to inspect a
material, or `(attr material "NAME")` to give a body its properties.
`material-sources.json` holds source IDs with title, maker, revision, URL,
retrieval date and hash, plus the kind and conditions of each property.
This follows the catalog's shared source-manifest convention. Source
files themselves are not vendored. Every scalar's adjacent comment either
cites its source ID or explicitly says `UNSOURCED`.

E and strengths are MPa, density is g/mm³ (literals use kg/m³), thermal
conductivity is W/(m K), specific heat is J/(kg K), viscosity is Pa s.
Physical properties are typical reference values; temperatures can differ
between keys and are stated in the evidence. Ranges and coupon spreads
are described in comments, not propagated as numerical uncertainty.

Yield defaults use the specification minimum for the named condition when
one is available. Product-form and thickness limits are part of the
condition; these numbers are not qualification of arbitrary purchased
stock, welds, heat-affected zones, temperature, fatigue or creep.
Maker typical strengths remain explicitly `datasheet typical`; these are
not guaranteed minimums. A strength read by `get` whose `strength-kind` is
not a specification minimum makes any claim that reads it an estimate,
naming it (`typical_strengths`): it can say a design meets a typical
coupon, never that it is guaranteed to.
Typical material data cannot establish
a guaranteed real-world safety factor.

PLA and PETG name the maker, horizontal coupon direction and FFF process.
Other polymer names specify stock shapes, not printed parts. Read the
source's coupon process and the per-property spread; do not invent an
across-layer knockdown.

Unsupported yield/tensile/compressive/bending fields are absent: a reader
that needs one must ask for a qualified input rather than getting an
optimistic default. `glass` and borosilicate expose `bending`, not `tensile`.
Remaining unsourced physical inputs are marked `UNSOURCED` in
std/materials.ink.

## Names and compatibility

Each old name below remains an accepted, silent alias. Material results
retain the requested `name`, and add the resolved `canonical-name`,
`condition`, and `evidence` manifest reference. No alias warning is emitted:
existing catalog models must continue to load without new warnings.
Case/spacing variants and `other-names` remain spelling errors with hints.
Aliases target canonical definitions directly, without chains.

Each alias, its canonical name and its condition are the `(def ALIAS CANONICAL)` lines and `condition` strings of [materials.ink](materials.ink).
