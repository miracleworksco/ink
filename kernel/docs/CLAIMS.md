# Requirements, facts and ratings

How a design says what must be true of it ("these two never touch",
"this motor is within its ratings") and gets an answer that shows its
work. This page says what the pieces are and what each may know. It describes the code as it is; change the code and this page
together. The words a model author writes are in the reference
([IR.md](IR.md), Facts and requirements).

## What a design writes

Two words. `fact` states a number about a part, in the instance that
places it. `require` asks a question, at the top level.

```lisp
(assembly opener
  (instance motor motor.drive
    (fact radial-load shaft door-weight)                 ; on its frame `shaft`
    (fact torque      shaft (* door-weight pitch-radius))
    (fact speed       shaft 22rpm)))
(require motor-ok (within-ratings opener/motor))
```

A fact's value is any number of params and units. Its unit is what ink's
arithmetic makes of it (a force times a length is a torque). `require`
asks one of two things: `(within-ratings PATH)`, every rating of every
part at or under PATH (less the parts left out on purpose:
`(within-ratings PATH (excluding P ...))`); or a claim (the reference's
Claims), with its clauses:

```lisp
(require clear (disjoint (shape-of arm) (shape-of post)) (accuracy 0.01))
```

A claim is answered by the prover (below), and no two requirements share
a name. `ink require` asks what the file itself requires, of ratings and
of claims, not what the files it imports do (`--imports`: those too).

A fact's name is not a free word: a place takes the facts its ratings
and formulas read, and no others. A fact stated where nothing takes its
name (`sideload` for `radial-load`) changes no answer, so it is said: a
`warning` under each requirement over its place, naming what that place
does take. `ink check` shows the same before anything is asked: `facts`
(each fact stated, its value and unit as read, its line, and whether
anything reads it) and `takes` (by place, the names taken there).

A fact need not be typed where the model works it out: a motor's torque
and speed, the load on its shaft, and each bearing's load and speed once
the design names what holds what (`(supports PATH S ...)`; sim's
DYNAMICS.md, Ratings, opens with the whole machine). `ink check` lists the loads a design applies and
the supports it names.

## The command

`ink require FILE` answers each requirement and, under it, what it rests
on, down to the line that states each number: in JSON, as every `ink`
command answers (a row each: `verdict`, `why`, `strength`, `covers`,
`from`, `warnings`, `notes`, `left_out`, and `rests_on`, the rows under it), or in
words with `--plain`:

```
motor-ok: NO, as stated  (garage-door.ink:13)
  the 4 ratings its parts carry, each checked:
  side-load on opener/motor/rotor: NO, as stated  (17HS19-1684S-PG27 gearbox shaft-radial-load-N-max)
    the size of radial-load on opener/motor/rotor shaft: 222.411 N  (garage-door.ink:9)
    at most
    100 N
  thrust on opener/motor/rotor: not settled: the design lacks something  (17HS19-1684S-PG27 gearbox shaft-axial-load-N-max)
    the size of axial-load on opener/motor/rotor shaft: not settled: the design lacks something  (ink effort on opener)
      note: nothing states it, and the model does not work it out: the file declares no (gravity gx gy gz): ...
  ...
```

**The headline.** A requirement of ratings never answers a bare yes
while anything is uncovered, assumed or unchecked. Its first line is its
verdict with what the verdict is of, in words and as the row's
`headline`: `YES, estimated, at 142 values of t from 0 to 2: for the 10
ratings checked; 1 thing not covered; placement not checked; not all
yes`. An open one says why at the top: how many of its ratings have no
answer, and the first fact wanted.

**Placement.** A requirement answers for the parts as they are joined,
not for where they sit, so where a file asks ratings and states a motion
`ink require` first runs `ink interfere`'s check at the start of the
motion (boxes, then the prover on the pairs that may meet: about a
second on a small machine). Parts that overlap in solid are not covered,
said first and by pair, and the file is not all yes (`placement`); none
found, the first line says so and the headline reads `no clash at the
start pose`. Then the parts that touch nothing are named (`LOOSE`;
`placement.loose`, by group): those no part outside their group comes
within 1 mm of, proven by a clearance, and no relation joins (a contact,
a mesh, a belt). A pin through a running fit touches; a plate lifted off
its standoffs, or a lid whose pin was deleted, does not. A loose part
changes no verdict (a rig may hang a load where its arm would be), and
what holds a part is not seen: a plate where one of its screws was still
touches. Left to `ink interfere FILE` and a render: the rest of the
motion, and what holds each part where it touches (`not_checked_here`).
It is skipped, saying why (`PLACEMENT NOT CHECKED`), for a file with
over 60 parts (a robot of 800 takes minutes), one that states no motion
(parts asked as stated, not a machine assembled), and one with no root.

A yes over a path says what it covered and, first, what it did not
(`NOT COVERED` in words; fields of its row and of the answer):

- **Parts.** `parts`: those under the path that carry a rating
  (`rated`); those the solver's model loads that carry none
  (`not_rated`, each with why it counts), which nothing checked; those
  of them whose record says its maker states no limit (`no_limit_stated`,
  in the record's words: a body's `(attr no-rating "...")`), which
  nothing could check, said apart ("NOT COVERED, and nothing could");
  those the design left out on purpose (`left_out_on_purpose`); and
  how many other bodies carry none. The line: a part counts where the
  design's model puts a demand on it that a maker rates: it drives the
  motion (an actuator's joint), a coupling passes force through it (a
  screw and its nut, a belt's pulleys, meshing gears), or the design
  names it as holding something (`(supports ...)`). A printed bracket, a
  plain shaft, a carriage are none of these: the model works out no
  load a datasheet could be held against, and a study answers for them.
  With no solver (the kernel's own `ink`) no part is known to be loaded.
  A rated part that holds a coupling's pull with no rating that reads a
  load across or along it (a motor rated for its torque, under its
  belt) is said too (`loads_no_rating_reads`); so is the bending and
  side force what an actuator carries puts across its output, with the
  numbers, where the bending is over a quarter of the torque it turns
  with and its record rates none.
- **The whole design.** `parts_under_no_requirement`, first in the
  answer: each part that carries a rating or that the model loads,
  under the assemblies the requirements name, that no requirement of
  the file asks over and none leaves out. Requirements asked part by
  part are each answered, and the file is not all yes while one is
  left. `(excluding P ...)` leaves out what it names, a rated part
  too, and is said as written.
- **Time.** `span`: the span of t the motion's run was over, said with
  "at N values of t". `motion_goes_on`: the written motion is not over
  at the span's end: a coordinate is not at rest there and does not
  repeat the span (a cycle repeats); or it is at rest and its drive's
  formula moves it again later, found by the kernel's bounds of its
  rate over all t after (said with when); or nothing proves it at rest
  for good, and it is said that this is not known. A profile with a
  duration, a clamp or a select on t is read to its end.
- **Driven and held.** `not_driven_or_held`: what the motion needs
  that nothing in the design supplies or holds, where the solver's model
  finds it: an actuator that moves through the whole motion with no
  effort to speak of (it drives nothing: its load is on its case, not
  its output); a named support whose box is clear of every other part
  (it is in mid-air, and passes its load to nothing); a coordinate needing an effort no actuator gives (a weight
  on a slide nothing joins to the motor that should lift it: the efforts
  found leave it out), a part on a joint that does not move and is no
  actuator's, loaded along it (it would not stay), a body moved by a
  formula of t through no joint, a load on a part on no joint (it goes
  to ground). Each says what to write.
- **The books.** `books_do_not_balance`: where the solver's own
  accounting does not add up (what a joint's supports were given
  against what the joint needs held, a support's offered load against
  every share it was given, the actuators' power against the work the
  motion takes, each joint's reaction against the weight, loads, pulls
  and mass times acceleration of what it carries), saying which part,
  how much and which way. A
  requirement resting on it is not settled, whatever its ratings said.
- **Not checked here.** `not_checked_here`: where the file relates
  parts (a gear mesh), its rows are the relation's spacing and alignment
  alone; whether the parts clash where they meet (teeth at their phase)
  is `ink interfere`'s to prove, and the answer says so.
- **Assumptions.** `assuming` ("ASSUMING" in words): what the answers
  assume because the design does not say, each once: a coupling that
  passes power and states no efficiency loses nothing, so the efforts
  through it are the least the motion can need. An answer resting on
  one is an estimate. Stating the thing (one line, which the assumption
  names) sharpens the answer.
- **Loads.** `loads_not_applied`: a load an imported file states at its
  top level on a part this design places too. It is on that file's
  root, which the design does not place; one that travels with its
  assembly is written inside it.

The verdict stays what the rated parts say. While a requirement has a
part not rated (and not left out on purpose), one with no limit stated,
something not driven or held, a motion that goes on or a load not
applied, the file is not all yes, and the answer lists them
(`parts_not_rated`, `parts_with_no_limit_stated`, `not_driven_or_held`,
`motion_goes_on`, `loads_not_applied`).

It exits 0 where every answer is yes, 3 where one is not, 1 where the
file does not load. A yes shown at one value of something it reads (one
moment of a motion no `(t-span)` covers) says so and is not counted
as all yes (`yes_at_one_value`). Where a param declares cases (`(cases v
...)`), every requirement is asked at the params' values and at each
combination of their cases: its answer is the worst of them, says which
case that was and at how many values of the param it was shown; a
claim's row has how it came out at each (`cases`), and the answer lists
the combinations asked. `--name N` asks one, `--imports` what its
imports require too, `--res`, `--t-span`, `--param` and `--keys` as the
reference says (IR.md, Claims: Results). A requirement of ratings is
answered over a span, so `--t` is refused where one is asked.

The command is the kernel's (`ink-cli`), so every build of `ink` has it.
The kernel's own asks the core, the tables, the formulas and the prover.
A distribution hands `ink_cli::main` what makes its solver tools
(`ink_cli::Solvers`), and they are asked after those.

The claims nobody writes are requirements too: those a limit clause or
a relation makes (`limits`, `speed-limits`, `effort-limits`, `relate`),
each under its name (`effort-limits:root/j/arm`). A file with nothing
required answers nothing and runs nothing (`all_yes` null): a rated
motor placed in a sketch makes no requirement by being there. Loading
and `ink check` work out no fact or requirement; `ink check` lists them.

A claim's row carries the prover's evidence beside its verdict: its
`status` (`holds`, `fails`, `unknown`), a witness point, an undecided
cell, its hints, its value and accuracy, how it held over a motion (the
reference's Claims: Results says each). A call stopped at its time limit
gives the rows done.

## A fact

A fact is a value with its story (`ink_claims::Fact`):

- what it is about: a name and a place (`radial-load` on
  `opener/motor/rotor shaft`). A part's place is its path; a frame's is
  the part's path, a space, and the frame's name;
- a number with its unit, or a verdict, when it is a statement;
- a band: the low and high its source stands behind. A number with no
  band is exact;
- how strong it is, over how much of each variable it was shown, and
  what it left out (below);
- where it came from (a line, a record's field, a tool and its model),
  the facts it rests on, and notes.

## The four verdicts

| Verdict | Means | Reasons |
|---|---|---|
| yes | it holds | |
| no | it does not | |
| not settled | something could still decide it | the design lacks something; a source could not decide; the band reaches across the limit; sources disagree |
| not available | nothing here can answer it as asked | nothing states it or offers it; outside what its source covers; it rests on itself; it checks nothing; it compares different units |

Which reason an open answer has is said by the code that could not
answer, never read out of its words: a solver refusing says whether the
design lacks something (a moving body with no mass:
`ink_lang::Refused::lacks`) or the question is outside what it covers (a
body its model leaves out: `Refused::outside`); one that says neither
could not decide.

A yes or a number that rests on a fact that is not settled or not
available is that too, for the same reason. A no stands on the one
thing that fails.

**Strength**, weakest first: estimated (a method that can be wrong),
as stated (it rests only on numbers the design or a record states;
nothing was proven), proven to its accuracy, proven. An answer is as
strong as the weakest thing it rests on.

**Coverage**: for each variable a fact reads (time, a ranged param),
whether it was shown at one value, at so many, or at every value.

## How a requirement is answered

One loop (`ink_claims::Claims::answer`), for any fact asked:

1. **Stated.** Every fact the design states about it. Two that differ
   leave it not settled: sources disagree, each shown with its number
   and its line. A stated fact is also held beside what a tool that
   measures the design's own model works out about the same thing, in
   the same unit (`Tool::measures`: the motion solver, the structural
   solver; each says which way its facts are the worse, more for a
   demand), so that a typed number never makes an answer better than the
   model's. By size: stated at or past the model's most, it stands (an
   allowance is the designer's to state); stated short of the model's
   least, the model's fact is a second source and they disagree; stated
   inside the model's band, the fact is the model's, its band from the
   worse of the two least to the model's most, and it says what was
   typed.
2. **Worked out by a rule.** A rating is a rule: each term at most the
   next (`x <= limit`), each term a number, a fact, or a fact's size
   (its number without its sign). `within-ratings` is a rule over the
   ratings under a path: no if one fails, else open if one is open, else
   yes; with no rating there, not available.
3. **Tools, in the order given**, the cheapest first, until one settles
   it. Nothing runs until a fact needs it. Where a comparison is too
   close to call on what one tool gave, the next tool that can settle
   that fact gives it, and the comparison is worked again.

For a requirement of ratings `ink require` plugs, in this order: the
tables parts carry, the formulas parts carry, the motion solver, the
structural solver. A claim does not go through the loop: the prover
checks every claim of the file in one run, and each answer is made a
fact of the core's kind (`ink_ratings::fact`: its verdict and reason,
how strong, over how much of what it reads, what it left out).

## A tool

A tool offers facts (`ink_claims::Tool`, one function): asked about a
fact, it says nothing, or gives it, or gives it open with a reason. It
may ask the core for other facts while it works, and must say which it
rested on.

| Tool | Crate | Offers |
|---|---|---|
| tables | `ink-formula` | a limit read off a maker's curve at another fact (a stepper's pull-out torque at a speed). A drawn curve's points (`"drawn": E`, their reading error) are read as the line between them, E either side, an estimate that says so; points that are measurements each alone are, between two, anywhere between the two. Below its first point it is that point's value, an estimate that says so, where the part says the curve may be read so (`"below": "first"`: a stepper whose holding torque is no less); elsewhere outside the curve it is not available, saying where the table begins or ends and why |
| formulas | `ink-ratings` (`Derived`) | a fact a part works out from other facts about it (a bearing's rated life) |
| duty | `ink-ratings` (`Duty`) | what is of a whole motion and of no moment: `mean-speed` and `duty-load` at a place (ISO 281's equivalent load of a varying duty), from `speed` and `equivalent-load` there at every moment; `rms-torque`, the root mean square of `torque` there over the motion's time (what a continuous rating is of); with no motion, what is stated, and no speed where none is |
| motion | `ink-effort` (in sim) | at each moment of the motion: an actuator's torque or force and its speed; the load across and along a joint's axis on the part it moves; each named support's share of a reaction, and its speed. Each from the least it needs: a speed from the written motion alone, with no mass and no gravity; a torque, force or load where the moving bodies it rests on have a mass and a gravity is declared, and otherwise "the design lacks something", naming the bodies, for that one alone. It measures the model: a stated fact is held beside it |
| structure | `ink-study` (in sim) | a part's peak stress and deflection under a study's loads |

**Tools do not import tools.** A solver tool depends on the core and on
its solver, never on another tool or on the reader. One exception
stands: the reader (`ink-ratings`) imports `ink-formula`, to hand it the
tables it read.

## What a part carries

**Ratings**, as data in its body's `ratings` attribute, written by the
part's generator: by name, the terms in order and where the limit is
from.

```json
"side-load": { "ordered": [{ "read": "radial-load", "on": "shaft" }, { "number": 100, "unit": "N" }],
               "sizes": true, "from": "17HS19-1684S-PG27 gearbox shaft-radial-load-N-max" }
```

A term is a number with its unit, a fact to read (on the part, or `on`
one of its frames), or a table read at a fact. `sizes` reads each fact
without its sign. A fact a rating compares as it is (a required life
against a rated one) may say the least it can be, `{"read":
"required-life", "least": 0, "unit": "h"}`: one stated under it is no
requirement (a slip of sign would pass any limit), and is not
available, saying so. A number, a table's points or one coordinate of a
point may be `{"record": [KEY, ...]}`: read from the part's record (its
body's `spec`) at that path, so a limit the record holds is said once,
not copied beside it. A table may say `"below": "first"` (read below its
first point as that point's value) or `"below": "none: WHY"` (not
available there, and WHY is said); and `"drawn": E` where its points are
read off the maker's drawn curve to within E of its y.

**Formulas**, as ink functions in its body's `derives` attribute: the
fact named is the function of the facts its parameters name, each read
in the unit written after it.

```lisp
(vec "rated-life"
  (lambda (equivalent-load speed) ...)   ; ISO 281, in the bearing family's file
  1N 1rpm)
```

The function is ink's own arithmetic, so the fact has ink's units and
bands. It is worked out when asked, from what the core says of each
fact it reads, wherever that was stated. A fact in another unit (a
speed in Hz, a bare number) is refused, never read as if it fitted. A
formula says where it holds by giving `(vec VALUE OK WHY)`; elsewhere
the fact is refused and says WHY.

```
life on root/stage/out-bearing: YES, as stated  (its record's load-ratings C-N, ISO 281 basic rating life)
  required-life on root/stage/out-bearing: 36000000 s  (planetary-gearhead.ink:115)
  at most
  rated-life on root/stage/out-bearing: 3318599215293.431 s  (bearing-catalog.ink:72)
    equivalent-load on root/stage/out-bearing: 66 N  (bearing-catalog.ink:62)
      radial-load on root/stage/out-bearing: 66 N  (planetary-gearhead.ink:112)
      ...
```

Which parts carry ratings is a part library's to say. A family whose
ratings are the same for every member may say them once, as data its
family file reads, with each limit a formula of nothing that reads the
member's record.

## Refused, on purpose

Each is an error or a refusal that says why, never a guess:

- a fact with no unit, or another quantity's, against a rating;
- a frequency against a turning rate (a cycle is not a radian);
- a fact whose value is no number (a division by zero, a number outside
  ink's range), or that reads a value ink works out at a run;
- a rating with the name of a fact at its place, and a design that
  states a fact its part works out (a formula's, or a duty's:
  `mean-speed`, `duty-load`, `rms-torque`): one would answer for the other;
- a stated fact under the least what reads it takes (a negative
  required life);
- a stated fact the model works out more than: sources disagree;
- a life at a still pose: a life is of a duty.

## What each crate may depend on

| Crate | Is | May depend on |
|---|---|---|
| `ink-claims` | the core: facts, rules, the loop | nothing |
| `ink-formula` | the table tool | `ink-claims` |
| `ink-lang` | the language, which reads `fact` and `require` as data; and the prover, which answers claims | not on any of these |
| `ink-ratings` | the reader: a loaded design as the core's facts and rules; the formula tool; a prover's answer as a fact (`fact`) | `ink-claims`, `ink-formula`, `ink-lang` |
| `ink-effort`, `ink-study` (in sim) | solver tools | `ink-claims`, `ink-lang`, their solver |
| `ink-cli` | `ink require`: plugs the tools and answers | `ink-claims`, `ink-ratings`, the kernel; no solver |
| `ink-ext` (in sim) | hands `ink require` its solver tools (`solvers`), as it hands `ink` its packages. A distribution with a part library does the same | all of the above |

The core knows no geometry, physics or units: a unit is a word it only
compares. The language knows nothing of the core.

## Not built

- A fact about no single part (a bolted joint's preload, a chain's
  tension between two parts): `fact` is written only in an instance.
- A claim asked at the moments a motion's run read (the core's `over`):
  a claim that reads t is the prover's, over the file's `(t-span)`, and
  with none it is asked at one moment and says so.
- A formula on anything but a body: one is not found through an
  assembly's frame, as a rating read `on` a frame is.
- A formula that gives an estimate or a band of its own (a table does).
- A fact that reads a value ink works out at a run (a solve, an
  estimate).
- Ratings on most parts: a part library says which of its families
  carry them.
