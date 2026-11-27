# Dynamics: validation and workings

How ink-dynamics' answers (DYNAMICS.md) are checked against references,
and the numerical choices behind them. Tests named here are in
ink-dynamics/src (validation.rs, answers.rs, modes.rs, springs.rs),
ink-ext/tests/dynamics.rs and ink-ext/tests/cli/dynamics.rs (the command line, with masses ink
measures). The library tests set their masses, so they check the
dynamics apart from mass measurement.

## Inverse dynamics

- validation.rs: P1 compound pendulum, P3 four-bar torque table, P5 lead
  screw motor torque, P6 50:1 geared arm (an ideal ratio with reflected
  rotor and payload inertia), P8 a wheel spinning while it steers (the
  gyroscopic moment I ω_s Ω on its axle bearings), P9 a plate turning
  about its diagonal (its bearing couple, ±31.153846 N over the bearings'
  500 mm), P10 Modern Robotics' 3-link chain with a tip wrench as a load
  (τ to its published digits), P4 belt force. Each to its reference's
  printed digits.
- answers.rs, by hand: a belt drive's bearings locked and held
  (`a_belt_drives_pull_on_its_shafts_by_hand`), a spur pair's
  (`a_spur_gear_pairs_bearings_by_hand`), a bevel pair's and a rack's
  (`a_bevel_pairs_forces_by_hand`, `a_rack_and_pinions_force_by_hand`), a
  screw's thrust (`a_lead_screws_thrust_by_hand`); a lossy belt
  (validation.rs `a_lossy_belt_drives_its_wheel_by_its_efficiency`).
- On the catalog's CoreXY each motor shaft carries its two spans, square
  to each other, √2 × 15 lbf = 94.361 N with the head moving steadily,
  by hand: ink's claims on it are unknown, the belts' own mass summed into
  each shaft's reaction (DYNAMICS.md, Reactions).
- Carried bodies: Gauss–Legendre takes two nodes on a line (exact) and
  four per quarter turn of an arc (within 5e-6 of the arc's term). Left
  out, the section's offset from the path and turning round an arc is
  about 2|e|/r + (e² + k²)/r² of a wrap's share (e the profile's centre
  across the path, k its section's radius of gyration): 0.1% of a CoreXY
  belt's.
- The midpoint RMS is 0.22% under a cycloid's exact RMS on its 1 s
  stroke.

## Steps

- Rotation is held within 1% of its unit entries, about 0.28 rad a step.
  Resolving a spin needs steps of at most √(0.08 / (ω² + |α|)); halving
  lands each between that and half.
- A group's departure is bounded over a step by the frame's motion and the
  spin's rate, acceleration and partials. On a reaction, the mass off its
  axis h moves the force by at most h (|L| + |W|²), and the moment about
  the axis by h |a₀ − g| plus its inertia not round across the axis times
  |L across| + |W| |W across|, and the rest of it times |L| + |W|² (W, L
  the spin's and its frame's turning rate and acceleration over the step,
  a₀ the frame's acceleration). A member's departure below 1e-9 of its own
  measure is none.
- A round roller, or a keyed rotor given its datasheet inertia, at 6000
  rpm over 2 s takes 16 steps. The catalog CoreXY's 12 s Lissajous takes
  196 steps (1675 resolving its rotors, 0.59 g mm off their axes;
  validation.rs `p16_corexy_lissajous`). A level ½" hex shaft spun up to
  6000 rpm carrying a hub, a sprocket and a collar (29 g mm off the axis,
  0.29 N mm of a 47 N mm peak) takes 19; held at 6000 rpm against ±400 N
  mm, 16, and its claim holds (answers.rs
  `a_hub_sprocket_and_collar_riding_a_hex_shaft_decide_at_6000_rpm`).
- A 6000 rpm disc with a 5 mm hole 30 mm off its axis, whose bearings see
  73.0 N turning, decides its bearing load claim in 18 steps, and its
  tilting moment claim (bound 1.81 N mm) as well.
- An actuator limit dropping to 5 N for 0.2 ms, 2e-10 s or 2e-12 s under
  a 10 N push fails at the outage (answers.rs
  `a_limits_own_jumps_are_steps`), and a sprung drive's limit dropping to
  0.1 N for 0.2 ms fails by its spring's force there
  (`a_sprung_drives_limit_jumps_are_steps_too`; in `ink effort`'s main
  answer, ink-ext/tests/cli/dynamics.rs `a_sprung_drives_limit_outage_is_in_its_main_answer`).
- A 1 um, 1 kHz vibration of 1 kg, its steps at the budget, reads a peak
  of m a = 39.48 N and its RMS that over sqrt 2, `partial`
  (`a_loose_motion_is_read_by_every_aggregate_and_says_partial`).

## Modes

modes.rs, each to 1e-9 against its closed form: a belt axis (a carriage
between a stepper's pulley and an idler without inertia) against the
two-mass closed form (rotor J / R² on N_r T_h / R², the carriage on k_b =
EA (1/L1 + 1/L2)), and with the rotor held rigidly (1/2π) √(k_b / m)
(`a_belt_axis_rings_at_its_two_mass_frequencies`); with an idler of its
own inertia, the three roots of its characteristic cubic
(`a_belt_axis_with_an_idler_rings_at_its_three_mass_frequencies`); a
rotor geared rigidly 10:20 to a shaft a stiff belt joins to a flywheel
(`a_rigid_belt_on_a_released_shaft_reflects_its_ratio`); masses tied
unequally by a rigid line and a rack's carriage and pinion, (1/2π) √(k /
(m + J / r²)) (`a_rigid_line_ties_unequal_masses_together`,
`a_rigid_ratio_reflects_its_inertia`); a carriage bolted at 0 carrying a
pulley between two spans (`a_joint_supporting_a_springs_end_is_released`).

The parts catalog's CoreXY (not in this repository; the head held at
(hx, hy); Gates' LL2MGT06 EA, which the model states; LDO-42STH40-1684AC
motors at 22.5 N m/rad): the gantry along y, the head riding, 50.1 Hz at
the centre, 47.9 at the front, 54.0 at the rear corners; the head along x
108.5 to 126.6 Hz; the rotors about 335 Hz (a rotor alone on its holding
stiffness, 322.5); the idlers 730 to 3570 Hz. The catalog's tests
(ignored for time, ~150 s cold) check
all twelve, at the centre and the front, against a model written by hand:
the head along x, the gantry along y, both rotors and the eight idlers,
each span a spring of EA / L with L from its belts' route, stretched by
its length's change less the belt fed through it; the motors' record;
the masses `ink mass` gives. They agree to about 1e-6 (the test holds
them to 1e-5). With massless idlers (a
4×4) they run 0.1 to 1.7% higher. The masses are measured to 2%, so the
frequencies are known from them to 1%, and move within that as ink's
measuring changes.

On belt data: a generic 6 mm GT2 belt measured 14.9 kN at 10–20 N (RepRap
wiki, Belt Compliance Measurement); Wang et al. (Robotics 7:75, 2018)
give 21.9 to 24.2 kN if their specimen was one 6 mm strand (its width is
not stated), their belts stiffening with tension. Hence -12% to +43% of
Gates' figure. A stepper's holding torque measured with both phases on,
driven at a per-phase peak of the rated current, gives 1/√2 of ink's
stiffness; measured with one phase on, at a per-phase RMS (a TMC driver's
run_current under Klipper), √2 of it.

## Simulation

Against closed forms (validation.rs):

- P4's belt carriage at full voltage, x(0.1 s) = 57.520173 mm and 300 mm
  at 0.467399050 s, to 1e-6 (`p4_simulated`); P7's hoist, 163.892229 mm
  at 1 s (`p7_simulated`); P1's bar let go at 60° reaching −60 at its
  elliptic half period, its energy held to 1e-7 (`p1_released`).
- A lead screw of efficiency E raised by a constant torque matches α =
  (τ − m g L / 2πE) / (m (L/2π)² / E) at E = 0.49, 0.45 and 0.35
  (`a_self_locking_screw_lifts_by_its_efficiency_and_holds_by_friction`).
- 1 kg blocks on `th` and `(* 2 th)` along gravity accelerate at 3/5 g
  (`a_simulated_def_read_by_several_arguments_moves_them_all`).
- A PD command lifting 1 kg 360° on a 1 mm lead screw of efficiency 0.9
  comes to rest held by friction at 0.314 s; at efficiency 1 it settles
  (`a_simulation_names_friction_switches_and_negative_inertia`).
- 1 N for 0.1 s on 1 kg leaves 100 mm/s (`a_pulse_is_not_stepped_over`);
  a 100 Hz PWM command's 200 switches cost a few stages each.
- examples/simulated-arm.ink: a servo steered by a PD command and a
  gate on its own spring; a clearance claim over it takes 66 steps.

Integrator: Dormand and Prince's 5(4) pair (Hairer, Nørsett and Wanner,
Solving ODEs I, II.5) until its stiffness test (h ‖f(y₇) − f(y₆)‖ /
‖y₇ − y₆‖ past 3.25 on 15 steps running, DOPRI5's) finds the motion
stiff; then Radau IIA of order 5 (Solving ODEs II, IV.8, as RADAU5); back
once its Jacobian's spectral radius, found every 15 steps, falls to a
tenth of where Radau took over and the step times it would keep
Dormand-Prince stable by half (LSODA's and AutoSwitch's switch). A
non-stiff motion costs what Dormand-Prince alone does (1 kg let go 1 mm
out on a 40 Hz spring: a second in 14036 states,
`a_non_stiff_ringing_costs_no_more_than_it_did`); a stiff one takes no
tiny steps (1 kg on a 1.6 kHz spring pushed by 500 (1 − cos 2πt) N
follows the closed form to 2e-7 mm, a second in 761 steps, where
Dormand-Prince alone stopped at its stage budget at 0.234 s,
`a_stiff_spring_the_motion_does_not_ring_on_takes_few_steps`); a 40 Hz
ringing damped by 100 N s/mm for 10 ms runs its second within the stage
budget (`a_moment_of_stiffness_is_stepped_past`). Radau IIA alone fails
the last two (both stop at the 20000-stage budget at 0.87 s) and costs
1.4 to 1.9 times the states on the non-stiff references (p1_released,
claims over a simulated motion), so both methods stay. Springs are
stepped by Radau throughout: Dormand-Prince's stability function exceeds
1 on the imaginary axis (|R(1.55i)| = 1.004), so an undamped ringing near
the tolerance grows at its steps (0.1 µm at 335 Hz to 0.23 µm). Each
speed is held to its position's tolerance times its own frequency (from
how its acceleration changes with it at the start), or, once stiff, by
what it moves over the step. A step over which the simulated
coordinates' generalized inertia changes more than twofold is halved, so
a dead point is closed in on. A switch in t is a rate with no bound over
the step whose bounds stay wider than 10⁻⁶ and than its rounding (10⁻¹²
of its size: a 1.8° step after 5556 turns is a jump); it is found by
bisection, stepped across at 10⁻⁹ of the span.

Springs (springs.rs):

- A 1 kg slide on `(stiffness 100)` whose command steps from 5 to 0 mm
  rings as 5 cos ωt to 4e-5 mm, its energy held to 1.5e-5, at its modes'
  50.33 Hz (`a_mass_on_a_spring_rings_at_its_mode_and_keeps_its_energy`).
- modes.rs's printer axis commanded one full step, 1.8° (0.2 mm), follows
  the closed form of its two masses to 9e-6 mm, each joint's ringing given
  its mode (223.72 or 362.67 Hz) and that mode's share
  (`a_stepped_belt_axis_follows_its_two_mass_response`).
- A commanded motor dragging its sprung carriage lags it by −(m a /
  k_b)(1 − cos ωt), 3.8 µm at 278 Hz, to 2e-6 mm
  (`a_commanded_motor_drags_its_sprung_carriage_behind`).
- A 1000 mm arm stepped 1e-5 rad rings by its 10 µm at its tip for a
  second (`a_long_arms_small_ringing_is_kept`); a 1.8° step rings alike
  after 5556 turns (`a_step_rings_alike_whatever_its_commands_origin`).
- The catalog CoreXY, a 20 mm move along y at 5000 mm/s² and 200 mm/s
  from (0, −10), over 0.2 s: the gantry's following error peaks at 0.104
  mm 10 ms into the acceleration (the hand estimate 2a / ω² at 50 Hz:
  0.101), the head's along x at 0.006 mm; at rest the gantry rings by 6.4
  µm in its 50.39 Hz mode (98% of it) and the head by 2.3 µm at 114.29 Hz
  (its accelerations last 40 ms, two periods of the 50 Hz mode, so they
  nearly cancel: at 150 mm/s, 30 ms, it rings by 0.17 mm); its spans stay
  between 63.3 and 70.1 N (15 lbf pretension); its steppers' effort claims
  read up to 47 N mm, within their 24 V curve by 330 N mm. The stage budget
  carries it about a second.

The following error's band, measured: 1 kg released from a swing A, its
error against A cos ωt stayed within 68 tolerances at 335 Hz over a
second (A = 0.1 µm, the worst; 49 at 1 µm; 31 at 10 µm and 13 at 0.1 mm,
where the stage budget stopped them at 0.88 and 0.51 s) and 11 at 50 Hz.
Over more periods Radau IIA damps an undamped ringing about 0.2
tolerances a period: a 0.15 µm ringing at 335 Hz is gone by 2.4 s (150
tolerances off, at most 0.45 of its band), one at 50 Hz 147 off over 15 s
(`a_long_ringing_stays_within_its_widening_band`). Measured over up to
1000 periods; a longer run may pass it. A ringing of 10 tolerances is
smoothed away. The band's growth counts the fastest mode's periods
whether that mode rings or not, so a light, stiff part widens it far past
the error on a long run.
