# Examples

Small models, each showing one thing ink does. Each holds every claim it
states except those it marks as meant to fail (`; expect-status: NAME
fails`). Each line names the section of
[the worked example](../docs/guide/walkthrough.md) it goes with.

- [build-order.ink](build-order.ink): an assembly order checked step by step: each part comes in clear of those placed before it, and a screw too long for its blind tap fails. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [chain-drive.ink](chain-drive.ink): a roller chain on two sprockets, every link a part, checked for interference as the chain runs. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [curves.ink](curves.ink): profiles, curves and sweeps from the standard library. ([3. A part](../docs/guide/walkthrough.md#3-ask-for-a-part))
- [desk-lamp.ink](desk-lamp.ink): a lamp whose centre of mass stays over its base on a 10° slope with the arm swung any way, and a long-armed copy that tips: measured claims, holding and failing. ([4. Answers](../docs/guide/walkthrough.md#4-read-the-answer))
- [frames.ink](frames.ink): a pin and a bearing seated by frames, with claims that fail if the pin floats or the bearing floats or sinks. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [geared-drive.ink](geared-drive.ink): a motor, a pulley on its shaft and a 15:45 gear pair, the teeth checked over the motor's angle. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [gyroid-heat-exchanger.ink](gyroid-heat-exchanger.ink): one gyroid wall splitting a barrel into two channels that never meet. ([3. A part](../docs/guide/walkthrough.md#3-ask-for-a-part))
- [linkage.ink](linkage.ink): a four-bar linkage whose angles are solved from its loop at every moment. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [machining-access.ink](machining-access.ink): what a 3-axis mill reaches from its setups, and a T-slot it can't. ([6. Can it be made?](../docs/guide/walkthrough.md#6-can-it-be-made))
- [motors-on-radius.ink](motors-on-radius.ink): four motors on an X frame's arms, with claims that fail if one is placed wrong. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [mould-undercut.ink](mould-undercut.ink): an injection-moulded half checked for undercuts along its pull, and the side action that frees one. ([6. Can it be made?](../docs/guide/walkthrough.md#6-can-it-be-made))
- [name-ornament.ink](name-ornament.ink): a name ornament made to order from the customer's text, its letters in a script font sized to fit, with a claim for every name it allows. ([3. A part](../docs/guide/walkthrough.md#3-ask-for-a-part))
- [nameplate.ink](nameplate.ink): nameplates sized from their text, engraved and raised. ([3. A part](../docs/guide/walkthrough.md#3-ask-for-a-part))
- [phone-stand.ink](phone-stand.ink): a part placed on a tilted face of another. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [print-checks.ink](print-checks.ink): an FDM part checked against published rules for walls, overhangs and moving clearances. ([6. Can it be made?](../docs/guide/walkthrough.md#6-can-it-be-made))
- [sheet-metal.ink](sheet-metal.ink): sheet-metal parts designed flat, their folded sizes checked against hand calculation. ([3. A part](../docs/guide/walkthrough.md#3-ask-for-a-part))
- [sheet-metal-lap-joint.ink](sheet-metal-lap-joint.ink): a bought bracket recreated from its maker's drawing. ([3. A part](../docs/guide/walkthrough.md#3-ask-for-a-part))
- [sheet-metal-varying-holes.ink](sheet-metal-varying-holes.ink): holes placed by a toleranced size, checked over its whole range. ([7. Sizing](../docs/guide/walkthrough.md#7-sizing-it))
- [sliding-drawer.ink](sliding-drawer.ink): a drawer on a pair of slides, following the slide as it extends. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [sls-hollow.ink](sls-hollow.ink): a hollowed SLS knob, its escape holes and enclosed voids checked. ([6. Can it be made?](../docs/guide/walkthrough.md#6-can-it-be-made))
- [spring-press.ink](spring-press.ink): a spring whose shape compresses as the press moves. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
- [two-link-arm.ink](two-link-arm.ink): a two-link arm on servos, each placed on the one before, its tip checked where the angles put it. ([5. Motion](../docs/guide/walkthrough.md#5-an-assembly-that-moves))
