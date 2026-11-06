#!/usr/bin/env python3
"""Reference implementation of the ink core's operations (docs/CORE.md), and
the generator of the conformance cases every implementation is tested
against (conformance/cases.json).

This file is a second, independent implementation of the spec: small,
written from CORE.md rather than from the Rust code, so when the two agree
the spec means what it says.

Value rules: IEEE-754 double arithmetic, with three changes that make every
operation total and every implementation agree:
- one zero: -0 is read as +0 (so atan2, division and pow never depend on
  the sign of a zero);
- division by zero gives 0;
- any other undefined result (NaN in IEEE: inf - inf, 0 * inf, sin of
  infinity, ln of a negative, ...) gives 0; so does sqrt of a negative,
  which is the continuous sqrt(max(a, 0)).
Infinities are values: overflow gives +-inf, ln(0) is -inf, pow(0, b < 0)
is +inf.

Derivative rules: the derivative of the branch the value took, with each
formula evaluated under the value rules above (so a division by zero in a
formula gives 0). At an exact tie between branches (min, max, abs at 0) the
derivative is the average of the two sides, so it does not depend on the
order of arguments. Through an expression, the chain rule op by op: each
op's partials (from its operands' values) times its operands' derivatives,
summed in order (`jet`), so a tie's average is da/2 + db/2.

Run: python3 conformance/reference.py > conformance/cases.json
"""

import json
import math
import random
import struct
import itertools
import re

INF = math.inf


def canon(v):
    """The value rules' representative: NaN -> 0, -0 -> 0."""
    if v != v or v == 0:
        return 0.0
    return v


def div(a, b):
    a, b = canon(a), canon(b)
    if b == 0:
        return 0.0
    return canon(a / b)


def mul(a, b):
    return canon(canon(a) * canon(b))


def is_int(b):
    return math.isfinite(b) and b == math.floor(b)


def pow_(a, b):
    a, b = canon(a), canon(b)
    if a < 0 and not is_int(b):
        return 0.0
    if a == 0:
        if b < 0:
            return INF
        if b == 0:
            return 1.0
        return 0.0
    try:
        return canon(math.pow(a, b))
    except OverflowError:
        # |a|^b too large: the sign is negative only for a negative base and
        # an odd integer exponent.
        neg = a < 0 and is_int(b) and int(b) % 2 == 1
        return -INF if neg else INF


def exp_(a):
    try:
        return math.exp(a)
    except OverflowError:
        return INF


def trig(f, a):
    if not math.isfinite(a):
        return 0.0
    return canon(f(a))


def ln_(a):
    if a == 0:
        return -INF
    if a < 0:
        return 0.0
    return math.log(a)


def floor_(a):
    if not math.isfinite(a):
        return a
    return canon(float(math.floor(a)))


def hypot_(a, b):
    return math.hypot(a, b)


# ---------------------------------------------------------------- values

UNARY = {
    "neg": lambda a: canon(-a),
    "abs": lambda a: abs(a),
    "sqrt": lambda a: math.sqrt(a) if a > 0 else 0.0,
    "square": lambda a: mul(a, a),
    "exp": exp_,
    "ln": ln_,
    "sin": lambda a: trig(math.sin, a),
    "cos": lambda a: trig(math.cos, a),
    "tan": lambda a: trig(math.tan, a),
    "floor": floor_,
}

BINARY = {
    "+": lambda a, b: canon(a + b),
    "-": lambda a, b: canon(a - b),
    "*": mul,
    "/": div,
    "min": lambda a, b: min(a, b),
    "max": lambda a, b: max(a, b),
    "atan2": lambda a, b: canon(math.atan2(a, b)),
    "pow": pow_,
    "emax": lambda a, b: hypot_(a, b) if a > 0 and b > 0 else max(a, b),
}


def value(op, args):
    args = [canon(v) for v in args]
    if op in UNARY:
        return UNARY[op](*args)
    if op in BINARY:
        return BINARY[op](*args)
    if op == "select":
        c, a, b = args
        return a if c < 0 else b
    raise KeyError(op)


# ----------------------------------------------------------- derivatives


def d_unary(op, a):
    if op == "neg":
        return -1.0
    if op == "abs":
        return 1.0 if a > 0 else -1.0 if a < 0 else 0.0
    if op == "sqrt":
        return div(0.5, math.sqrt(a)) if a > 0 else 0.0
    if op == "square":
        return mul(2.0, a)
    if op == "exp":
        return exp_(a)
    if op == "ln":
        return div(1.0, a) if a > 0 else 0.0
    if op == "sin":
        return trig(math.cos, a)
    if op == "cos":
        return canon(-trig(math.sin, a))
    if op == "tan":
        t = trig(math.tan, a)
        return canon(1.0 + mul(t, t))
    if op == "floor":
        return 0.0
    raise KeyError(op)


def d_binary(op, a, b):
    if op == "+":
        return [1.0, 1.0]
    if op == "-":
        return [1.0, -1.0]
    if op == "*":
        return [b, a]
    if op == "/":
        if b == 0:
            return [0.0, 0.0]
        # -(a / b) / b, not -a / b^2: b^2 underflows to 0 for tiny b.
        return [div(1.0, b), canon(-div(div(a, b), b))]
    if op == "min":
        return [1.0, 0.0] if a < b else [0.0, 1.0] if b < a else [0.5, 0.5]
    if op == "max":
        return [1.0, 0.0] if a > b else [0.0, 1.0] if b > a else [0.5, 0.5]
    if op == "atan2":
        r2 = canon(mul(a, a) + mul(b, b))
        if r2 == 0:
            return [0.0, 0.0]
        return [div(b, r2), canon(-div(a, r2))]
    if op == "pow":
        v = pow_(a, b)
        # 0 where the value is the constant 0 (a < 0, b not whole).
        fa = 0.0 if b == 0 or (a < 0 and not is_int(b)) else mul(b, pow_(a, b - 1.0))
        fb = mul(v, math.log(a)) if a > 0 else 0.0
        return [fa, fb]
    if op == "emax":
        if a > 0 and b > 0:
            e = hypot_(a, b)
            return [div(a, e), div(b, e)]
        return d_binary("max", a, b)
    raise KeyError(op)


def derivative(op, args):
    """Partials with respect to each argument (padded to three)."""
    args = [canon(v) for v in args]
    if op in UNARY:
        d = [d_unary(op, args[0])]
    elif op in BINARY:
        d = d_binary(op, *args)
    elif op == "select":
        d = [0.0, 1.0, 0.0] if args[0] < 0 else [0.0, 0.0, 1.0]
    else:
        raise KeyError(op)
    d = [canon(v) for v in d]
    return d + [0.0] * (3 - len(d))


def expr(op):
    if op in UNARY:
        return f"({op} x)"
    if op in BINARY:
        return f"({op} x y)"
    return "(select x y z)"


def arity(op):
    return 1 if op in UNARY else 2 if op in BINARY else 3


# ------------------------------------------------------------------ cases

# Inputs that probe every rule: zeros of both signs, units, kinks, poles,
# overflow thresholds, tiny and huge magnitudes, infinities.
SPECIAL = [0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 3.0, -3.0, 1.5, -1.5, 0.1, 1e-8, -1e-8,
           1e-300, -1e-300, 1e300, -1e300, math.pi / 2, -math.pi / 2, math.pi, 709.0, 710.0,
           -745.0, 100.0, 1e15, INF, -INF]
PAIRS = [0.0, -0.0, 1.0, -1.0, 0.5, -2.0, 3.0, 2.5, -1.5, 1e-300, 1e300, -1e300, INF, -INF]


def enc(v):
    if v == INF:
        return "inf"
    if v == -INF:
        return "-inf"
    return v


def point_cases(rng):
    out = []
    for op in list(UNARY) + list(BINARY) + ["select"]:
        n = arity(op)
        pts = []
        if n == 1:
            pts = [[a] for a in SPECIAL] + [[rng.uniform(-10, 10)] for _ in range(40)]
        elif n == 2:
            pts = [[a, b] for a in PAIRS for b in PAIRS]
            pts += [[rng.uniform(-10, 10), rng.uniform(-10, 10)] for _ in range(60)]
            # Exact ties away from zero.
            pts += [[v, v] for v in (1.0, -2.0, 3.5, 1e300)]
            # A negative base with an exponent just off a whole number: the
            # value is the constant 0, so is the slope (b - 1 rounds to -1).
            pts += [[-1.0, 2.0 ** -54], [-2.0, 3.0 + 2.0 ** -51]]
        else:
            pts = [[c, a, b] for c in (-1.0, 0.0, -0.0, 1.0, -INF, INF) for a in (2.0, -3.0) for b in (5.0, INF)]
        for p in pts:
            args = p + [0.0] * (3 - n)
            out.append({"op": op, "expr": expr(op), "at": [enc(v) for v in args],
                        "value": enc(value(op, p)), "d": [enc(v) for v in derivative(op, p)]})
    return out


# Boxes for range rules: each argument's interval.
BOXES1 = [(-1, 1), (0, 1), (-2, -1), (1, 2), (-3, 3), (0.5, 1.5), (-0.5, 0.5), (0, 0), (2, 2),
          (-1e300, 1e300), (0, 1e300), (-10, -1e-300), (1.4, 1.7), (-4, 4), (3.1, 3.2),
          (700, 720), (-0.1, 0.1), (2.5, 7.5)]
BOXES2 = [(-1, 1), (0, 2), (-3, -1), (1, 3), (-0.5, 0.5), (0, 0), (2, 2), (-1e300, 1e300), (0.5, 4)]


def samples(lo, hi, n, extra=()):
    """Endpoints, a uniform sweep, and special points inside [lo, hi]."""
    if lo == hi:
        return [lo]
    xs = {lo, hi}
    for k in range(1, n):
        xs.add(lo + (hi - lo) * k / n)
    for v in extra:
        if lo <= v <= hi:
            xs.add(v)
    # Integers (floor's jumps, pow's integer exponents) and multiples of pi/2
    # (tan's poles) when the box is not huge.
    if hi - lo <= 20:
        for k in range(math.ceil(lo), math.floor(hi) + 1):
            xs.add(float(k))
            xs.add(k - 1e-12)
            xs.add(k + 1e-12)
        for k in range(math.ceil(lo / (math.pi / 2)), math.floor(hi / (math.pi / 2)) + 1):
            xs.add(k * math.pi / 2)
    return sorted(v for v in xs if lo <= v <= hi)


KINKS = (0.0, 1e-300, -1e-300, 1.0, -1.0)


def range_cases():
    out = []
    for op in list(UNARY) + list(BINARY) + ["select"]:
        n = arity(op)
        if n == 1:
            boxes = [[b] for b in BOXES1]
        elif n == 2:
            boxes = [[a, b] for a in BOXES2 for b in BOXES2]
        else:
            boxes = [[c, (1, 2), (3, 5)] for c in ((-2, -1), (-1, 1), (0, 1), (0, 0))]
        for box in boxes:
            axes = [samples(lo, hi, 400 if n == 1 else 40, KINKS) for lo, hi in box]
            vals, ders = [], [[], [], []]
            for p in _product(axes):
                vals.append(value(op, p))
                for k, d in enumerate(derivative(op, p)[:n]):
                    ders[k].append(d)
            full = [list(b) for b in box] + [[0.0, 0.0]] * (3 - n)
            out.append({"op": op, "expr": expr(op), "box": [[enc(lo), enc(hi)] for lo, hi in full],
                        "hull": [enc(min(vals)), enc(max(vals))],
                        "dhull": [[enc(min(d)), enc(max(d))] if d else [0.0, 0.0] for d in ders]})
    return out


def _product(axes):
    if not axes:
        yield []
        return
    for v in axes[0]:
        for rest in _product(axes[1:]):
            yield [v] + rest


F32_MAX = 3.4028234663852886e38


def f32(v):
    """v rounded to the nearest f32 (infinite beyond its range)."""
    try:
        return struct.unpack("f", struct.pack("f", v))[0]
    except OverflowError:
        return math.copysign(INF, v)


# Small expression trees, interpreted independently of the kernel. Each
# node carries a value and five derivatives: x, y, z, t, p. In particular,
# warp substitutes *all* inputs simultaneously and composes their Jacobian.
VARS = ("x", "y", "z", "t", "p")


def parse(src):
    tokens = iter(re.findall(r"\(|\)|[^\s()]+", src))

    def read(token):
        if token == "(":
            items = []
            while (token := next(tokens)) != ")":
                items.append(read(token))
            return items
        try:
            return float(token)
        except ValueError:
            return token

    tree = read(next(tokens))
    assert next(tokens, None) is None
    return tree


def jet(tree, env):
    if isinstance(tree, (int, float)):
        return float(tree), [0.0] * 5
    if isinstance(tree, str):
        return env[tree]
    op, *children = tree
    if op == "warp":
        inner = dict(env)
        for name, child in zip(VARS[:4], children[1:]):
            inner[name] = jet(child, env)
        return jet(children[0], inner)
    args = [jet(c, env) for c in children]
    vals = [v for v, _ in args]
    if op in ("<", "<=", ">", ">=", "=="):
        a, b = vals
        return float({"<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b, "==": a == b}[op]), [0.0] * 5
    partials = derivative(op, vals)[:len(args)]
    d = [canon(sum(mul(w, a[1][k]) for w, a in zip(partials, args))) for k in range(5)]
    return value(op, vals), d


def evaluate(tree, at, context=(0.0, 0.0)):
    inputs = list(at) + list(context)
    return jet(tree, {v: (inputs[k], [float(j == k) for j in range(5)]) for k, v in enumerate(VARS)})


COMPOSED = [
    "(sqrt (+ (+ (square x) (square y)) (square z)))",
    "(- (sqrt (+ (square (- x 1)) (square (+ y 2)))) z)",
    "(min (min x y) z)", "(max x (max y z))",
    "(min (max x y) (max (neg x) z))",
    "(max (min x y) (min (neg y) z))",
    "(abs (- (square x) (square y)))",
    "(square (+ x (* y z)))",
    "(/ (+ x y) (+ 1 (square z)))",
    "(sqrt (max (- (* x y) z) 0))",
    "(select (- x y) (+ x z) (* y z))",
    "(select (neg (< x y)) (square x) (square y))",
    "(select (neg (<= x y)) (+ x z) (- y z))",
    "(select (neg (== x y)) (* x z) (+ y z))",
    "(sin (+ x (* y z)))", "(cos (- (* x y) z))",
    "(tan (* 0.125 (+ x y)))", "(atan2 (- x (* y z)) 1)",
    "(exp (* 0.25 (- x y)))", "(ln (+ 1 (square (+ x y))))",
    "(pow (+ 1 (square x)) (+ y (* 0.25 z)))",
    "(pow (- x y) 3)", "(pow (+ 1 (abs y)) (sin z))",
    "(atan2 (- x y) (+ 1 (square z)))",
    "(sin (atan2 y x))", "(emax (- x y) (+ y z))",
    "(+ (floor (* 0.5 x)) (* y z))",
    "(select z (ln (+ 1 (square x))) (sqrt (+ 1 (square y))))",
    "(warp (+ (square x) (* y z)) (+ x y) (- y z) (* 0.5 z))",
    "(warp (atan2 y x) (- x z) (+ x y) z)",
    "(warp (warp (+ x (* y z)) y z x) (+ x 1) (* 2 y) (- z 1))",
    "(warp (sqrt (+ (square x) (square y))) (- x y) (+ x y) z)",
]
CONTEXTUAL = [
    "(+ (* p x) (* t y))",
    "(sin (+ x (* p t)))",
    "(pow (+ 1 (square (- x p))) (+ 1 (* 0.25 t)))",
    "(select (- t p) (+ x (* t z)) (- y (* p z)))",
    "(warp (+ x (* t y)) (- x p) (+ y t) z)",
    "(warp (+ x (* t y)) (- x p) (+ y t) z (+ t (* p x)))",
    "(warp (warp (+ x t) y z x (+ t z)) (+ x p) y z (* t x))",
    "(+ (square p) (- t z))",
]


def expression_cases():
    points, ranges = [], []
    ats = [(0, 0, 0), (1, 1, 1), (-1, -1, -1), (0.25, -0.5, 0.75),
           (-0.75, 0.5, -0.25), (2, -1, 0.5)]
    boxes = [[[-1, 1]] * 3, [[0.25, 0.75], [-0.75, -0.25], [0.25, 0.5]]]
    for category, exprs, contexts in [("composed", COMPOSED, [(0, 0)]),
                                      ("context", CONTEXTUAL, [(0.75, -1.25), (-1.5, 2)])]:
        for i, src in enumerate(exprs):
            tree = parse(src)
            for context in contexts:
                base = {"op": f"{category}/{i:02}", "category": category, "expr": src, "context": context, "ulps": 64}
                for at in ats:
                    v, d = evaluate(tree, at, context)
                    points.append(dict(base, at=at, value=enc(v), d=list(map(enc, d))))
                for box in boxes:
                    axes = [samples(lo, hi, 8, (0,)) for lo, hi in box]
                    results = [evaluate(tree, at, context) for at in _product(axes)]
                    ranges.append(dict(base, box=box,
                                       hull=[enc(min(v for v, _ in results)), enc(max(v for v, _ in results))],
                                       dhull=[[enc(min(d[k] for _, d in results)), enc(max(d[k] for _, d in results))] for k in range(3)]))
            # Over a range of t and of p too (a time span, a tolerance): the
            # range rules must hold for every value of them. Spans crossing
            # 0 and the select's switch (t = p), and one narrowed to a value.
            if category == "context":
                for spans in CONTEXT_SPANS:
                    base = {"op": f"{category}/{i:02}/spans", "category": category, "expr": src, "context_spans": spans, "ulps": 64}
                    for box in boxes:
                        axes = [samples(lo, hi, 4, (0,)) for lo, hi in list(box) + list(spans)]
                        results = [evaluate(tree, at[:3], at[3:]) for at in _product(axes)]
                        ranges.append(dict(base, box=box,
                                           hull=[enc(min(v for v, _ in results)), enc(max(v for v, _ in results))],
                                           dhull=[[enc(min(d[k] for _, d in results)), enc(max(d[k] for _, d in results))] for k in range(3)]))
    return points, ranges


CONTEXT_SPANS = [[(0.5, 1.0), (-1.5, -1.0)], [(-0.25, 0.75), (-0.5, 1.5)], [(-2.0, -1.0), (0.375, 0.375)]]


def extreme_cases():
    """Deliberate boundary and pole probes at f32's limits (a float's
    largest and smallest normal values, division and ln at zero, either
    side of tan's pole), with the inputs rounded to f32 first.
    """
    inputs = {
        "neg": [(F32_MAX,), (-F32_MAX,)], "abs": [(-F32_MAX,)],
        "sqrt": [(F32_MAX,), (2**-126,)],
        "ln": [(0,), (-1,), (2**-126,), (F32_MAX,), (1-2**-23,), (1+2**-23,)],
        "exp": [(88,), (89,), (-87,)],
        "tan": [(f32(math.pi/2) + d,) for d in (-2**-22, -2**-23, 0, 2**-23, 2**-22)],
        "/": [(1, 0), (-1, 0), (F32_MAX, 0.5), (-F32_MAX, 0.5), (2**-126, 2**-126), (1, 2**-126), (1e-12, -1)],
        "pow": [(0, -0.5), (-1, 0.5), (1e20, 2.5), (2**-126, -0.5)],
        "atan2": [(F32_MAX, F32_MAX), (F32_MAX, -F32_MAX), (2**-126, 2**-126)],
        "emax": [(F32_MAX/4, F32_MAX/4)],
        "+": [(F32_MAX/4, F32_MAX/4)], "*": [(F32_MAX, 0.5)],
        "sin": [(1e-12,), (-1e-12,)],
    }
    points = []
    for op, pts in inputs.items():
        for i, args in enumerate(pts):
            args = list(map(f32, args))
            at = args + [0.0] * (3-len(args))
            v = value(op, args)
            base = {"op": f"extreme/{op}/{i}", "category": "extreme", "expr": expr(op), "at": at}
            points.append(dict(base, value=enc(v), d=list(map(enc, derivative(op, args))), ulps=4 if op in ("sqrt", "ln", "exp", "tan", "pow", "atan2", "emax", "sin") else 0))
    return points


# There is no grid literal in the IR. Bake this bounded polynomial once.
# The layout is a fixture, not a reference implementation of the bounds
# service: the harness asserts it, including every independently computed
# sample, before testing the interpolant. A changed bake layout must update
# this fixture; it must never silently change what the reference tests.
GRID_ORIGIN = [-4.098632812499973] * 3
GRID_DIMS = [10] * 3
GRID_H = 1.0
PRELUDE = """(param p 0)
(def sampled (grid-bake (max (max (- (abs x) 6) (max (- (abs y) 6) (- (abs z) 6)))
  (+ (+ (+ (square x) (square y)) (square z))
     (- (* 0.125 (* x (* y z))) 1))) (res 1)))
"""


def polynomial(x, y, z):
    return x*x + y*y + z*z + 0.125*x*y*z - 1


def grid_samples():
    return [f32(polynomial(*(GRID_ORIGIN[a] + GRID_H * c[a] for a in range(3))))
            for k in range(10) for j in range(10) for i in range(10) for c in [(i, j, k)]]


def grid_component(at, mask, cell=None):
    """Tensor-product hat functions; a differentiated axis replaces its
    two weights by -1/h and +1/h. This directly gives all eight components.
    At faces floor selects the upper cell; the outermost face uses the last.
    """
    q = [max(0, min(GRID_DIMS[a]-1, (at[a]-GRID_ORIGIN[a])/GRID_H)) for a in range(3)]
    cell = cell or [min(math.floor(q[a]), GRID_DIMS[a]-2) for a in range(3)]
    weights = [[-1/GRID_H, 1/GRID_H] if mask & (1 << a) else [1-(q[a]-cell[a]), q[a]-cell[a]] for a in range(3)]
    terms = []
    for bits in itertools.product((0, 1), repeat=3):
        p = [GRID_ORIGIN[a] + GRID_H * (cell[a]+bits[a]) for a in range(3)]
        terms.append(f32(polynomial(*p)) * math.prod(weights[a][bits[a]] for a in range(3)))
    return math.fsum(terms)


def grid_jet(at, mask, cell=None):
    return grid_component(at, mask, cell), [0.0 if mask & (1 << a) else grid_component(at, mask | (1 << a), cell) for a in range(3)]


def baked_jet(at, mask):
    # grid-bake exposes the continuation, written in plain ops: clamp each
    # coordinate, interpolate, then add distance to the box. At a box face
    # min/max's tie rule gives the clamp derivative 1/2; the distance term
    # and its derivatives evaluate to zero there under the total rules.
    factors = [0.5 if at[a] in (GRID_ORIGIN[a], GRID_ORIGIN[a] + 9*GRID_H) else 1.0 for a in range(3)]
    def component(m):
        return grid_component(at, m) * math.prod(factors[a] for a in range(3) if m & (1 << a))
    return component(mask), [0.0 if mask & (1 << a) else component(mask | (1 << a)) for a in range(3)]


def sampled_cases():
    points, ranges = [], []
    def world(q):
        return [GRID_ORIGIN[a] + GRID_H * q[a] for a in range(3)]
    # Interiors, each axis's face (and just below/above), first/last cells,
    # all three faces meeting, and the grid's final vertex.
    qs = [(0, 0, 0), (0.25, 0.5, 0.75), (3.25, 4.5, 5.75), (4, 4, 4),
          (8.25, 8.5, 8.75), (9, 9, 9)]
    for a in range(3):
        for offset in (-2**-10, 0, 2**-10):
            q = [3.25, 3.5, 3.75]
            q[a] = 4 + offset
            qs.append(q)
    boxes = [([3.125, 3.25, 3.5], [3.75, 3.875, 3.75]),
             ([3.5]*3, [4.5]*3), ([8.25]*3, [9]*3), ([4]*3, [4]*3)]
    for mask in range(8):
        src = "sampled"
        for a, axis in enumerate(VARS[:3]):
            if mask & (1 << a):
                src = f"(diff {src} {axis})"
        base = {"op": f"sampled/{mask}", "category": "sampled", "expr": src, "ulps": 64}
        for q in qs:
            at = world(q)
            v, d = baked_jet(at, mask)
            points.append(dict(base, at=at, value=v, d=d))
        for qlo, qhi in boxes:
            lo, hi = world(qlo), world(qhi)
            vals, ders = [], [[], [], []]
            # Every intersection of a cell and the box is multilinear.
            # Its extrema occur at its corners. Include both one-sided
            # limits at an interior face, but a box starting on a face
            # only includes the upper cell (the spec's ownership rule).
            for cell in itertools.product(*(range(min(8, math.floor(qlo[a])), min(8, math.floor(qhi[a]))+1) for a in range(3))):
                corners = [[max(lo[a], world(cell)[a]), min(hi[a], world([c+1 for c in cell])[a])] for a in range(3)]
                for at in itertools.product(*corners):
                    v, d = grid_jet(at, mask, cell)
                    vals.append(v)
                    for a in range(3):
                        ders[a].append(d[a])
            for at in itertools.product(*zip(lo, hi)):
                v, d = baked_jet(at, mask)
                vals.append(v)
                for a in range(3):
                    ders[a].append(d[a])
            ranges.append(dict(base, box=list(map(list, zip(lo, hi))), hull=[min(vals), max(vals)],
                               dhull=[[min(d), max(d)] for d in ders]))
    return points, ranges


def self_checks():
    # Hand calculations pin the chain rule and simultaneous substitution.
    v, d = evaluate(parse("(sqrt (+ (square x) (square y)))"), [3, 4, 0])
    assert v == 5.0 and all(math.isclose(a, b, rel_tol=1e-15) for a, b in zip(d, [0.6, 0.8, 0.0, 0.0, 0.0]))
    assert evaluate(parse("(warp (+ x (* t y)) y x z (+ t p))"), [2, 3, 4], (5, 7)) == (27.0, [12.0, 1.0, 0.0, 2.0, 2.0])


def main():
    self_checks()
    rng = random.Random(20261007)
    cases = {
        "note": "Generated by conformance/reference.py; do not edit by hand. Infinities are the strings inf and -inf.",
        "points": point_cases(rng),
        "ranges": range_cases(),
    }
    for generate in (expression_cases, sampled_cases):
        points, ranges = generate()
        cases["points"].extend(points)
        cases["ranges"].extend(ranges)
    cases["points"].extend(extreme_cases())
    cases["prelude"] = PRELUDE
    cases["grid"] = {"origin": GRID_ORIGIN, "dims": GRID_DIMS, "h": GRID_H, "values": grid_samples()}
    print(json.dumps(cases, indent=None, separators=(",", ":")))


if __name__ == "__main__":
    main()
