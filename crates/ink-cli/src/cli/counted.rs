//! The CPU instructions a piece of work retires, counted by the processor
//! (Linux perf events, user space only; elsewhere none is counted): what the view budget's ratios
//! and Tape::work are calibrated on, by hand
//! (scripts/view-budget-calibrate.sh runs the measurements, ignored tests,
//! which print what the tests then hold the kernel's counts to). Counted,
//! not timed, they do not move with the machine's load or clock, only with
//! its instruction set and the build.

/// Run `f` on a pool of one thread: its parallel work runs on that thread,
/// so `instructions` there counts all of it.
pub fn on_one_thread<R: Send>(f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new().num_threads(1).build().expect("a thread").install(f)
}

/// What a piece of work cost the thread that ran it: the instructions it
/// retired, and the cycles and CPU seconds they took (which move with the
/// machine, and its load).
#[derive(Clone, Copy, Debug)]
pub struct Cost {
    pub instructions: u64,
    pub cycles: u64,
    pub seconds: f64,
}

/// The calling thread's CPU time.
fn cpu() -> f64 {
    let mut t = libc::timespec { tv_sec: 0, tv_nsec: 0 };
    // SAFETY: writes the calling thread's CPU clock into `t`.
    unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut t) };
    t.tv_sec as f64 + t.tv_nsec as f64 * 1e-9
}

/// `f`'s result and what the calling thread spent running it (all of its
/// work on a pool of one: `on_one_thread`). None where the system does not
/// let a process count its own (perf_event_paranoid above 2, or no such
/// counter).
pub fn instructions<R>(f: impl FnOnce() -> R) -> Option<(R, Cost)> {
    let first = open(1)?;
    let Some(second) = open(0) else {
        // SAFETY: closes the descriptor `open` returned.
        unsafe { libc::close(first) };
        return None;
    };
    let fds = [first, second];
    let t = cpu();
    let r = f();
    let seconds = cpu() - t;
    let n = fds.map(|fd| {
        let mut n = 0u64;
        // SAFETY: reads 8 bytes into `n` from a counter this thread opened.
        let got = unsafe { libc::read(fd, (&raw mut n).cast(), 8) };
        // SAFETY: closes the descriptor `open` returned.
        unsafe { libc::close(fd) };
        (got == 8).then_some(n)
    });
    Some((r, Cost { instructions: n[0]?, cycles: n[1]?, seconds }))
}

/// A counter of the calling thread's hardware event `config` in user space
/// (0 cycles, 1 instructions retired), counting from now.
#[cfg(target_os = "linux")]
fn open(config: u64) -> Option<i32> {
    // struct perf_event_attr as of its first version (64 bytes): type,
    // size, config, then (after sample_period, sample_type, read_format)
    // the flags: exclude_kernel (bit 5) and exclude_hv (bit 6).
    let mut attr = [0u64; 8];
    attr[0] = 64u64 << 32; // PERF_TYPE_HARDWARE (0), size
    attr[1] = config;
    attr[5] = (1 << 5) | (1 << 6);
    // SAFETY: perf_event_open with a valid attr for this thread (pid 0,
    // any cpu, no group, no flags).
    let fd = unsafe { libc::syscall(libc::SYS_perf_event_open, attr.as_ptr(), 0, -1, -1, 0) };
    (fd >= 0).then_some(fd as i32)
}

/// No counter outside Linux: perf events are Linux's.
#[cfg(not(target_os = "linux"))]
fn open(_config: u64) -> Option<i32> {
    None
}

mod tests {
    /// A loop's instructions are counted, and twice the loop about twice.
    #[test]
    #[ignore = "a calibration's, run by hand"]
    fn counts_instructions() {
        let spin = |n: u64| (0..n).fold(0u64, |a, i| std::hint::black_box(a.wrapping_mul(31).wrapping_add(i)));
        let Some((_, a)) = super::instructions(|| spin(1_000_000)) else {
            eprintln!("this system does not let a process count its instructions");
            return;
        };
        let (_, b) = super::instructions(|| spin(2_000_000)).unwrap();
        let (a, b) = (a.instructions, b.instructions);
        assert!(a > 1_000_000 && (b as f64 / a as f64 - 2.0).abs() < 0.1, "{a} then {b}");
    }
}

/// How far, as a factor either way, Tape::work and Tape::simplify_work
/// may be from what evaluating and simplifying a tape retired when last
/// counted (`ops::COUNTED`).
const OP_TOLERANCE: f64 = 1.3;

mod ops {
    use ink::eval::{Ctx, eval_grads, eval_interval, eval_points};
    use ink::interval::I;
    use ink::tape::{BOTH, Pass, Tape};

    /// Tapes of 64 terms of each kind of instruction, summed.
    fn cases() -> Vec<(&'static str, Tape, Ctx)> {
        let terms = |f: &dyn Fn(usize) -> String| (0..64).map(f).collect::<Vec<_>>().join(" ");
        let c = |i: usize| 1.0 + i as f64 / 97.0;
        let cases: Vec<(&str, String)> = vec![
            ("arithmetic", terms(&|i| format!("(* x {})", c(i)))),
            ("mul", terms(&|i| format!("(* (* x {}) y)", c(i)))),
            ("div", terms(&|i| format!("(/ (* x {}) (+ 2 (abs y)))", c(i)))),
            ("sqrt", terms(&|i| format!("(sqrt (* x {}))", c(i)))),
            ("sin", terms(&|i| format!("(sin (* x {}))", c(i)))),
            ("cos", terms(&|i| format!("(cos (* x {}))", c(i)))),
            ("tan", terms(&|i| format!("(tan (* x {}))", c(i)))),
            ("exp", terms(&|i| format!("(exp (* x {}))", c(i)))),
            ("ln", terms(&|i| format!("(ln (* x {}))", c(i)))),
            ("atan2", terms(&|i| format!("(atan2 y (* x {}))", c(i)))),
            ("pow", terms(&|i| format!("(pow (abs (* x {})) 1.7)", c(i)))),
            ("floor", terms(&|i| format!("(floor (* x {}))", c(i)))),
            ("min", terms(&|i| format!("(min (* x {}) (- y {}))", c(i), i))),
            ("select", terms(&|i| format!("(select (* x {}) y z)", c(i)))),
            ("rad", terms(&|i| format!("(deg->rad (* x {}))", c(i)))),
            ("emax", terms(&|i| format!("(emax (* x {}) y)", c(i)))),
            ("grid", terms(&|i| format!("(move g {} 0 0)", c(i)))),
        ];
        cases
            .into_iter()
            .map(|(name, src)| {
                let src = format!("(import \"std\")\n(def g (grid-bake (sphere 5) (res 0.5)))\n(body b (shape (+ {src})))");
                let mut doc = ink::load_str(&src, "ops", std::path::Path::new(".")).unwrap();
                let ctx = doc.ctx(&[], 0.0).unwrap();
                let it = doc.items("b").unwrap().remove(0);
                (name, Tape::new(&doc.g, it.shape), ctx)
            })
            .collect()
    }

    fn points() -> Vec<[f64; 3]> {
        (0..4096).map(|k| [0.1 + (k % 16) as f64 * 0.37, 0.2 + (k / 16 % 16) as f64 * 0.29, 0.3 + (k / 256) as f64 * 0.31]).collect()
    }

    fn boxes() -> Vec<[I; 3]> {
        points().iter().take(512).map(|p| p.map(|v| I::new(v, v + 0.05))).collect()
    }

    /// Each case's tape simplified over each box.
    fn simplified(tape: &Tape, ctx: &Ctx) -> Vec<Tape> {
        let simplify = |b: &[I; 3]| {
            let mut ch = vec![BOTH; tape.len()];
            eval_interval(tape, ctx, *b, Some(&mut ch));
            tape.simplify(&ch)
        };
        boxes().iter().map(simplify).collect()
    }

    /// What each case's tape retired, as `measure_op_work` printed it: a
    /// point's evaluation, a gradient's, a box's, and simplifying it over
    /// a box beyond that box's pass, on one thread of a Zen 5 (AVX-512),
    /// the test build.
    const COUNTED: &[(&str, [f64; 4])] = &[
        ("arithmetic", [189.4, 14922.7, 32056.5, 156636.7]),
        ("mul", [274.0, 22428.9, 55568.9, 210746.0]),
        ("div", [308.4, 23257.1, 48341.9, 214321.9]),
        ("sqrt", [249.5, 19409.7, 42123.4, 205973.6]),
        ("sin", [6076.6, 28872.4, 61712.7, 205980.3]),
        ("cos", [6185.9, 28936.4, 64903.7, 206485.6]),
        ("tan", [7955.5, 25866.8, 63339.6, 205986.4]),
        ("exp", [4971.4, 22755.7, 55220.1, 205972.9]),
        ("ln", [5336.4, 23631.7, 56691.1, 205984.1]),
        ("atan2", [7910.7, 29926.6, 76190.2, 213150.5]),
        ("pow", [17650.2, 64514.1, 94874.2, 261812.3]),
        ("floor", [232.5, 19089.7, 49675.2, 207903.7]),
        ("min", [384.4, 30007.8, 58717.6, 171139.0]),
        ("select", [264.5, 17725.3, 38240.8, 72611.9]),
        ("rad", [1302.5, 19153.8, 54275.8, 205984.5]),
        ("emax", [4635.9, 26012.9, 65175.2, 210526.0]),
        ("grid", [101204.0, 1115657.5, 923281.9, 795908.9]),
    ];

    /// What a pass of a tape retires (Tape::work: at points, their
    /// gradients and over boxes) and simplifying it (simplify_work), over
    /// 64 terms of each kind of instruction, is within OP_TOLERANCE of what
    /// was last counted (COUNTED). Counted by the kernel, so it checks the
    /// same on any machine.
    #[test]
    fn a_tapes_work_is_what_it_retires() {
        let (mut off, cases) = (Vec::new(), cases());
        assert_eq!(cases.len(), COUNTED.len(), "a case for each counted");
        for ((name, tape, ctx), &(counted_name, counted)) in cases.into_iter().zip(COUNTED) {
            assert_eq!(name, counted_name);
            let simplify = simplified(&tape, &ctx).iter().map(|t| tape.simplify_work(t)).sum::<usize>() as f64 / boxes().len() as f64;
            let said = [tape.work(Pass::Point) as f64, tape.work(Pass::Gradient) as f64, tape.work(Pass::Box) as f64, simplify];
            let ratio: [f64; 4] = std::array::from_fn(|k| counted[k] / said[k]);
            off.extend(
                ratio
                    .iter()
                    .enumerate()
                    .filter(|r| !(r.1 < &super::OP_TOLERANCE && *r.1 > 1.0 / super::OP_TOLERANCE))
                    .map(|(k, r)| (name, ["point", "gradient", "box", "simplify"][k], *r)),
            );
        }
        assert!(off.is_empty(), "counted over said: {off:?}");
    }

    /// What each case retires, counted on one thread, as COUNTED's rows.
    #[test]
    #[ignore = "a calibration's, run by hand"]
    fn measure_op_work() {
        let (pts, boxes) = (points(), boxes());
        for (name, tape, ctx) in cases() {
            let counted = super::on_one_thread(|| {
                let (_, p) = super::instructions(|| eval_points(&tape, &ctx, &pts, &mut Vec::new()))?;
                let (_, g) = super::instructions(|| eval_grads(&tape, &ctx, &pts))?;
                let (_, b) = super::instructions(|| boxes.iter().map(|b| eval_interval(&tape, &ctx, *b, None).lo).sum::<f64>())?;
                let (_, s) = super::instructions(|| simplified(&tape, &ctx))?;
                let n = |c: super::Cost, k: usize| c.instructions as f64 / k as f64;
                Some([n(p, pts.len()), n(g, pts.len()), n(b, boxes.len()), n(s, boxes.len()) - n(b, boxes.len())])
            });
            let c = counted.expect("this system does not let a process count its instructions");
            println!("        (\"{name}\", [{:.1}, {:.1}, {:.1}, {:.1}]),", c[0], c[1], c[2], c[3]);
        }
    }
}
