use crate::{Dynamics, model};
use ink_lang::Doc;
use ink_lang::view::ModelView;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering::SeqCst};
use std::time::Duration;

thread_local! {
    /// Run where a simulation's model is built (Run::sim), once: what a
    /// thread waiting in rayon there may run meanwhile.
    pub(super) static BUILDING: std::cell::Cell<Option<Box<dyn FnOnce()>>> = const { std::cell::Cell::new(None) };
}

static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&Dynamics]);

const PENDULUM: &str = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (def swing (simulate 60))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
        (assembly root (instance p (pendulum swing)))"#;

fn load(name: &str) -> Result<(), String> {
    Doc::load_str_in(&REGISTRY, PENDULUM, name, Path::new("."), &[]).map(drop).map_err(|e| e.message)
}

/// A simulation's model is built (its masses measured) on rayon's
/// workers, and a worker waiting there runs whatever else is queued:
/// loading many files side by side, other files' loads (a CoreXY
/// printer's ringing ran 167 of them while its model was built). Here, forced in a pool of two: where the model
/// is built, a join's second half goes to the other worker, and the
/// first, waiting for it, runs another simulated file's load queued on
/// its own thread. Built under the lock on the simulations made (or
/// inside the run's initialiser), that load waits for the lock its own
/// thread holds, and every thread waits; built
/// outside, both load.
#[test]
fn a_load_run_while_a_simulation_is_built_loads() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        let inner = Arc::new(std::sync::Mutex::new(None));
        let got = inner.clone();
        let outer = pool.install(move || {
            BUILDING.set(Some(Box::new(move || {
                let begun = AtomicBool::new(false);
                rayon::join(
                    || {
                        // Queued on this thread, after the other half: a
                        // thief takes the oldest, the half.
                        rayon::spawn(move || *got.lock().unwrap() = Some(load("inner.ink")));
                        while !begun.load(SeqCst) {
                            std::thread::yield_now();
                        }
                    },
                    || {
                        begun.store(true, SeqCst);
                        std::thread::sleep(Duration::from_millis(200));
                    },
                );
            })));
            load("outer.ink")
        });
        // The inner load ran inside the outer's join, or runs on now.
        let inner = loop {
            if let Some(r) = inner.lock().unwrap().take() {
                break r;
            }
            std::thread::sleep(Duration::from_millis(10));
        };
        _ = tx.send((outer, inner));
    });
    let (outer, inner) = rx.recv_timeout(Duration::from_secs(60)).expect("the loads deadlocked: a simulation's load run while another's model is built waits for itself");
    assert_eq!((outer, inner), (Ok(()), Ok(())));
}

/// `ink check` asks each solve on demand whether it refuses
/// (ModelView::mass_checks): an effort says what refuses it at every
/// moment without stepping its motion, which its claim steps.
#[test]
fn check_does_not_step_an_effort_motion() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
            (body arm (shape (move (box 200 20 10) 100 0 0)) (attr density 2.7g/cm3))
            (assembly (joint (angle 0)) (instance arm arm (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))
            (assembly root (instance j (joint (* 45 (sin (* 180deg/s t))))))"#;
    let doc = Doc::load_str_in(&REGISTRY, src, "arm.ink", Path::new("."), &[]).unwrap();
    let (ctx, _) = doc.ctx_at(&[], 0.0).unwrap();
    let mut v = ModelView::new(doc, ctx);
    let before = model::CHOSEN.get();
    assert!(v.mass_checks().is_empty());
    assert_eq!(model::CHOSEN.get(), before, "ink check stepped the motion");
}
