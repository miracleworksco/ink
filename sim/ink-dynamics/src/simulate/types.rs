use crate::model::Model;

/// How far a simulation was stepped over its span (Simulation::stepped):
/// to its end, or where it stopped short, when and why; its steps,
/// rejections and states evaluated.
pub struct Stepped {
    pub end: f64,
    pub steps: usize,
    pub rejected: usize,
    pub evaluations: usize,
    pub stopped: Option<(f64, String)>,
}

/// What a simulation's springs did (Simulation::springs): each released
/// joint carrying inertia's following error, the joints held where the
/// springs leave them, each spring line's least and greatest tension (N,
/// at a knot) with when, and the moment the written motion comes to rest
/// for good (None: it moves at the end).
pub struct Springs {
    pub joints: Vec<Following>,
    pub held: Vec<String>,
    pub spans: Vec<Span>,
    pub still: Option<f64>,
}

/// A released joint's following error over a simulation (mm a slide, deg a
/// turn): its largest and when, its value at the end, and its ringing once
/// the written motion is at rest.
pub struct Following {
    pub path: String,
    pub turn: bool,
    pub peak: f64,
    pub at: f64,
    pub end: f64,
    pub ringing: Option<Ringing>,
    /// Where it states a stiffness: its spring's largest force (N a slide,
    /// N mm a turn), what an actuator so held supplies.
    pub effort: Option<f64>,
}

/// A following error's ringing from moment `from`: half its swing there
/// (peak to peak), and the mode carrying most of it at the end: its
/// frequency (Hz) and its share of the ringing there (None where the
/// simulation stops short).
pub struct Ringing {
    pub from: f64,
    pub amplitude: f64,
    pub frequency: Option<f64>,
    pub share: Option<f64>,
}

/// A spring line's least and greatest tension over a simulation (t, N).
pub struct Span {
    pub name: String,
    pub site: String,
    pub least: (f64, f64),
    pub most: (f64, f64),
}

/// Why a state's accelerations are not found.
pub enum Fail {
    /// The state is refused (its forces not known, a dead point, a command
    /// not finite): why.
    Refused(String),
    /// No accelerations agree with the direction the couplings' losses turn
    /// with (a self-locking coupling driven back, `locked`, or Newton's
    /// method not converging): at rest, friction holds it still.
    Unbalanced { locked: bool, why: String },
}

/// What a simulation steps: the simulated coordinates (Model::simulated's
/// order), then the deviations of the released joints carrying inertia
/// along themselves (`w`, by index among Model::released), each a position
/// and a speed; those carrying none (`y`) held where the springs leave them.
pub struct Layout {
    pub ns: usize,
    pub w: Vec<usize>,
    pub y: Vec<usize>,
    pub nr: usize,
    /// Each stepped component's radius of gyration, how far the parts it
    /// moves move by a unit of it (a turn's, mm a radian: the square root
    /// of its inertia over their mass; 1 a slide or a coordinate), and its
    /// own frequency (a simulated coordinate's from how its acceleration
    /// changes with it at the start, a deviation's its springs' stiffness
    /// along it over its inertia; at least once over the span): what its
    /// tolerance reads (Simulation::scales).
    pub gyration: Vec<f64>,
    pub own: Vec<f64>,
    /// The springs' fastest mode at the start (Hz): the periods of it
    /// stepped widen a following error's band (Simulation::band).
    pub fastest: f64,
}

impl Layout {
    pub(super) fn n(&self) -> usize {
        self.ns + self.w.len()
    }

    /// Every released joint stepped (no massless one known).
    pub(super) fn of(model: &Model) -> Layout {
        let (ns, nr) = (model.simulated.len(), model.released().len());
        Layout { ns, w: (0..nr).collect(), y: Vec::new(), nr, gyration: vec![1.0; ns + nr], own: vec![1.0; ns + nr], fastest: 0.0 }
    }
}
