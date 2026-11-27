//! The shape package: two metrics of a shape a claim can bound,
//! `(min-thickness part > t)` and `(cavities part < 1)`. It has no forms.

use ink_core::graph::NodeId;
use ink_lang::claims::Held;
use ink_lang::sexp::Sexp;
use ink_lang::view::ModelView;
use ink_lang::{Excluded, InkError, Measured, Package, ReadError, Reader};
use std::sync::Arc;

pub struct Shape;

impl Package for Shape {
    fn name(&self) -> &'static str {
        "shape"
    }

    fn forms(&self) -> &'static [&'static str] {
        &[]
    }

    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[
            ("min-thickness", "the thinnest wall of a shape (mm), by rays from surface samples: (min-thickness part > t); of (air-around part m), the narrowest gap"),
            ("cavities", "how many enclosed voids a shape has (air no path leads out of) with air deeper than the accuracy, by proof: (cavities part < 1)"),
        ]
    }

    fn metric_quantity(&self, metric: &str) -> Option<ink_lang::units::Quantity> {
        (metric == "min-thickness").then_some(ink_lang::units::Quantity::Length)
    }

    fn shape_metrics(&self) -> &'static [&'static str] {
        &["min-thickness", "cavities"]
    }

    fn exclusions(&self) -> &'static [(&'static str, &'static str)] {
        &[("undecided-cells", "voids a count's finest cells could not rule out")]
    }

    fn shape_accuracy(&self, metric: &str) -> bool {
        metric == "cavities"
    }

    /// A prism's walls are its profile's in the plane, and it has no
    /// cavity: both are the same asked there.
    fn plane_metric(&self, _metric: &str) -> bool {
        true
    }

    /// The thinnest wall, by rays every way, and the voids are the shape's
    /// own: the same however it is turned or moved.
    fn rigid_metric(&self, _metric: &str) -> bool {
        true
    }

    fn shape_metric(&self, v: &mut ModelView, metric: &str, f: NodeId, res: Option<f64>, accuracy: f64) -> Result<Measured, InkError> {
        shape_metric(v, metric, f, res, accuracy)
    }

    fn shape_res(&self, v: &mut ModelView, metric: &str, f: NodeId, res: Option<f64>) -> Option<f64> {
        (metric == "min-thickness").then(|| v.thickness_res(f, res)).flatten()
    }

    fn read(&self, r: &mut Reader, form: &Sexp, head: &str, _args: &[Sexp]) -> Result<Arc<dyn ink_lang::Form>, InkError> {
        Err(r.err(form, format!("'{head}' is no form of shape")))
    }
}

/// The metrics of a shape (Package::shape_metrics): `min-thickness`,
/// an estimate, the thinnest wall by rays from surface samples spaced
/// `res` (ModelView::thickness: 128 along its box), where it is (`at`);
/// `cavities`, the enclosed voids with air deeper than `accuracy`, by
/// proof (ModelView::voids): those it proved, and, where its cells leave
/// more possible, the most there may be (`undecided-cells`, excluded);
/// decided, how it held (the thin air it set aside). An empty shape has no
/// wall and no void.
fn shape_metric(v: &mut ModelView, metric: &str, f: u32, res: Option<f64>, accuracy: f64) -> Result<Measured, InkError> {
    let measured = |value, at, res| Measured {
        value,
        at,
        res,
        finest: None,
        warnings: Vec::new(),
        excluded: Vec::new(),
        unseen: None,
        unseen_finer: true,
        resolves: 0.0,
        unbounded: None,
        extent: None,
        read_error: None,
        proven: false,
        held: None,
        band: None,
    };
    Ok(match metric {
        "min-thickness" => match v.thickness(f, res) {
            None => measured(Ok(f64::INFINITY), None, None),
            Some(Err((h, why))) => measured(Err(why.into()), None, Some(h)),
            // Each ray runs through material to where it ends: a wall no
            // ray crosses, or one crossed slant, only reads thicker.
            Some(Ok((at, t, h))) => Measured {
                read_error: Some(ReadError { below: None, above: Some(0.0), proven: true, how: "its rays run through real walls, so a wall they miss only reads thicker" }),
                ..measured(Ok(t), Some(at), Some(h))
            },
        },
        _ => match v.voids(f, res, accuracy)? {
            None => Measured { proven: true, held: Some(Held::Outright), ..measured(Ok(0.0), None, None) },
            Some(c) => {
                let held = c.count().map(|_| if c.thin == 0 { Held::Outright } else { Held::SetAside { cells: c.thin, at: c.thin_at } });
                let mut m = Measured { proven: true, held, ..measured(Ok(c.at_least as f64), c.count().is_none().then_some(c.at).flatten(), None) };
                if let Some((what, decide)) = c.why_undecided() {
                    m.excluded.push(Excluded {
                        possible: false,
                        key: "undecided-cells",
                        value: c.at_most.map_or(f64::INFINITY, |n| n as f64),
                        at: c.at,
                        what,
                        decide,
                        singular: false,
                    });
                }
                m
            }
        },
    })
}
