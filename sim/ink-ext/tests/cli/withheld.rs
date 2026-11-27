//! What this release holds back (src/withheld.rs) is refused by name: a
//! model using it, or its command, reads that it is not in this release,
//! never that it is unknown.

use crate::common::*;

#[test]
fn flow_is_refused_as_not_in_this_release() {
    let s = Scratch::new("withheld-flow");
    let f = s.model(
        "pipe.ink",
        "(import \"std\")\n(def bore (max (- (len2 y z) 2) (- (abs x) 10)))\n\
         (flow feed (shape bore) (fluid oil) (inlet in (region (+ x 9)) (pressure 20Pa)) (outlet out (region (- 9 x))))\n",
    );
    for args in [vec!["check", f.as_str()], vec!["flow", f.as_str()]] {
        let e = s.ink(&args).err();
        let says = e["error"]["message"].as_str().unwrap_or_default();
        assert!(says.contains("is not in this release of ink") && says.contains("flow"), "ink {args:?}: {}", pretty(&e));
    }
}

#[test]
fn least_time_is_refused_as_not_in_this_release() {
    let s = Scratch::new("withheld-least-time");
    let f = s.model(
        "cart.ink",
        "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(body cart (shape (box 100 100 100)) (attr density 1g/cm3))\n\
         (assembly (axis (x 0)) (instance c cart (at) (slide x x) (effort-limits -10N 10N)))\n\
         (assembly root (instance a (axis (* 300 (/ t 1s)))))\n(require quick (value (least-time) < 0.7s))\n",
    );
    // Its option of `ink effort` too, on a file that loads.
    let plain = s.model("plain.ink", &std::fs::read_to_string(&f).unwrap().replace("(require quick (value (least-time) < 0.7s))", ""));
    for args in [vec!["check", f.as_str()], vec!["effort", plain.as_str(), "--least-time"]] {
        let e = s.ink(&args).err();
        let says = e["error"]["message"].as_str().unwrap_or_default();
        assert!(says.contains("is not in this release of ink") && says.contains("least time"), "ink {args:?}: {}", pretty(&e));
    }
}
