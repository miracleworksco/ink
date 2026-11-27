use super::point::GRAD_INSTRUCTIONS;
use super::{Ctx, L, eval_grads, eval_grads_at};
use crate::graph::{Graph, Var};
use crate::tape::Tape;

#[test]
fn intermediate_gradients_share_one_pass_per_batch() {
    let mut graph = Graph::new();
    let x = graph.var(Var::X);
    let mut f = x;
    for i in 1..16 {
        let k = graph.constant(i as f64);
        let shifted = graph.add(x, k);
        f = graph.mul(f, shifted);
    }
    let tape = Tape::new(&graph, f);
    let instructions: Vec<_> = (0..tape.len()).step_by(3).collect();
    let pts = vec![[0.5; 3]; L + 1];
    let ctx = Ctx::new(0.0, vec![]);
    GRAD_INSTRUCTIONS.set(0);
    let rows = eval_grads_at(&tape, &ctx, &pts, &instructions);
    assert_eq!(GRAD_INSTRUCTIONS.get(), 2 * (instructions.last().unwrap() + 1));
    for (&i, row) in instructions.iter().zip(rows) {
        assert_eq!(row, eval_grads(&tape.cone(i), &ctx, &pts));
    }
}
