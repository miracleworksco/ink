//! `ink` with ext's packages (sim, dynamics, topology): the kernel's command line
//! (ink-cli) given this distribution's registry and language reference.
fn main() {
    ink_cli::main(ink_ext::registry(), ink_ext::SPEC, concat!("ink-ext ", env!("CARGO_PKG_VERSION")), Some(ink_ext::solvers))
}
