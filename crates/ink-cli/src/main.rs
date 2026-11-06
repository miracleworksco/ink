//! The kernel's `ink`: the command line (this crate's library) with the
//! kernel's packages.
fn main() {
    ink_cli::main(ink::registry(), ink::SPEC)
}
