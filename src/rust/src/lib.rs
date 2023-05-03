use extendr_api::prelude::*;

#[extendr]
fn placeholder() {
    rprintln!("placeholder");
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod APackOfTheClones;
    fn placeholder;
}