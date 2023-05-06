use extendr_api::prelude::*;

// just a test
#[extendr]
fn rust_sum(component_vector: Integers) -> Rint {
    component_vector.elt(0) + component_vector.elt(1)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod APackOfTheClones;
    fn rust_sum;
}
