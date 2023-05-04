use extendr_api::prelude::*;
use std::collections::HashMap;

fn initialize_hashmap_vector(n: usize) -> Vec<HashMap<String, i32>> {
    let mut vector = Vec::with_capacity(n);
    for _ in 1..n {
        vector.push(HashMap::new());
    }
    vector
}

#[extendr]
fn rust_get_clone_sizes(clusters: Integers, clonotype_ids: Strings) {
    let mut arr = initialize_hashmap_vector(clusters.len());
    rprintln!("unfinished");
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod APackOfTheClones;
    //fn rust_get_clone_sizes;
}
