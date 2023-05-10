use extendr_api::prelude::*;
use extendr_api::wrapper::List;

#[extendr]
fn distV(c1: List, c2: List) -> Doubles {
    // extract centroid vectors (4th element)
    let c1_centroid = c1.elt(3).unwrap();
    let c2_centroid = c2.elt(3).unwrap();

    // every element of a List is an Robj which you have
    // to convert into the appropriate type
    let c1_centroid_vec = Doubles::try_from(c1_centroid).unwrap();
    let c2_centroid_vec = Doubles::try_from(c2_centroid).unwrap();

    // there is no vector operations in rust; work in "iterators"
    c1_centroid_vec
        .iter()
        // combine x and y into one iterator to iterate over it together
        .zip(c2_centroid_vec.iter())
        // take only the first 2 iterator items
        .take(2)
        // apply some operation to each
        .map(|(xi, yi)| {
            xi - yi
        })
        // collect into a Double vec
        .collect::<Doubles>()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod APackOfTheClones;
    fn distV;
}
