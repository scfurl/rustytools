use extendr_api::prelude::*;

mod getconsensus;


/// @export
#[extendr]
fn getconsensus(rstring: Robj, index_add: Robj)  {
    // let rstring = "NNNNNNNNATGCGTANNNNNNNTTGACNNNNNNNN".to_string();
    let rstring = rstring.as_str().unwrap().to_string();
    let index_add = *index_add.as_real_vector().unwrap().first().unwrap() as usize;
    getconsensus::getconsensus(rstring, index_add);
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rustytools;
    fn getconsensus;
}
