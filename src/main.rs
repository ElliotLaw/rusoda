// use crate::array_::DArray;

mod arrops;
mod blas;
mod common;
mod ddot_test;
mod lsoda;
mod lsoda_test;

// mod linalg_func;
// use std::ptr::NonNull;

fn main() {
    lsoda_test::oral_1cpt_test()
    // ddot_test::ddot_test()
}
