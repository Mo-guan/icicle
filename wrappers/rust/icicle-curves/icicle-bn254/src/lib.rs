#[macro_use(impl_curve)]
extern crate icicle_core;
#[macro_use]
#[cfg(test)]
extern crate lazy_static;

pub mod curve;
pub mod msm;
pub mod ntt;
