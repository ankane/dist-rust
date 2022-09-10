#![allow(unsafe_code)]

mod c {
    extern "C" {
        pub fn erf(x: f64) -> f64;
        pub fn tgamma(x: f64) -> f64;
    }
}

#[inline]
pub fn erf(x: f64) -> f64 {
    unsafe { c::erf(x) }
}

#[inline]
pub fn tgamma(x: f64) -> f64 {
    unsafe { c::tgamma(x) }
}
