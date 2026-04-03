#![cfg_attr(not(feature = "nightly"), allow(unsafe_code))]

#[cfg(not(feature = "nightly"))]
mod c {
    extern "C" {
        pub fn erf(x: f64) -> f64;
        pub fn tgamma(x: f64) -> f64;
    }
}

#[inline]
pub fn atan(x: f64) -> f64 {
    x.atan()
}

#[inline]
pub fn cos(x: f64) -> f64 {
    x.cos()
}

#[cfg(feature = "nightly")]
#[inline]
pub fn erf(x: f64) -> f64 {
    x.erf()
}

#[cfg(not(feature = "nightly"))]
#[inline]
pub fn erf(x: f64) -> f64 {
    unsafe { c::erf(x) }
}

#[inline]
pub fn exp(x: f64) -> f64 {
    x.exp()
}

#[inline]
pub fn fabs(x: f64) -> f64 {
    x.abs()
}

#[inline]
pub fn floor(x: f64) -> f64 {
    x.floor()
}

#[inline]
pub fn log(x: f64) -> f64 {
    x.ln()
}

#[inline]
pub fn pow(x: f64, y: f64) -> f64 {
    x.powf(y)
}

#[inline]
pub fn sin(x: f64) -> f64 {
    x.sin()
}

#[inline]
pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

#[cfg(feature = "nightly")]
#[inline]
pub fn tgamma(x: f64) -> f64 {
    x.gamma()
}

#[cfg(not(feature = "nightly"))]
#[inline]
pub fn tgamma(x: f64) -> f64 {
    unsafe { c::tgamma(x) }
}
