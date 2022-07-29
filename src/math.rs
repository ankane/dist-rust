extern "C" {
    fn erf(x: f64) -> f64;
    fn tgamma(x: f64) -> f64;
}

#[inline]
pub fn erf2(x: f64) -> f64 {
    unsafe { erf(x) }
}

#[inline]
pub fn gamma(x: f64) -> f64 {
    unsafe { tgamma(x) }
}
