// Winitzki, S. (2008).
// A handy approximation for the error function and its inverse.
// https://drive.google.com/file/d/0B2Mt7luZYBrwZlctV3A3eF82VGM/view?resourcekey=0-UQpPhwZgzP0sF4LHBDlLtg
// from https://sites.google.com/site/winitzki

use core::f64::consts::PI;
use libm::{log, sqrt};

pub fn inverse_erf(x: f64) -> f64 {
    let (sign, x) = if x < 0.0 {
        (-1.0, -x)
    } else {
        (1.0, x)
    };

    let a = 0.147;
    let ln = log(1.0 - x * x);
    let f1 = 2.0 / (PI * a);
    let f2 = ln / 2.0;
    let f3 = f1 + f2;
    let f4 = 1.0 / a * ln;
    sign * sqrt(-f1 - f2 + sqrt(f3 * f3 - f4))
}
