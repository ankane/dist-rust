// Winitzki, S. (2008).
// A handy approximation for the error function and its inverse.
// https://drive.google.com/file/d/0B2Mt7luZYBrwZlctV3A3eF82VGM/view?resourcekey=0-UQpPhwZgzP0sF4LHBDlLtg
// from https://sites.google.com/site/winitzki

use std::f64::consts::PI;

pub fn erf(x: f64) -> f64 {
    let (sign, x) = if x < 0.0 {
        (-1.0, -x)
    } else {
        (1.0, x)
    };

    let a = 0.14;
    let x2 = x * x;
    sign * (1.0 - (-x2 * (4.0 / PI + a * x2) / (1.0 + a * x2)).exp()).sqrt()
}

pub fn inverse_erf(x: f64) -> f64 {
    let (sign, x) = if x < 0.0 {
        (-1.0, -x)
    } else {
        (1.0, x)
    };

    let a = 0.147;
    let ln = (1.0 - x * x).ln();
    let f1 = 2.0 / (PI * a);
    let f2 = ln / 2.0;
    let f3 = f1 + f2;
    let f4 = 1.0 / a * ln;
    sign * (-f1 - f2 + (f3 * f3 - f4).sqrt()).sqrt()
}
