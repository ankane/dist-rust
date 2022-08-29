use crate::math::tgamma;
use crate::Normal;
use std::f64::consts::PI;

pub struct StudentsT;

impl StudentsT {
    pub fn pdf<T: Into<f64>>(x: f64, n: T) -> f64 {
        let n = n.into();

        if n.is_nan() || n <= 0.0 {
            return f64::NAN;
        }

        if n == f64::INFINITY {
            return Normal::pdf(x, 0.0, 1.0);
        }

        tgamma((n + 1.0) / 2.0) / ((n * PI).sqrt() * tgamma(n / 2.0)) * (1.0 + x * x / n).powf(-(n + 1.0) / 2.0)
    }

    // Hill, G. W. (1970).
    // Algorithm 395: Student's t-distribution.
    // Communications of the ACM, 13(10), 617-619.
    pub fn cdf<T: Into<f64>>(x: f64, n: T) -> f64 {
        let n = n.into();

        // TODO support n > 0.0
        if x.is_nan() || n.is_nan() || n < 1.0 {
            return f64::NAN;
        }

        if x == f64::NEG_INFINITY {
            return 0.0;
        }

        if x == f64::INFINITY {
            return 1.0;
        }

        if n == f64::INFINITY {
            return Normal::cdf(x, 0.0, 1.0);
        }

        let (start, sign) = if x < 0.0 {
            (0.0, 1.0)
        } else {
            (1.0, -1.0)
        };

        let mut z = 1.0;
        let t = x * x;
        let mut y = t / n;
        let mut b = 1.0 + y;

        if n > n.floor() || (n >= 20.0 && t < n) || n > 200.0 {
            // asymptotic series for large or noninteger n
            if y > 10e-6 {
                y = b.ln();
            }
            let a = n - 0.5;
            b = 48.0 * a * a;
            y *= a;
            y = (((((-0.4 * y - 3.3) * y - 24.0) * y - 85.5) / (0.8 * y * y + 100.0 + b) + y + 3.0) / b + 1.0) * y.sqrt();
            return start + sign * Normal::cdf(-y, 0.0, 1.0);
        }

        // make n mutable and int
        // n is int between 1 and 200 if made it here
        let mut n = n as u8;

        if n < 20 && t < 4.0 {
            // nested summation of cosine series
            y = y.sqrt();
            let mut a = y;
            if n == 1 {
                a = 0.0;
            }

            // loop
            if n > 1 {
                n -= 2;
                while n > 1 {
                    a = (n - 1) as f64 / (b * n as f64) * a + y;
                    n -= 2;
                }
            }
            a = if n == 0 { a / b.sqrt() } else { (y.atan() + a / b) * (2.0 / PI) };
            return start + sign * (z - a) / 2.0;
        }

        // tail series expanation for large t-values
        let mut a = b.sqrt();
        y = a * n as f64;
        let mut j = 0;
        while a != z {
            j += 2;
            z = a;
            y = y * (j - 1) as f64 / (b * j as f64);
            a += y / (n + j) as f64;
        }
        z = 0.0;
        y = 0.0;
        a = -a;

        // loop (without n + 2 and n - 2)
        while n > 1 {
            a = (n - 1) as f64 / (b * n as f64) * a + y;
            n -= 2;
        }
        a = if n == 0 { a / b.sqrt() } else { (y.atan() + a / b) * (2.0 / PI) };
        start + sign * (z - a) / 2.0
    }

    // Hill, G. W. (1970).
    // Algorithm 396: Student's t-quantiles.
    // Communications of the ACM, 13(10), 619-620.
    pub fn ppf<T: Into<f64>>(p: f64, n: T) -> f64 {
        let n = n.into();

        // TODO support n > 0.0
        if !(0.0..=1.0).contains(&p) || n < 1.0 {
            return f64::NAN;
        }

        if n == f64::INFINITY {
            return Normal::ppf(p, 0.0, 1.0);
        }

        // distribution is symmetric
        let (sign, p) = if p < 0.5 {
            (-1.0, 1.0 - p)
        } else {
            (1.0, p)
        };

        // two-tail to one-tail
        let p = 2.0 * (1.0 - p);

        if n == 2.0 {
            return sign * (2.0 / (p * (2.0 - p)) - 2.0).sqrt();
        }

        let half_pi = PI / 2.0;

        if n == 1.0 {
            let p = p * half_pi;
            return sign * p.cos() / p.sin();
        }

        let a = 1.0 / (n - 0.5);
        let b = 48.0 / (a * a);
        let mut c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
        let d = ((94.5 / (b + c) - 3.0) / b + 1.0) * (a * half_pi).sqrt() * n;
        let mut x = d * p;
        let mut y = x.powf(2.0 / n);
        if y > 0.05 + a {
            // asymptotic inverse expansion about normal
            x = Normal::ppf(p * 0.5, 0.0, 1.0);
            y = x * x;
            if n < 5.0 {
                c += 0.3 * (n - 4.5) * (x + 0.6);
            }
            c += (((0.05 * d * x - 5.0) * x - 7.0) * x - 2.0) * x + b;
            y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
            y = a * y * y;
            y = if y > 0.002 { y.exp() - 1.0 } else { 0.5 * y * y + y };
        } else {
            y = ((1.0 / (((n + 6.0) / (n * y) - 0.089 * d - 0.822) * (n + 2.0) * 3.0) + 0.5 / (n + 4.0)) * y - 1.0) * (n + 1.0) / (n + 2.0) + 1.0 / y;
        }
        sign * (n * y).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::StudentsT;
    use std::f64::{INFINITY, NEG_INFINITY};

    fn assert_in_delta(act: f64, exp: f64, delta: f64) {
        if exp.is_finite() {
            assert!((exp - act).abs() < delta, "{} != {}", act, exp);
        } else {
            assert_eq!(act, exp);
        }
    }

    #[test]
    fn test_pdf_one() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.03183, 0.06366, 0.15915, 0.31831, 0.15915, 0.06366, 0.03183, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::pdf(*input, 1), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_two() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.02741, 0.06804, 0.19245, 0.35355, 0.19245, 0.06804, 0.02741, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::pdf(*input, 2), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_thirty() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.00678, 0.05685, 0.23799, 0.39563, 0.23799, 0.05685, 0.00678, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::pdf(*input, 30), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_non_integer() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.02504, 0.06796, 0.2008, 0.36181, 0.2008, 0.06796, 0.02504, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::pdf(*input, 2.5), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_infinity() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.00443, 0.05399, 0.24197, 0.39894, 0.24197, 0.05399, 0.00443, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::pdf(*input, INFINITY), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_less_than_one() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.02963, 0.0519, 0.1183, 0.26968, 0.1183, 0.0519, 0.02963, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::pdf(*input, 0.5), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_nan() {
        assert!(StudentsT::pdf(f64::NAN, 1).is_nan());
        assert!(StudentsT::pdf(0.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_pdf_zero_n() {
        assert!(StudentsT::pdf(0.5, 0).is_nan());
    }

    #[test]
    fn test_cdf_one() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.10242, 0.14758, 0.25, 0.5, 0.75, 0.85242, 0.89758, 1.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::cdf(*input, 1), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf_two() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.04773, 0.09175, 0.21132, 0.5, 0.78868, 0.90825, 0.95227, 1.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::cdf(*input, 2), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf_thirty() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.00269, 0.02731, 0.16265, 0.5, 0.83735, 0.97269, 0.99731, 1.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::cdf(*input, 30), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf_non_integer() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.03629, 0.0787, 0.20203, 0.5, 0.79797, 0.9213, 0.96371, 1.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::cdf(*input, 2.5), exp, 0.00005);
        }
    }

    #[test]
    fn test_cdf_infinity() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.00135, 0.02275, 0.15866, 0.5, 0.84134, 0.97725, 0.99865, 1.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::cdf(*input, INFINITY), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf_nan() {
        assert!(StudentsT::cdf(f64::NAN, 1.0).is_nan());
        assert!(StudentsT::cdf(0.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_cdf_zero_n() {
        assert!(StudentsT::cdf(0.5, 0).is_nan());
    }

    #[test]
    fn test_ppf_one() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -3.07768, -1.37638, -0.72654, -0.32492, 0.0, 0.32492, 0.72654, 1.37638, 3.07768, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::ppf(*input, 1), exp, 0.00001);
        }
    }

    #[test]
    fn test_ppf_two() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.88562, -1.06066, -0.61721, -0.28868, 0.0, 0.28868, 0.61721, 1.06066, 1.88562, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::ppf(*input, 2), exp, 0.00001);
        }
    }

    #[test]
    fn test_ppf_thirty() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.31042, -0.85377, -0.53002, -0.25561, 0.0, 0.25561, 0.53002, 0.85377, 1.31042, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::ppf(*input, 30), exp, 0.00001);
        }
    }

    #[test]
    fn test_ppf_non_integer() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.73025, -1.01016, -0.59731, -0.28146, 0.0, 0.28146, 0.59731, 1.01016, 1.73025, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::ppf(*input, 2.5), exp, 0.0002);
        }
    }

    #[test]
    fn test_ppf_infinity() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.28155, -0.84162, -0.5244, -0.25335, 0.0, 0.25335, 0.5244, 0.84162, 1.28155, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(StudentsT::ppf(*input, INFINITY), exp, 0.00001);
        }
    }

    #[test]
    fn test_ppf_nan() {
        assert!(StudentsT::ppf(f64::NAN, 1.0).is_nan());
        assert!(StudentsT::ppf(0.5, f64::NAN).is_nan());
    }

    #[test]
    fn test_ppf_negative_p() {
        assert!(StudentsT::ppf(-1.0, 1).is_nan());
    }

    #[test]
    fn test_ppf_zero_n() {
        assert!(StudentsT::ppf(0.5, 0).is_nan());
    }
}
