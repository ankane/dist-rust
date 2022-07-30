use crate::erf::inverse_erf;
use crate::math::erf;
use std::f64::consts::{E, PI, SQRT_2};

pub struct Normal;

impl Normal {
    pub fn pdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        // TODO uncomment in 0.2.0
        // if std_dev <= 0.0 {
        //     return f64::NAN;
        // }

        let n = (x - mean) / std_dev;
        (1.0 / (std_dev * (2.0 * PI).sqrt())) * E.powf(-0.5 * n * n)
    }

    pub fn cdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        // TODO uncomment in 0.2.0
        // if std_dev <= 0.0 {
        //     return f64::NAN;
        // }

        0.5 * (1.0 + erf((x - mean) / (std_dev * SQRT_2)))
    }

    pub fn ppf(p: f64, mean: f64, std_dev: f64) -> f64 {
        assert!(p >= 0.0 && p <= 1.0);
        // TODO uncomment in 0.2.0
        // if p < 0.0 || p > 1.0 || std_dev <= 0.0 {
        //     return f64::NAN;
        // }

        mean + std_dev * SQRT_2 * inverse_erf(2.0 * p - 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::Normal;
    use std::f64::{INFINITY, NEG_INFINITY};

    fn assert_in_delta(act: f64, exp: f64, delta: f64) {
        if exp.is_finite() {
            assert!((exp - act).abs() < delta, "{} != {}", act, exp);
        } else {
            assert_eq!(act, exp);
        }
    }

    #[test]
    fn test_pdf() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.00443, 0.05399, 0.24197, 0.39894, 0.24197, 0.05399, 0.00443, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::pdf(*input, 0.0, 1.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_mean_std_dev() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.027, 0.06476, 0.12099, 0.17603, 0.19947, 0.17603, 0.12099, 0.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::pdf(*input, 1.0, 2.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_nan() {
        assert!(Normal::pdf(f64::NAN, 0.0, 1.0).is_nan());
        assert!(Normal::pdf(0.0, f64::NAN, 1.0).is_nan());
        assert!(Normal::pdf(0.0, 0.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_pdf_zero_std_dev() {
        assert!(Normal::pdf(0.0, 0.0, 0.0).is_nan());
    }

    #[test]
    fn test_pdf_negative_std_dev() {
        // TODO return NAN in 0.2.0
        assert!(Normal::pdf(0.0, 0.0, -1.0).is_sign_negative());
    }

    #[test]
    fn test_cdf() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.00135, 0.02275, 0.15866, 0.5, 0.84134, 0.97725, 0.99865, 1.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::cdf(*input, 0.0, 1.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf_mean_std_dev() {
        let inputs = [NEG_INFINITY, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, INFINITY];
        let expected = [0.0, 0.02275, 0.06681, 0.15866, 0.30854, 0.5, 0.69146, 0.84134, 1.0];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::cdf(*input, 1.0, 2.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf_nan() {
        assert!(Normal::cdf(f64::NAN, 0.0, 1.0).is_nan());
        assert!(Normal::cdf(0.0, f64::NAN, 1.0).is_nan());
        assert!(Normal::cdf(0.0, 0.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_cdf_zero_std_dev() {
        assert!(Normal::cdf(0.0, 0.0, 0.0).is_nan());
    }

    #[test]
    fn test_cdf_negative_std_dev() {
        // TODO return NAN in 0.2.0
        assert_in_delta(Normal::cdf(0.0, 0.0, -1.0), 0.5, 0.00001);
    }

    #[test]
    fn test_ppf() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.28155, -0.84162, -0.5244, -0.25335, 0.0, 0.25335, 0.5244, 0.84162, 1.28155, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::ppf(*input, 0.0, 1.0), exp, 0.0002);
        }
    }

    #[test]
    fn test_ppf_mean_std_dev() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.5631, -0.68324, -0.0488, 0.49331, 1.0, 1.50669, 2.0488, 2.68324, 3.5631, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::ppf(*input, 1.0, 2.0), exp, 0.0004);
        }
    }

    #[test]
    fn test_ppf_nan() {
        // TODO uncomment in 0.2.0
        // assert!(Normal::ppf(f64::NAN, 0.0, 1.0).is_nan());
        assert!(Normal::ppf(0.0, f64::NAN, 1.0).is_nan());
        assert!(Normal::ppf(0.0, 0.0, f64::NAN).is_nan());
    }

    #[test]
    #[should_panic(expected = "assertion failed: p >= 0.0 && p <= 1.0")]
    fn test_ppf_negative_p() {
        // TODO return NAN in 0.2.0
        Normal::ppf(-1.0, 0.0, 1.0);
    }

    #[test]
    fn test_ppf_zero_std_dev() {
        // TODO return NAN in 0.2.0
        assert_in_delta(Normal::ppf(0.5, 0.0, 0.0), 0.0, 0.00001);
    }

    #[test]
    fn test_ppf_negative_std_dev() {
        // TODO return NAN in 0.2.0
        assert_in_delta(Normal::ppf(0.5, 0.0, -1.0), 0.0, 0.00001);
    }
}
