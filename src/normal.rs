use crate::erf::inverse_erf;
use crate::math::erf;
use std::f64::consts::{E, PI};

pub struct Normal;

impl Normal {
    pub fn pdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        (1.0 / (std_dev * (2.0 * PI).sqrt())) * E.powf(-0.5 * ((x - mean) / std_dev).powf(2.0))
    }

    pub fn cdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        0.5 * (1.0 + erf((x - mean) / (std_dev * 2.0_f64.sqrt())))
    }

    pub fn ppf(p: f64, mean: f64, std_dev: f64) -> f64 {
        assert!(p >= 0.0 && p <= 1.0);

        mean + std_dev * 2.0_f64.sqrt() * inverse_erf(2.0 * p - 1.0)
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
        let inputs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let expected = [0.00443, 0.05399, 0.24197, 0.39894, 0.24197, 0.05399, 0.00443];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::pdf(*input, 0.0, 1.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_pdf_mean_std_dev() {
        let inputs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let expected = [0.027, 0.06476, 0.12099, 0.17603, 0.19947, 0.17603, 0.12099];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::pdf(*input, 1.0, 2.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf() {
        let inputs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let expected = [0.00135, 0.02275, 0.15866, 0.5, 0.84134, 0.97725, 0.99865];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::cdf(*input, 0.0, 1.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_cdf_mean_std_dev() {
        let inputs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let expected = [0.02275, 0.06681, 0.15866, 0.30854, 0.5, 0.69146, 0.84134];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::cdf(*input, 1.0, 2.0), exp, 0.00001);
        }
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
    #[should_panic(expected = "assertion failed: p >= 0.0 && p <= 1.0")]
    fn test_ppf_negative_p() {
        Normal::ppf(-1.0, 0.0, 1.0);
    }
}
