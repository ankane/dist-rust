use core::f64::consts::{E, PI, SQRT_2};
use crate::math::{erf, fabs, log, pow, sqrt};

pub struct Normal;

impl Normal {
    /// Returns the probability density function (PDF) of the normal distribution.
    pub fn pdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        if std_dev <= 0.0 {
            return f64::NAN;
        }

        let n = (x - mean) / std_dev;
        (1.0 / (std_dev * sqrt(2.0 * PI))) * pow(E, -0.5 * n * n)
    }

    /// Returns the cumulative distribution function (CDF) of the normal distribution.
    pub fn cdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        if std_dev <= 0.0 {
            return f64::NAN;
        }

        0.5 * (1.0 + erf((x - mean) / (std_dev * SQRT_2)))
    }

    /// Returns the percent-point/quantile function (PPF) of the normal distribution.
    // Wichura, M. J. (1988).
    // Algorithm AS 241: The Percentage Points of the Normal Distribution.
    // Journal of the Royal Statistical Society. Series C (Applied Statistics), 37(3), 477-484.
    #[allow(clippy::excessive_precision)]
    pub fn ppf(p: f64, mean: f64, std_dev: f64) -> f64 {
        if !(0.0..=1.0).contains(&p) || std_dev <= 0.0 || mean.is_nan() || std_dev.is_nan() {
            return f64::NAN;
        }

        if p == 0.0 {
            return f64::NEG_INFINITY;
        }

        if p == 1.0 {
            return f64::INFINITY;
        }

        let q = p - 0.5;
        if fabs(q) < 0.425 {
            let r = 0.180625 - q * q;
            mean + std_dev * q *
                (((((((2.5090809287301226727e3 * r + 3.3430575583588128105e4) * r + 6.7265770927008700853e4) * r + 4.5921953931549871457e4) * r + 1.3731693765509461125e4) * r + 1.9715909503065514427e3) * r + 1.3314166789178437745e2) * r + 3.3871328727963666080e0) /
                (((((((5.2264952788528545610e3 * r + 2.8729085735721942674e4) * r + 3.9307895800092710610e4) * r + 2.1213794301586595867e4) * r + 5.3941960214247511077e3) * r + 6.8718700749205790830e2) * r + 4.2313330701600911252e1) * r + 1.0)
        } else {
            let mut r = if q < 0.0 { p } else { 1.0 - p };
            r = sqrt(-log(r));
            let sign = if q < 0.0 { -1.0 } else { 1.0 };
            if r < 5.0 {
                r -= 1.6;
                mean + std_dev * sign *
                    (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r + 2.41780725177450611770e-1) * r + 1.27045825245236838258e0) * r + 3.64784832476320460504e0) * r + 5.76949722146069140550e0) * r + 4.63033784615654529590e0) * r + 1.42343711074968357734e0) /
                    (((((((1.05075007164441684324e-9 * r + 5.47593808499534494600e-4) * r + 1.51986665636164571966e-2) * r + 1.48103976427480074590e-1) * r + 6.89767334985100004550e-1) * r + 1.67638483018380384940e0) * r + 2.05319162663775882187e0) * r + 1.0)
            } else {
                r -= 5.0;
                mean + std_dev * sign *
                    (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 1.24266094738807843860e-3) * r + 2.65321895265761230930e-2) * r + 2.96560571828504891230e-1) * r + 1.78482653991729133580e0) * r + 5.46378491116411436990e0) * r + 6.65790464350110377720e0) /
                    (((((((2.04426310338993978564e-15 * r + 1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r + 7.86869131145613259100e-4) * r + 1.48753612908506148525e-2) * r + 1.36929880922735805310e-1) * r + 5.99832206555887937690e-1) * r + 1.0)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Normal;
    use core::f64::{INFINITY, NEG_INFINITY};

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
    fn test_pdf_infinite_mean() {
        assert_in_delta(Normal::pdf(0.0, NEG_INFINITY, 1.0), 0.0, 0.00001);
        assert_in_delta(Normal::pdf(0.0, INFINITY, 1.0), 0.0, 0.00001);
    }

    #[test]
    fn test_pdf_infinite_std_dev() {
        assert_in_delta(Normal::pdf(0.0, 0.0, INFINITY), 0.0, 0.00001);
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
        assert!(Normal::pdf(0.0, 0.0, -1.0).is_nan());
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
    fn test_cdf_infinite_mean() {
        assert_in_delta(Normal::cdf(1.0, NEG_INFINITY, 1.0), 1.0, 0.00001);
        assert_in_delta(Normal::cdf(1.0, INFINITY, 1.0), 0.0, 0.00001);
    }

    #[test]
    fn test_cdf_infinite_std_dev() {
        assert_in_delta(Normal::cdf(1.0, 0.0, INFINITY), 0.5, 0.00001);
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
        assert!(Normal::cdf(0.0, 0.0, -1.0).is_nan());
    }

    #[test]
    fn test_ppf() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.28155, -0.84162, -0.5244, -0.25335, 0.0, 0.25335, 0.5244, 0.84162, 1.28155, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::ppf(*input, 0.0, 1.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_ppf_test_data() {
        // test data from paper
        assert_in_delta(Normal::ppf(0.25, 0.0, 1.0), -0.6744897501960817, 0.0000000000000001);
        assert_in_delta(Normal::ppf(0.001, 0.0, 1.0), -3.090232306167814, 0.000000000000001);
        assert_in_delta(Normal::ppf(1e-20, 0.0, 1.0), -9.262340089798408, 0.000000000000004);
    }

    #[test]
    fn test_ppf_mean_std_dev() {
        let inputs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let expected = [NEG_INFINITY, -1.5631, -0.68324, -0.0488, 0.49331, 1.0, 1.50669, 2.0488, 2.68324, 3.5631, INFINITY];
        for (input, exp) in inputs.iter().zip(expected) {
            assert_in_delta(Normal::ppf(*input, 1.0, 2.0), exp, 0.00001);
        }
    }

    #[test]
    fn test_ppf_nan() {
        assert!(Normal::ppf(f64::NAN, 0.0, 1.0).is_nan());
        assert!(Normal::ppf(0.0, f64::NAN, 1.0).is_nan());
        assert!(Normal::ppf(0.0, 0.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_ppf_negative_p() {
        assert!(Normal::ppf(-1.0, 0.0, 1.0).is_nan());
    }

    #[test]
    fn test_ppf_zero_std_dev() {
        assert!(Normal::ppf(0.5, 0.0, 0.0).is_nan());
    }

    #[test]
    fn test_ppf_negative_std_dev() {
        assert!(Normal::ppf(0.5, 0.0, -1.0).is_nan());
    }
}
