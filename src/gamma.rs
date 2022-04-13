//! Provides an approximation of the gamma function.
//!
//! Based on the Python implementation shown at <https://en.wikipedia.org/wiki/Lanczos_approximation>
//! which is based on the Lanczos approximation of the gamma function.

use std::f64::consts::PI;

static COEFFICIENTS: &'static [f64] = &[
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
];

/// Returns the approximate value of the [gamma function](https://en.wikipedia.org/wiki/Gamma_function)
/// for the given argument, or `None` if the gamma function has no defined value for the given argument.
///
/// # Examples
///
/// This function will return `Some(y)` for arguments which can be mapped to a value by the
/// gamma function, otherwise it will return `None`.
/// ```
/// use distrs::gamma;
/// # let x = 1.0;
/// match gamma::calculate(x) {
///     Some(y) => println!("gamma({}) gives {}", x, y),
///     None => println!("gamma({}) is not defined", x)
/// };
/// ```
/// The gamma function does not have defined values for zero nor for negative integers, so for these
/// arguments `None` will be returned:
/// ```
/// # use distrs::gamma;
/// assert_eq!(None, gamma::calculate(0.0));
/// assert_eq!(None, gamma::calculate(-1.0));
/// assert_eq!(None, gamma::calculate(-2.0));
/// assert_eq!(None, gamma::calculate(-3.0));
/// ```
///
/// # Warning
///
/// This function is based on an approximation (specifically the
/// [Lanczos approximation](https://en.wikipedia.org/wiki/Lanczos_approximation)) and does not
/// calculate exact (perfectly accurate) values. Take a look at the unit tests built into the
/// gamma.rs source code file and you'll notice that the the absolute accuracy (against precise
/// values reported in
/// [this Wikipedia page](https://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function))
/// is good for arguments between zero and ten, but begins to decrease after that. The proportional
/// error (discrepancy over expected value) never seems to become high, however. But it's possible
/// that accuracy may not be good for arguments in-between those tested. So approach with caution,
/// and contribute your own unit tests if you need to investigate accuracy in a particular area of
/// the argument domain.
///
/// Also, note that a trick has been used to check for zero and negative integer values (within the
/// given f64 argument) and it's possible that this will fail to detect negative integers in some
/// cases. So apply your own pre-checks if you need to catch cases where negative integer arguments
/// must not be passed to this function.
pub fn calculate<T: Into<f64>>(x: T) -> Option<f64> {
    let x = x.into();
    return if gamma_undefined_for(x) {
        Option::None
    } else if x < 0.5 {
        Some(reflection_formula(x))
    } else {
        Some(lanczos_gamma(x))
    };
}

fn gamma_undefined_for(x: f64) -> bool {
    // The gamma function is not defined for zero nor for negative integer arguments.
    // Use `x.trunc() == x` to check whether the f64 argument is a whole number. (Unit testing
    // suggests that this works, but it's possible that the floating point representation may
    // cause this trick to break for some values.)
    x < 0.5 && x.trunc() == x
}

// The reflection formula must be used when (the real number) x is less than 0.5.
fn reflection_formula(x: f64) -> f64 {
    return PI / ((PI * x).sin() * lanczos_gamma(1.0 - x));
}

// This function only works for arguments of 0.5 or greater.
fn lanczos_gamma(z: f64) -> f64 {
    let z = z - 1.0;
    let mut x = 0.99999999999980993;
    for i in 0..COEFFICIENTS.len() {
        let pval = COEFFICIENTS[i];
        let i = i as f64;
        x += pval / (z + i + 1.0);
    }
    let t = z - 0.5 + COEFFICIENTS.len() as f64;
    let y = (2.0 * PI).sqrt() * (t).powf(z + 0.5) * (-t).exp() * x;
    return y;
}

#[cfg(test)]
mod tests {
    use crate::gamma;
    use std::f64::consts::PI;

    // The greatest acceptable error proportion (between actual and expected values) acceptable by
    // the unit tests.
    const ACCEPTABLE_PROPORTIONAL_ERROR: f64 = 0.00000000000001;

    fn assert_within(expected: f64, actual: f64, delta: f64) {
        let diff = (actual - expected).abs();
        assert!(
            diff <= delta,
            "Absolute discrepancy too large: ({} - {}).abs() = {} which exceeds max delta of {}",
            expected,
            actual,
            diff,
            delta
        );
        let prop_error = (diff / expected).abs();
        assert!(
            prop_error < ACCEPTABLE_PROPORTIONAL_ERROR,
            "Proportional error ({}/{}).abs() = {} which exceeds {}",
            diff,
            expected,
            prop_error,
            ACCEPTABLE_PROPORTIONAL_ERROR
        );
    }

    // Gives an exact answer, but cannot take arguments greater than 20 before overflow.
    fn factorial(a: u64) -> u64 {
        if a == 0 {
            return 1;
        }
        if a <= 2 {
            return a;
        }
        let mut product: u64 = 1;
        for i in 2..=a {
            product *= i;
        }
        return product;
    }

    // Gives an approximate answer, but can handle much larger numbers.
    fn factorial_f64(a: u64) -> f64 {
        if a == 0 {
            return 1.0;
        }
        if a <= 2 {
            return a as f64;
        }
        let mut product: f64 = 1.0;
        for i in 2..=a {
            product *= i as f64;
        }
        return product;
    }

    // Expected values are based on https://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function

    #[test]
    fn positive_integers_tiny() {
        for i in 1..=5 {
            assert_within(
                factorial(i - 1) as f64,
                gamma::calculate(i as f64).unwrap(),
                0.00000000000001,
            );
        }
    }

    #[test]
    fn positive_integers_small() {
        for i in 6..=10 {
            assert_within(
                factorial(i - 1) as f64,
                gamma::calculate(i as f64).unwrap(),
                0.00000001,
            );
        }
    }

    #[test]
    fn positive_integers_low_teens() {
        for i in 11..=15 {
            assert_within(
                factorial(i - 1) as f64,
                gamma::calculate(i as f64).unwrap(),
                0.001,
            );
        }
    }

    #[test]
    fn positive_integers_high_teens() {
        assert_within(factorial(15) as f64, gamma::calculate(16.0).unwrap(), 0.01);
        assert_within(factorial(16) as f64, gamma::calculate(17.0).unwrap(), 0.1);
        assert_within(factorial(17) as f64, gamma::calculate(18.0).unwrap(), 0.5);
        assert_within(factorial(18) as f64, gamma::calculate(19.0).unwrap(), 10.0);
    }

    #[test]
    fn positive_integers_twenties() {
        assert_within(factorial(19) as f64, gamma::calculate(20.0).unwrap(), 224.0);
        assert_within(
            factorial(20) as f64,
            gamma::calculate(21.0).unwrap(),
            2048.0,
        );
        // Factorial(21) is too large a number to hold in a u64, so we need to use the factorial_f64 version.
        assert_within(
            factorial_f64(21) as f64,
            gamma::calculate(22.0).unwrap(),
            40960.0,
        );
        assert_within(
            factorial_f64(22) as f64,
            gamma::calculate(23.0).unwrap(),
            655360.0,
        );
        assert_within(
            factorial_f64(23) as f64,
            gamma::calculate(24.0).unwrap(),
            50331648.0,
        );
        assert_within(
            factorial_f64(24) as f64,
            gamma::calculate(25.0).unwrap(),
            402653184.0,
        );
        assert_within(
            factorial_f64(25) as f64,
            gamma::calculate(26.0).unwrap(),
            8589934592.0,
        );
        assert_within(
            factorial_f64(26) as f64,
            gamma::calculate(27.0).unwrap(),
            412316860416.0,
        );
        assert_within(
            factorial_f64(27) as f64,
            gamma::calculate(28.0).unwrap(),
            37383395344384.0,
        );
        assert_within(
            factorial_f64(28) as f64,
            gamma::calculate(29.0).unwrap(),
            // This error looks huge, but it's actually less than 1 in 2.9E14 of the expected result.
            1020346790576128.0,
        );
    }

    #[test]
    fn positive_half_integers_tiny() {
        for n in 0..=5 {
            let expected =
                PI.sqrt() * factorial(2 * n) as f64 / (4_u64.pow(n as u32) * factorial(n)) as f64;
            let actual = gamma::calculate(n as f64 + 0.5).unwrap();
            assert_within(expected, actual, 0.0000000000001);
        }
    }

    #[test]
    fn positive_half_integers_small() {
        for n in 6..=10 {
            let expected =
                PI.sqrt() * factorial(2 * n) as f64 / (4_u64.pow(n as u32) * factorial(n)) as f64;
            let actual = gamma::calculate(n as f64 + 0.5).unwrap();
            assert_within(expected, actual, 0.00000001);
        }
    }

    #[test]
    fn negative_half_integers_tiny() {
        assert_within(
            -3.5449077018110320546,
            gamma::calculate(-0.5).unwrap(),
            0.00000000000001,
        );
        assert_within(
            2.3632718012073547031,
            gamma::calculate(-1.5).unwrap(),
            0.00000000000001,
        );
        assert_within(
            -0.9453087204829418812,
            gamma::calculate(-2.5).unwrap(),
            0.000000000000001,
        );
    }

    #[test]
    fn positive_fractions() {
        assert_within(
            2.6789385347077476337,
            gamma::calculate(1.0 / 3.0).unwrap(),
            0.000000000000001,
        );
        assert_within(
            3.6256099082219083119,
            gamma::calculate(1.0 / 4.0).unwrap(),
            0.000000000000001,
        );
        assert_within(
            4.5908437119988030532,
            gamma::calculate(1.0 / 5.0).unwrap(),
            0.000000000000001,
        );
        assert_within(
            5.5663160017802352043,
            gamma::calculate(1.0 / 6.0).unwrap(),
            0.000000000000001,
        );
        assert_within(
            6.5480629402478244377,
            gamma::calculate(1.0 / 7.0).unwrap(),
            0.00000000000001,
        );
        assert_within(
            7.5339415987976119047,
            gamma::calculate(1.0 / 8.0).unwrap(),
            0.000000000000001,
        );
    }

    #[test]
    fn positive_real_axis_local_minimum() {
        assert_within(
            0.885603194410888,
            gamma::calculate(1.461632144968362341262).unwrap(),
            0.000000000000001,
        );
    }

    #[test]
    fn negative_real_axis_local_minima() {
        assert_within(
            -3.5446436111550050891219639933,
            gamma::calculate(-0.5040830082644554092582693045).unwrap(),
            0.00000000000001,
        );
        assert_within(
            2.3024072583396801358235820396,
            gamma::calculate(-1.5734984731623904587782860437).unwrap(),
            0.00000000000001,
        );
        assert_within(
            -0.8881363584012419200955280294,
            gamma::calculate(-2.6107208684441446500015377157).unwrap(),
            0.00000000000001,
        );
        assert_within(
            0.2451275398343662504382300889,
            gamma::calculate(-3.6352933664369010978391815669).unwrap(),
            0.00000000000001,
        );
        assert_within(
            -0.0527796395873194007604835708,
            gamma::calculate(-4.6532377617431424417145981511).unwrap(),
            0.00000000000001,
        );
        assert_within(
            0.0093245944826148505217119238,
            gamma::calculate(-5.6671624415568855358494741745).unwrap(),
            0.00000000000001,
        );
        assert_within(
            -0.0013973966089497673013074887,
            gamma::calculate(-6.6784182130734267428298558886).unwrap(),
            0.00000000000001,
        );
        assert_within(
            0.0001818784449094041881014174,
            gamma::calculate(-7.6877883250316260374400988918).unwrap(),
            0.00000000000001,
        );
        assert_within(
            -0.0000209252904465266687536973,
            gamma::calculate(-8.6957641638164012664887761608).unwrap(),
            0.00000000000001,
        );
        assert_within(
            0.0000021574161045228505405031,
            gamma::calculate(-9.7026725400018637360844267649).unwrap(),
            0.00000000000001,
        );
    }

    #[test]
    fn zero() {
        assert!(gamma::calculate(0_u32).is_none());
        assert!(gamma::calculate(-0_i32).is_none());
        assert!(gamma::calculate(0_i32).is_none());
        assert!(gamma::calculate(0.0).is_none());
        assert!(gamma::calculate(-0.0).is_none());
    }

    #[test]
    fn negative_integers_f64() {
        assert!(gamma::calculate(-1.0).is_none());
        assert!(gamma::calculate(-1.00).is_none());
        assert!(gamma::calculate(-1.000).is_none());
        assert!(gamma::calculate(-1.0000).is_none());
        assert!(gamma::calculate(-1.00000).is_none());
        assert!(gamma::calculate(-1.000000).is_none());
        assert!(gamma::calculate(-2.0).is_none());
        assert!(gamma::calculate(-3.0).is_none());
        assert!(gamma::calculate(-4.0).is_none());
        assert!(gamma::calculate(-5.0).is_none());
        assert!(gamma::calculate(-6.0).is_none());
        assert!(gamma::calculate(-7.0).is_none());
        assert!(gamma::calculate(-8.0).is_none());
        assert!(gamma::calculate(-9.0).is_none());
        assert!(gamma::calculate(-10.0).is_none());
    }

    #[test]
    fn negative_integers_i32() {
        assert!(gamma::calculate(-1).is_none());
        assert!(gamma::calculate(-2).is_none());
        assert!(gamma::calculate(-3).is_none());
        assert!(gamma::calculate(-4).is_none());
        assert!(gamma::calculate(-5).is_none());
        assert!(gamma::calculate(-6).is_none());
        assert!(gamma::calculate(-7).is_none());
        assert!(gamma::calculate(-8).is_none());
        assert!(gamma::calculate(-9).is_none());
        assert!(gamma::calculate(-10).is_none());
    }
}
