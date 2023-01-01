//! PDF, CDF, and percent-point/quantile functions for the normal and Student’s t distributions
//!
//! [View the docs](https://github.com/ankane/dist-rust)

#![cfg_attr(feature = "libm", forbid(unsafe_code))]
#![cfg_attr(not(feature = "libm"), deny(unsafe_code))]

mod normal;
mod students_t;

#[cfg(not(feature = "libm"))]
mod math;

pub use normal::Normal;
pub use students_t::StudentsT;
