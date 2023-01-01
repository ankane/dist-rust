//! PDF, CDF, and percent-point/quantile functions for the normal and Student’s t distributions
//!
//! [View the docs](https://github.com/ankane/dist-rust)

#![cfg_attr(feature = "no_std", no_std)]
#![cfg_attr(feature = "no_std", forbid(unsafe_code))]
#![cfg_attr(not(feature = "no_std"), deny(unsafe_code))]

mod normal;
mod students_t;

#[cfg(feature = "no_std")]
use libm as math;

#[cfg(not(feature = "no_std"))]
mod math;

pub use normal::Normal;
pub use students_t::StudentsT;
