//! PDF, CDF, and percent-point/quantile functions for the normal and Student’s t distributions
//!
//! [View the docs](https://github.com/ankane/dist-rust)

#![no_std]

mod erf;
mod normal;
mod students_t;

pub use normal::Normal;
pub use students_t::StudentsT;
