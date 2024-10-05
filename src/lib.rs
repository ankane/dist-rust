#![doc = include_str!("../README.md")]
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
