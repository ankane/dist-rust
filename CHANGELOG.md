## 0.2.1 (unreleased)

- Added `no_std` feature

## 0.2.0 (2022-08-28)

- Improved accuracy of `Normal::ppf` and `StudentsT::ppf`
- Return `NAN` instead of panicking for invalid inputs

## 0.1.3 (2022-07-31)

- Added support for `df` between zero and one to `StudentsT::pdf`
- Fixed error with `StudentsT::cdf` when `x` is infinite or NaN
- Fixed bug with `StudentsT` functions when `df` is infinity

## 0.1.2 (2022-07-26)

- Improved accuracy of `Normal::cdf` and `StudentsT::cdf`
- Fixed bug with `Normal` functions when `std_dev` is not one

## 0.1.1 (2022-04-20)

- Added support for non-integer degrees of freedom

## 0.1.0 (2022-01-03)

- First release
