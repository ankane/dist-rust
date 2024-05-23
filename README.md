# Dist Rust

PDF, CDF, and percent-point/quantile functions for the normal and Student’s t distributions

🎉 Zero dependencies

[![Build Status](https://github.com/ankane/dist-rust/actions/workflows/build.yml/badge.svg)](https://github.com/ankane/dist-rust/actions)

## Installation

Add this line to your application’s `Cargo.toml` under `[dependencies]`:

```toml
distrs = "0.2"
```

## Getting Started

- [Normal](#normal)
- [Student’s t](#students-t)

### Normal

```rust
use distrs::Normal;

Normal::pdf(x, mean, std_dev);
Normal::cdf(x, mean, std_dev);
Normal::ppf(p, mean, std_dev);
```

### Student’s t

```rust
use distrs::StudentsT;

StudentsT::pdf(x, df);
StudentsT::cdf(x, df);
StudentsT::ppf(p, df);
```

## Features

- `no_std` - enable `no_std` support (requires [libm](https://github.com/rust-lang/libm))

## References

- [Algorithm AS 241: The Percentage Points of the Normal Distribution](https://www.jstor.org/stable/2347330)
- [Algorithm 395: Student’s t-distribution](https://dl.acm.org/doi/10.1145/355598.355599)
- [Algorithm 396: Student’s t-quantiles](https://dl.acm.org/doi/10.1145/355598.355600)

## History

View the [changelog](https://github.com/ankane/dist-rust/blob/master/CHANGELOG.md)

## Contributing

Everyone is encouraged to help improve this project. Here are a few ways you can help:

- [Report bugs](https://github.com/ankane/dist-rust/issues)
- Fix bugs and [submit pull requests](https://github.com/ankane/dist-rust/pulls)
- Write, clarify, or fix documentation
- Suggest or add new features

To get started with development:

```sh
git clone https://github.com/ankane/dist-rust.git
cd dist-rust
cargo test
```
