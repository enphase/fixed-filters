# `fixed-filters`

Biquad filters for fixed-point numbers.

`fixed-filters` is a `#![no_std]` library for creating [biquad filters](https://en.wikipedia.org/wiki/Digital_biquad_filter).  These are simple 2nd order IIR filters which can be configured to implement a variety of filter behaviors (i.e. lowpass, bandpass, highpass, notch, proportional-integral, etc).  It is designed specifically to work with 32-bit fixed-point numbers based on the [fixed](https://crates.io/crates/fixed) crate.

**Alpha:** This crate requires the alpha release of 2.0.0 of the [fixed](https://crates.io/crates/fixed) crate that makes use of const generics instead of the [*typenum*
crate](https://crates.io/crate/typenum). This version requires the nightly compiler with the [`generic_const_exprs` feature] enabled.

# How to use

The main data structure for this crate is `Biquad<X,Y>` where `X` and `Y` are `i32` constant generics which represent the number of fractional bits of the input and output signals respectively.  While the run functions are all designed use fixed-point arithmetic, the construction APIs accept floating point variables to simplify use.  The most generic way of creating a filter is to construct with `new()` and pass the filter coefficients.

```rust
#![feature(generic_const_exprs)]

use fixed_filters::Biquad;
use fixed::FixedI32;

// poles
let a0 = 1.;
let a1 = -1.97985135;
let a2 = 0.980052318;

// zeros
let b0 = 0.000050241;
let b1 = 0.000100482;
let b2 = 0.000050241;

// make a filter with 20 fractional bits on input and 16 fractional bits on output
let mut filter = Biquad::<20, 16>::new(a0, a1, a2, b0, b1, b2);

let signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
let mut y = FixedI32::<16>::ZERO;
for s in signal {
	let x = FixedI32::<20>::from_num(s); // convert to fixed-point number
	y = filter.update(x);
}
```

While passing the pole and zero coefficients is the most univeral way of creating a filter, `Biquad` includes a number of associated functions for creating standard filter forms with the use specifying just the sampling period, the center frequency, and the q factor.

```rust
// create a butterworth lowpass filter
let ts = 1./10e3;
let f0 = 100.;
let q = core::f32::consts::FRAC_1_SQRT_2;
let mut filter = Biquad::<16, 16>::lowpass(ts, f0, q);
```

The library includes associted functions for creating the following filter forms:

1.  lowpass
2.  bandpass
3.  highpass
4.  noth
5.  single_pole_lowpass
6.  proportional-integral behavior

By default, the filters do not limit their output and can in fact over-flow if you are not careful.  You can add limits to an existing filter using the `set_min`, `set_max`, or `set_limit` functions, or you can setup limits at construction using a builder pattern with the `with_min`, `with_max`, or `with_limit` functions.

```rust
// create a proportional-integral controller with a limit ot +-20
let kp = 0.1;
let ki = 0.01;
let mut filter = Biquad::<16, 16>::pi(kp, ki).with_limit(20.0);
```

# Benchmarking

The performance of this crate was benchmarked against other available `#![no_std]` biquad filter implementations by running a filter at 10kHz in a simple [RTIC](https://rtic.rs/1/book/en/) application on an STM32F411RE Nucelo board.  The RTIC monotonic API was used to quantify the number of ticks used for each filters' update method.


| Crate                         | Cycles  |
| ----------------------------- | -------:|
| fixed-filters                 | 90      |
| biquad::DirectForm1           | 70      |
| biquad::DirectForm2Transposed | 64      |
| idsp::iir_int                 |         |