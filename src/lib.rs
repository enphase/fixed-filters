// Copyright 2023 Enphase Energy, Inc.
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

/*!
# `fixed-filters`

Biquad filters for fixed-point numbers.

`fixed-filters` is a `#![no_std]` library for creating [biquad filters](https://en.wikipedia.org/wiki/Digital_biquad_filter).  These are simple 2nd order IIR filters which can be configured to implement a variety of filter behaviors (i.e. lowpass, bandpass, highpass, notch, proportional-integral, etc).  It is designed specifically to work with 32-bit fixed-point numbers based on the [fixed](https://crates.io/crates/fixed) crate.

**Alpha:** This crate requires the alpha release of 2.0.0 of the [fixed](https://crates.io/crates/fixed) crate that makes use of const generics instead of the [*typenum* crate](https://crates.io/crate/typenum). This version requires the nightly compiler with the [`generic_const_exprs` feature] enabled.

# How to use

The main data structure for this crate is [Biquad32<X,Y>](crate::biquad::Biquad32) where `X` and `Y` are `i32` constant generics which represent the number of fractional bits of the input and output signals respectively.  While the run functions are all designed use fixed-point arithmetic, the construction APIs accept floating point variables to simplify use.  The most generic way of creating a filter is to construct with `new()` and pass the filter coefficients.

```rust
#![feature(generic_const_exprs)]

use fixed_filters::Biquad32;
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
let mut filter = Biquad32::<20, 16>::new(a0, a1, a2, b0, b1, b2);

let signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
let mut y = FixedI32::<16>::ZERO;
for s in signal {
    let x = FixedI32::<20>::from_num(s); // convert to fixed-point number
    y = filter.update(x);
}
```

While passing the pole and zero coefficients is the most universal way of creating a filter, `Biquad32` includes a number of associated functions for creating standard filter forms with the user specifying just the sampling period, the center frequency, and the q factor.

```rust
use fixed_filters::Biquad32;

// create a butterworth lowpass filter
let ts = 1./10e3;
let f0 = 100.;
let q = core::f32::consts::FRAC_1_SQRT_2;
let mut filter = Biquad32::<16, 16>::lowpass(ts, f0, q);
```

The library includes associted functions for creating the following filter forms:

1.  [lowpass](crate::biquad::Biquad32::lowpass)
2.  [bandpass](crate::biquad::Biquad32::bandpass)
3.  [highpass](crate::biquad::Biquad32::highpass)
4.  [notch](crate::biquad::Biquad32::notch)
5.  [single_pole_lowpass](crate::biquad::Biquad32::single_pole_lowpass)
6.  [proportional-integral behavior](crate::biquad::Biquad32::pi)

By default, the filters do not limit their output and can in fact over-flow if you are not careful.  You can add limits to an existing filter using the `set_min`, `set_max`, or `set_limit` functions, or you can setup limits at construction using a builder pattern with the `with_min`, `with_max`, or `with_limit` functions.

```rust
use fixed_filters::Biquad32;

// create a proportional-integral controller with a limit ot +-20
let kp = 0.1;
let ki = 0.01;
let mut filter = Biquad32::<16, 16>::pi(kp, ki).with_limit(20.0);
```

# Benchmarking

The performance of this crate was benchmarked against other available `#![no_std]` biquad filter implementations by running a filter at 10kHz in a simple [RTIC](https://rtic.rs/1/book/en/) application on an [STM32F411RE Nucelo](https://www.st.com/en/evaluation-tools/nucleo-f411re.html) board.  The RTIC monotonic API was used to quantify the number of CPU cycles used for each filters' update method.

| Crate                         | Cycles  |
| :---------------------------- |:-------:|
| fixed-filters                 | 80      |
| biquad::DirectForm1           | 80      |
| biquad::DirectForm2Transposed | 72      |
| idsp::iir_int                 | 90      |

From these results, it can be seen that the crates all provide similar performance in terms of execution efficiency with [biquad](https://crates.io/crates/biquad/0.3.0) DirectForm2Transposed having a slight edge in microcontrollers with an FPU.  This crate however does not support internal limiting of the output signal which can be critical in certain applications such as implementing anti-windup and PI controllers.
*/
#![cfg_attr(not(test), no_std)]

pub mod iir;

pub use iir::{Biquad16, Biquad32};
