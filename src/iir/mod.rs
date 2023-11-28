pub mod biquad16;
pub mod biquad32;

pub use biquad16::Biquad16;
pub use biquad32::Biquad32;
use core::f32::consts::PI;
use fixed::types::I2F30;

#[allow(unused_imports)]
use micromath::F32Ext; // enable floating point sin/cos approximations in microcontrollers

pub type Coef = I2F30;

/// Normalize coefficients to 1/a0
pub fn normalize(
    a0: f32,
    a1: f32,
    a2: f32,
    b0: f32,
    b1: f32,
    b2: f32,
) -> (f32, f32, f32, f32, f32) {
    let x = a0.recip();
    (a1 * x, a2 * x, b0 * x, b1 * x, b2 * x)
}

/// Calculate intermediate variables needed for various filter coefficients
pub fn calculate_intermediate_variables(ts: f32, f0: f32, q: f32) -> (f32, f32) {
    let omega = 2. * PI * f0 * ts;
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();
    let alpha = sin_omega / (2. * q);
    (alpha, cos_omega)
}
