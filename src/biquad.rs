#[allow(unused_imports)]
use micromath::F32Ext; // enable floating point sin/cos approximations in microcontrollers

use core::f32::consts::PI;
use fixed::types::I2F30;
use fixed::FixedI32;
use heapless::HistoryBuffer;

type Coef = I2F30;

// normalize coefficients to 1/a0
fn normalize(a0: f32, a1: f32, a2: f32, b0: f32, b1: f32, b2: f32) -> (f32, f32, f32, f32, f32) {
    let x = a0.recip();
    (a1 * x, a2 * x, b0 * x, b1 * x, b2 * x)
}

fn calculate_intermediate_variables(ts: f32, f0: f32, q: f32) -> (f32, f32) {
    let omega = 2. * PI * f0 * ts;
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();
    let alpha = sin_omega / (2. * q);
    (alpha, cos_omega)
}

/// A fixed-point digital [Biquad filter](https://en.wikipedia.org/wiki/Digital_biquad_filter)
///
/// This is a Direct Form 1 fixed point implementation that uses the [fixed crate](https://crates.io/crates/fixed).
/// The filter works on fixed-point 32 bit signals where the fractional bits is passed as a constant generic
/// `Biquad<16, 20>` will create a 32-bit biquad filter for a 32-bit signal with 16 fractional bits on the
/// input and 20 fractional bits on the output.
///
/// While the run methods for the filter use fixed-point arithmetic, all of the generating functions use
/// floating-point arithmetic and inputs to simplify generation of filters.
///
/// While you can create Biqad filters by passing raw coefficients using the `Biquad::new()` method, it is
/// there are also a number of associated functions for generating Biquads of different standard types.
///
/// # Examples
///
/// ```
/// use fixed_filters::Biquad;
/// use fixed::FixedI32;
///
/// let ts = 1./10e3;  // sampling period
/// let f0 = 15.0;  // cutoff frequency
/// let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
/// let mut filter = Biquad::<20, 20>::lowpass(ts, f0, q);
///
/// let signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
/// let mut y = FixedI32::<20>::ZERO;
/// for s in signal {
///    let x = FixedI32::<20>::from_num(s); // convert to fixed-point number
///    y = filter.update(x);
/// }
/// ```
pub struct Biquad<const X: i32, const Y: i32> {
    // coefficients
    a1: Coef,
    a2: Coef,
    b0: Coef,
    b1: Coef,
    b2: Coef,

    // history
    x: HistoryBuffer<FixedI32<X>, 2>,
    y: HistoryBuffer<FixedI32<Y>, 2>,
}

impl<const X: i32, const Y: i32> Biquad<X, Y> {
    /// Create a new biquad filter from the biquad coefficients
    ///
    /// While this method can be used to create any biqauad filter, most use cases
    /// should use one of the associated methods which will automatically
    /// calculate the coefficients for the standard filter forms from the
    /// sampling frequency, center frequency, and the quality factor.
    pub fn new(a0: f32, a1: f32, a2: f32, b0: f32, b1: f32, b2: f32) -> Self {
        let (a1, a2, b0, b1, b2) = normalize(a0, a1, a2, b0, b1, b2);
        Self {
            a1: Coef::from_num(a1),
            a2: Coef::from_num(a2),
            b0: Coef::from_num(b0),
            b1: Coef::from_num(b1),
            b2: Coef::from_num(b2),
            x: HistoryBuffer::new_with(FixedI32::<X>::ZERO),
            y: HistoryBuffer::new_with(FixedI32::<Y>::ZERO),
        }
    }

    /// Add a new input sample and get the resulting output
    pub fn update(&mut self, x: FixedI32<X>) -> FixedI32<Y> {
        // Calculate the biquad output using Direct Form 1
        // TODO: Consider using a double wide accumulator for better precision
        let mut y = FixedI32::<Y>::ZERO;
        y.mul_acc(self.b0, x);
        y.mul_acc(self.b1, *self.x.recent().unwrap());
        y.mul_acc(self.b2, *self.x.oldest().unwrap());
        y.mul_acc(-self.a1, *self.y.recent().unwrap());
        y.mul_acc(-self.a2, *self.y.oldest().unwrap());

        // update the FIFOs
        self.x.write(x);
        self.y.write(y);

        y
    }

    /// Reset the filter
    ///
    /// This will clear the fifos in the biquad filter
    pub fn reset(&mut self) {
        self.x.clear_with(FixedI32::<X>::ZERO);
        self.y.clear_with(FixedI32::<Y>::ZERO);
    }

    /// Return the latest output of the filter without updating
    pub fn value(&self) -> FixedI32<Y> {
        *self.y.recent().unwrap()
    }

    /// Constructs a biquad filter with single pole lowpass (i.e.  "RC") behavior
    ///
    /// Effectively makes `y[n] = y[n-1] + alpha * (x[n]-y[n-1])` style filter
    ///
    /// # Examples
    ///
    /// ```
    /// use fixed_filters::Biquad;
    ///
    /// let ts = 1./10e3;  // sampling period
    /// let f0 = 15.0;  // cutoff frequency
    /// let mut filter = Biquad::<25, 25>::single_pole_lowpass(ts, f0);
    /// ```
    pub fn single_pole_lowpass(ts: f32, f0: f32) -> Self {
        let omega = 2. * 3.14159 * f0 * ts;
        let alpha = omega / (omega + 1.);
        let a0 = 1.;
        let b0 = alpha;
        let b1 = 0.;
        let b2 = 0.;
        let a1 = alpha - 1.;
        let a2 = 0.;
        Self::new(a0, a1, a2, b0, b1, b2)
    }

    /// Constructs a lowpass biquad filter
    ///
    /// Uses arithmetic from <https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html>
    ///
    /// # Examples
    ///
    /// ```
    /// use fixed_filters::Biquad;
    ///
    /// let ts = 1./10e3;  // sampling period
    /// let f0 = 15.0;  // cutoff frequency
    /// let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
    /// let mut filter = Biquad::<20, 20>::lowpass(ts, f0, q);
    /// ```
    pub fn lowpass(ts: f32, f0: f32, q: f32) -> Self {
        let (alpha, cos_omega) = calculate_intermediate_variables(ts, f0, q);
        let a0 = 1. + alpha;
        let a1 = -2. * cos_omega;
        let a2 = 1. - alpha;
        let b0 = (1. - cos_omega) / 2.;
        let b1 = 2. * b0;
        let b2 = b0;
        Self::new(a0, a1, a2, b0, b1, b2)
    }

    /// Constructs a highpass biquad filter
    ///
    /// Uses arithmetic from <https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html>
    ///
    /// # Examples
    ///
    /// ```
    /// use fixed_filters::Biquad;
    ///
    /// let ts = 1./10e3;  // sampling period
    /// let f0 = 15.0;  // cutoff frequency
    /// let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
    /// let mut filter = Biquad::<20, 20>::highpass(ts, f0, q);
    /// ```
    pub fn highpass(ts: f32, f0: f32, q: f32) -> Self {
        let (alpha, cos_omega) = calculate_intermediate_variables(ts, f0, q);
        let a0 = 1. + alpha;
        let a1 = -2. * cos_omega;
        let a2 = 1. - alpha;
        let b0 = (1. + cos_omega) / 2.;
        let b1 = -2. * b0;
        let b2 = b0;
        Self::new(a0, a1, a2, b0, b1, b2)
    }

    /// Constructs a bandpass biquad filter
    ///
    /// Uses arithmetic from <https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html>
    ///
    /// # Examples
    ///
    /// ```
    /// use fixed_filters::Biquad;
    ///
    /// let ts = 1./10e3;  // sampling period
    /// let f0 = 15.0;  // center frequency
    /// let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
    /// let mut filter = Biquad::<20, 20>::bandpass(ts, f0, q);
    /// ```
    pub fn bandpass(ts: f32, f0: f32, q: f32) -> Self {
        let (alpha, cos_omega) = calculate_intermediate_variables(ts, f0, q);
        let a0 = 1. + alpha;
        let a1 = -2. * cos_omega;
        let a2 = 1. - alpha;
        let b0 = alpha;
        let b1 = 0.;
        let b2 = -b0;
        Self::new(a0, a1, a2, b0, b1, b2)
    }

    /// Constructs a notch biquad filter
    ///
    /// Uses arithmetic from <https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html>
    ///
    /// # Examples
    ///
    /// ```
    /// use fixed_filters::Biquad;
    ///
    /// let ts = 1./10e3;  // sampling period
    /// let f0 = 15.0;  // notch frequency
    /// let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
    /// let mut filter = Biquad::<20, 20>::notch(ts, f0, q);
    /// ```
    pub fn notch(ts: f32, f0: f32, q: f32) -> Self {
        let (alpha, cos_omega) = calculate_intermediate_variables(ts, f0, q);
        let a0 = 1. + alpha;
        let a1 = -2. * cos_omega;
        let a2 = 1. - alpha;
        let b0 = 1.;
        let b1 = a1;
        let b2 = 1.;
        Self::new(a0, a1, a2, b0, b1, b2)
    }

    /// Constructs a biquad filter with proportional-integral behavior
    ///
    /// # Examples
    ///
    /// ```
    /// use fixed_filters::Biquad;
    ///
    /// let kp = 0.2;
    /// let ki = 0.1;
    /// let mut filter = Biquad::<25, 30>::pi(kp, ki);
    /// ```
    pub fn pi(kp: f32, ki: f32) -> Self {
        let a0 = 1.;
        let a1 = -1.;
        let a2 = 0.;
        let b0 = kp + ki;
        let b1 = -kp;
        let b2 = 0.;
        Self::new(a0, a1, a2, b0, b1, b2)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use approx::assert_abs_diff_eq;

    // simulation function
    fn simulate<const FRAC: i32>(
        f: f32,
        fs: f32,
        filter: &mut Biquad<FRAC, FRAC>,
    ) -> FixedI32<FRAC> {
        // simulate
        let mut t = 0.0;
        let dt = 1. / fs;
        let stop = 10. / f;
        let mut y = FixedI32::<FRAC>::ZERO;
        let mut max_y = y;
        while t < stop {
            t += dt;
            let x = (2.0 * PI * f * t).sin();
            y = filter.update(FixedI32::<FRAC>::from_num(x));

            // allow the filter to converge before we record the amplitude
            if t > 0.9 * stop && y > max_y {
                max_y = y
            };
        }
        max_y
    }

    #[test]
    fn fifo() {
        let mut x = HistoryBuffer::<FixedI32<16>, 2>::new_with(FixedI32::<16>::ZERO);
        x.write(FixedI32::<16>::from_num(0.2));
        x.write(FixedI32::<16>::from_num(0.4));
        assert_eq!(*x.recent().unwrap(), FixedI32::<16>::from_num(0.4));
        assert_eq!(*x.oldest().unwrap(), FixedI32::<16>::from_num(0.2));
    }

    #[test]
    fn pi() {
        let kp = 0.2;
        let ki = 0.1;
        let mut filter = Biquad::<25, 25>::pi(kp, ki);

        let mut y = filter.update(FixedI32::<25>::from_num(1.0));
        assert_abs_diff_eq!(f64::from(y), 0.3, epsilon = 0.0001);
        y = filter.update(FixedI32::<25>::from_num(1.0));
        assert_abs_diff_eq!(f64::from(y), 0.4, epsilon = 0.0001);
    }

    #[test]
    fn lowpass() {
        let fs = 10e3;
        let f0 = 60.0;
        let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
        let mut filter = Biquad::<25, 25>::lowpass(1. / fs, f0, q);

        // simulate
        filter.reset();
        let amplitude = simulate(0.5 * f0, fs, &mut filter);
        assert!(amplitude > 0.95);

        filter.reset();
        let amplitude = simulate(2.0 * f0, fs, &mut filter);
        assert!(amplitude < 0.5);
    }

    #[test]
    fn highpas() {
        let fs = 10e3;
        let f0 = 60.0;
        let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
        let mut filter = Biquad::<25, 25>::highpass(1. / fs, f0, q);

        // simulate
        filter.reset();
        let amplitude = simulate(0.5 * f0, fs, &mut filter);
        assert!(amplitude < 0.5);

        filter.reset();
        let amplitude = simulate(2.0 * f0, fs, &mut filter);
        assert!(amplitude > 0.95);
    }

    #[test]
    fn bandpass() {
        let fs = 10e3;
        let f0 = 60.0;
        let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
        let mut filter = Biquad::<25, 25>::bandpass(1. / fs, f0, q);

        // simulate
        filter.reset();
        let amplitude = simulate(0.5 * f0, fs, &mut filter);
        assert!(amplitude < 0.75);

        filter.reset();
        let amplitude = simulate(f0, fs, &mut filter);
        assert!(amplitude > 0.95);

        filter.reset();
        let amplitude = simulate(2.0 * f0, fs, &mut filter);
        assert!(amplitude < 0.75);
    }

    #[test]
    fn notch() {
        let fs = 1e3;
        let f0 = 60.0;
        let q = core::f32::consts::FRAC_1_SQRT_2; // butterworth
        let mut filter = Biquad::<25, 25>::notch(1. / fs, f0, q);

        // simulate
        filter.reset();
        let amplitude = simulate(0.5 * f0, fs, &mut filter);
        assert!(amplitude > 0.7);

        filter.reset();
        let amplitude = simulate(f0, fs, &mut filter);
        assert!(amplitude < 0.1);

        filter.reset();
        let amplitude = simulate(2.0 * f0, fs, &mut filter);
        assert!(amplitude > 0.7);
    }
}