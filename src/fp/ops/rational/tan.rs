use core::cmp::Ordering;

use num_traits::FloatConst;

use crate::{fp::NEWTON_TRIG, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Computes the tangent of a number (in radians).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::FRAC_PI_4();
    /// 
    /// let abs_difference = (x.tan() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn tan(mut self) -> Self
    {
        let pi = Self::PI();
        let half_pi = Self::FRAC_PI_2();

        self %= pi;
        if matches!(self.abs_partial_cmp(half_pi), Some(Ordering::Greater))
        {
            self = self.add_with_sign(false, pi, self.is_sign_positive())
        }

        let (sin, ss, cos, cs) = self.sin_cos_extra_sign();

        let mut s = ss ^ cs;
        let mut y = sin/cos;

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let x = y.atan_with_sign_extra_sign();
                let (dx, s) = x.sub_extra_sign(self);
                let dy = dx*(y.squared() + Self::one());
                let yy = y.add_with_sign(false, dy, !s);
                if !yy.is_finite()
                {
                    break
                }
                y = yy
            }
        }
        if y.is_nan() // Checks again
        {
            return Self::infinity().copysign(self)
        }

        y
    }
}