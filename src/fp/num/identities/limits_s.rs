use crate::{fp::Fps, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[must_use]
    #[inline]
    pub fn max_value() -> Self
    {
        Fp::max_value().extra_sign(false)
    }
    
    #[must_use]
    #[inline]
    pub fn min_value() -> Self
    {
        Fp::max_value().extra_sign(true)
    }
    
    #[must_use]
    #[inline]
    pub fn min_positive_value() -> Self
    {
        Fp::min_positive_value().extra_sign(false)
    }

    /// [Machine epsilon] value.
    ///
    /// This is the difference between `1.0` and the next larger representable number.
    ///
    /// [Machine epsilon]: https://en.wikipedia.org/wiki/Machine_epsilon
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::epsilon();
    ///
    /// assert_eq!(x, FpDouble::from(f64::EPSILON));
    /// ```
    #[must_use]
    #[inline]
    pub fn epsilon() -> Self
    {
        Fp::min_positive_value().extra_sign(false)
    }
}