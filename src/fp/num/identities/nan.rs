use core::num::FpCategory;

use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the `NaN` value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let nan = FpSingle::nan();
    ///
    /// assert!(nan.is_nan());
    /// ```
    #[must_use]
    #[inline]
    pub fn nan() -> Self
    {
        Self::snan()
    }

    /// Returns `true` if this value is NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let nan = FpDouble::nan();
    /// let f = FpDouble::from(7.0);
    ///
    /// assert!(nan.is_nan());
    /// assert!(!f.is_nan());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_nan(self) -> bool
    {
        matches!(self.classify(), FpCategory::Nan)
    }
    
    /// Returns the `qNaN` value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let qnan = FpSingle::qnan();
    ///
    /// assert!(qnan.is_nan());
    /// assert!(!qnan.is_snan());
    /// ```
    #[must_use]
    #[inline]
    pub fn qnan() -> Self
    {
        Self::from_bits(U::max_value() >> (util::bitsize_of::<U>() - Self::SIGN_POS))
    }
    
    /// Returns the `sNaN` value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let snan = FpSingle::snan();
    ///
    /// assert!(snan.is_nan());
    /// assert!(snan.is_snan());
    /// ```
    #[must_use]
    #[inline]
    pub fn snan() -> Self
    {
        let mut nan = Self::qnan();
        if Self::EXP_POS != 0
        {
            let mask = !(U::one() << (Self::EXP_POS - 1));
            nan = Self::from_bits(nan.to_bits() & mask);
        }
        nan
    }

    /// Returns `true` if the number is a signaling NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let snan = FpSingle::snan();
    /// let qnan = FpSingle::qnan();
    ///
    /// assert!(snan.is_snan());
    /// assert!(!qnan.is_snan());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_snan(self) -> bool
    {
        Self::MANTISSA_DIGITS != 0 && self.is_nan() && {
            let mask = U::one() << (Self::MANTISSA_DIGITS - 1);
            (self.to_bits() & mask).is_zero()
        }
    }
}