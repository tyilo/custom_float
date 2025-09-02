use core::num::FpCategory;

use crate::{fp::Fps, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[must_use]
    #[inline]
    pub fn infinity() -> Self
    {
        Fp::infinity().extra_sign(false)
    }

    #[must_use]
    #[inline]
    pub fn neg_infinity() -> Self
    {
        Fp::infinity().extra_sign(true)
    }

    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_infinite(self) -> bool
    {
        self.0.is_infinite()
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the infinite value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let infinity = FpSingle::infinity();
    ///
    /// assert!(infinity.is_infinite());
    /// assert!(!infinity.is_finite());
    /// assert!(infinity > FpSingle::max_value());
    /// ```
    #[must_use]
    #[inline]
    pub fn infinity() -> Self
    {
        Self::from_bits(Self::shift_exp(Self::max_exponent_bits()))
    }

    /// Returns the negative infinite value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let neg_infinity = FpSingle::neg_infinity();
    ///
    /// assert!(neg_infinity.is_infinite());
    /// assert!(!neg_infinity.is_finite());
    /// assert!(neg_infinity < FpSingle::min_value());
    /// ```
    #[must_use]
    #[inline]
    pub fn neg_infinity() -> Self
    {
        Self::infinity().with_sign(true)
    }

    /// Returns `true` if this value is positive infinity or negative infinity, and
    /// `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let f = FpSingle::from(7.0);
    /// let inf = FpSingle::infinity();
    /// let neg_inf = FpSingle::neg_infinity();
    /// let nan = FpSingle::nan();
    ///
    /// assert!(!f.is_infinite());
    /// assert!(!nan.is_infinite());
    ///
    /// assert!(inf.is_infinite());
    /// assert!(neg_inf.is_infinite());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_infinite(self) -> bool
    {
        matches!(self.classify(), FpCategory::Infinite)
    }
}