use core::num::FpCategory;

use num_traits::{ConstZero, Zero};

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns `0.0`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let inf = FpSingle::infinity();
    /// let zero = FpSingle::zero();
    ///
    /// assert_eq!(FpSingle::from(7.0)/inf, zero);
    /// assert_eq!(zero * FpSingle::from(10.0), zero);
    /// ```
    #[must_use]
    pub fn zero() -> Self
    {
        Self::from_bits(U::zero())
    }
    
    // TODO: Doctest
    /// Sets `self` to the additive identity element of `Self`, `0.0`.
    fn set_zero(&mut self)
    {
        *self = Self::zero()
    }

    // TODO: Doctest
    /// Returns `true` if `self` is equal to the additive identity.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_zero(self) -> bool
    {
        matches!(self.classify(), FpCategory::Zero)
    }

    /// Returns `-0.0`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let inf = FpSingle::infinity();
    /// let zero = FpSingle::zero();
    /// let neg_zero = FpSingle::neg_zero();
    ///
    /// assert_eq!(zero, neg_zero);
    /// assert_eq!(FpSingle::from(7.0)/inf, zero);
    /// assert_eq!(zero * FpSingle::from(10.0), zero);
    /// ```
    #[must_use]
    #[inline]
    pub fn neg_zero() -> Self
    {
        Self::zero().with_sign(true)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Zero for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn zero() -> Self
    {
        Self::zero()
    }

    #[inline]
    fn is_zero(&self) -> bool
    {
        (*self).is_zero()
    }

    #[inline]
    fn set_zero(&mut self)
    {
        self.set_zero()
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> ConstZero for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> + ConstZero
{
    const ZERO: Self = Self::ZERO;
}