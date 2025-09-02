use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Raw transmutation from bits.
    ///
    /// Note that this function is distinct from [`Fp::from_uint`], which attempts to
    /// preserve the *numeric* value, and not the bitwise value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let v = FpSingle::from_bits(0x41480000);
    /// assert_eq!(v, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub const fn from_bits(bits: U) -> Self
    {
        Self(bits)
    }

    /// Raw transmutation to bits.
    ///
    /// Note that this function is distinct from [`Fp::to_uint`], which attempts to
    /// preserve the *numeric* value, and not the bitwise value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// assert_ne!(FpSingle::from(1.0).to_bits(), FpSingle::from(1.0).to_uint().unwrap()); // to_bits() is not casting!
    /// assert_eq!(FpSingle::from(12.5).to_bits(), 0x41480000);
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub const fn to_bits(self) -> U
    {
        self.0
    }

    /// Returns the sign bit of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn sign_bit(self) -> U
    {
        if !SIGN_BIT
        {
            return U::zero()
        }
        (self.to_bits() & Self::sign_mask()) >> Self::SIGN_POS
    }
    /// Returns the exponent bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn exp_bits(self) -> U
    {
        if EXP_SIZE == 0
        {
            return U::zero()
        }
        (self.to_bits() & Self::exponent_mask()) >> Self::EXP_POS
    }
    /// Returns the integer bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn int_bits(self) -> U
    {
        if Self::IS_INT_IMPLICIT
        {
            return if self.is_normal() {U::one()} else {U::zero()}
        }
        (self.to_bits() & Self::int_mask()) >> Self::INT_POS
    }
    /// Returns the fractional bits of the custom floating-point number.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn frac_bits(self) -> U
    {
        if FRAC_SIZE == 0
        {
            return U::zero()
        }
        (self.to_bits() & Self::frac_mask()) >> Self::FRAC_POS
    }

    /// Returns the exponent bias
    #[must_use]
    #[inline]
    pub fn exp_bias() -> U
    {
        Self::max_exponent_bits() >> 1usize
    }
}