use num_traits::FromBytes;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> + FromBytes
{
    /// Create a floating point value from its representation as a byte array in big endian.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let value = FpSingle::from_be_bytes(&[0x41, 0x48, 0x00, 0x00]);
    /// assert_eq!(value, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn from_be_bytes(bytes: &U::Bytes) -> Self
    {
        Self::from_bits(U::from_be_bytes(bytes))
    }

    /// Create a floating point value from its representation as a byte array in little endian.
    ///
    /// See [`from_bits`](Self::from_bits) for some discussion of the
    /// portability of this operation (there are almost no issues).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let value = FpSingle::from_le_bytes(&[0x00, 0x00, 0x48, 0x41]);
    /// assert_eq!(value, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn from_le_bytes(bytes: &U::Bytes) -> Self
    {
        Self::from_bits(U::from_le_bytes(bytes))
    }

    /// Create a floating point value from its representation as a byte array in native endian.
    ///
    /// As the target platform's native endianness is used, portable code
    /// likely wants to use [`from_be_bytes`] or [`from_le_bytes`], as
    /// appropriate instead.
    ///
    /// [`from_be_bytes`]: f32::from_be_bytes
    /// [`from_le_bytes`]: f32::from_le_bytes
    ///
    /// See [`from_bits`](Self::from_bits) for some discussion of the
    /// portability of this operation (there are almost no issues).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let value = FpSingle::from_ne_bytes(if cfg!(target_endian = "big") {
    ///     &[0x41, 0x48, 0x00, 0x00]
    /// } else {
    ///     &[0x00, 0x00, 0x48, 0x41]
    /// });
    /// assert_eq!(value, FpSingle::from(12.5));
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn from_ne_bytes(bytes: &U::Bytes) -> Self
    {
        Self::from_bits(U::from_ne_bytes(bytes))
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> FromBytes for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> + FromBytes
{
    type Bytes = U::Bytes;

    #[inline]
    fn from_be_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_be_bytes(bytes)
    }

    #[inline]
    fn from_le_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_le_bytes(bytes)
    }

    #[inline]
    fn from_ne_bytes(bytes: &Self::Bytes) -> Self
    {
        Self::from_ne_bytes(bytes)
    }
}