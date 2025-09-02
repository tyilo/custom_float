use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the largest finite value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::max_value();
    /// assert_eq!(x, FpDouble::from(f64::MAX));
    /// ```
    #[must_use]
    #[inline]
    pub fn max_value() -> Self
    {
        Self::from_bits(Self::shift_exp(Self::max_exponent_bits()) - U::one())
    }
    
    /// Returns the smallest finite value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::min_value();
    ///
    /// assert_eq!(x, FpDouble::from(f64::MIN));
    /// ```
    #[must_use]
    #[inline]
    pub fn min_value() -> Self
    {
        if !SIGN_BIT
        {
            return Self::zero()
        }
        Self::max_value().with_sign(true)
    }

    /// Returns the smallest positive, normal value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::f64;
    ///
    /// let x = FpDouble::min_positive_value();
    ///
    /// assert_eq!(x, FpDouble::from(f64::MIN_POSITIVE));
    /// ```
    #[must_use]
    #[inline]
    pub fn min_positive_value() -> Self
    {
        let mut bits = Self::shift_exp(U::one());
        if !Self::IS_INT_IMPLICIT
        {
            bits = bits + Self::shift_int(U::one());
        }
        Self::from_bits(bits)
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
        let bias = Self::exp_bias();
        
        if !Self::IS_INT_IMPLICIT
        {
            return Self::from_bits(Self::shift_exp(bias) + Self::shift_frac(U::one()))
        }

        let exp_frac = U::from(util::exp2_ilog(FRAC_SIZE + INT_SIZE, EXP_BASE)).unwrap();
        if bias <= exp_frac
        {
            return Self::from_bits(util::powu::<U, usize>(U::from(EXP_BASE).unwrap(), bias.to_usize().unwrap() - 1))
        }
        Self::from_bits(Self::shift_exp(bias - exp_frac))
    }
}