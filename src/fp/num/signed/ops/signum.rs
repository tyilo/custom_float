use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns a number that represents the sign of `self`.
    ///
    /// - `1.0` if the number is positive, `+0.0` or `inf`
    /// - `-1.0` if the number is negative, `-0.0` or `-inf`
    /// - `NaN` if the number is `NaN`
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.5);
    ///
    /// assert_eq!(f.signum(), FpDouble::one());
    /// assert_eq!(FpDouble::neg_infinity().signum(), -FpDouble::one());
    ///
    /// assert!(FpDouble::nan().signum().is_nan());
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn signum(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        Self::one().copysign(self)
    }
}