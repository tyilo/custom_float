use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the logarithm of the number with respect to an arbitrary base.
    ///
    /// The result might not be correctly rounded owing to implementation details;
    /// `self.log2()` can produce more accurate results for base 2, and
    /// `self.log10()` can produce more accurate results for base 10.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let ten = FpDouble::from(10.0);
    /// let two = FpDouble::from(2.0);
    ///
    /// // log10(10) - 1 == 0
    /// let abs_difference_10 = (ten.log(FpDouble::from(10.0)) - FpDouble::one()).abs();
    ///
    /// // log2(2) - 1 == 0
    /// let abs_difference_2 = (two.log(FpDouble::from(2.0)) - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference_10 < FpDouble::from(1e-10));
    /// assert!(abs_difference_2 < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn log(self, base: Self) -> Self
    {
        self.ln()/base.ln()
    }
}