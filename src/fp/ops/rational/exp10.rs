use num_traits::FloatConst;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns `10^(self)`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(2.0);
    ///
    /// // 10^2 - 100 == 0
    /// let abs_difference = (f.exp10() - FpDouble::from(100.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-6));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn exp10(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.expb()
        }
        (self*Self::LN_10()).exp()
    }
}