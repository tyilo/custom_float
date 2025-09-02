use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// The positive difference of two numbers.
    ///
    /// * If `self <= other`: `0:0`
    /// * Else: `self - other`
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(3.0);
    /// let y = FpDouble::from(-3.0);
    ///
    /// let abs_difference_x = (x.abs_sub(FpDouble::one()) - FpDouble::from(2.0)).abs();
    /// let abs_difference_y = (y.abs_sub(FpDouble::one()) - FpDouble::zero()).abs();
    ///
    /// assert!(abs_difference_x < FpDouble::from(1e-10));
    /// assert!(abs_difference_y < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[deprecated(
        note = "you probably meant `(self - other).abs()`: \
                this operation is `(self - other).max(0.0)` \
                except that `abs_sub` also propagates NaNs (also \
                known as `fdimf` in C)."
    )]
    #[inline]
    pub fn abs_sub(self, other: Self) -> Self
    {
        if other >= self
        {
            return Self::zero()
        }
        self - other
    }
}