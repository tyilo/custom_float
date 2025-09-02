use crate::{fp::NEWTON_EXP, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    pub(super) fn exp_nonewton(self) -> Self
    {
        (self/Self::ln_exp_base()).expb()
    }

    /// Returns `e^(self)`, (the exponential function).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let one = FpDouble::one();
    /// // e^1
    /// let e = one.exp();
    /// 
    /// // ln(e) - 1 == 0
    /// let abs_difference = (e.ln() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn exp(self) -> Self
    {
        let mut y = self.exp_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_EXP;

            for _ in 0..NEWTON
            {
                let (x, x_s) = y.ln_nonewton_extra_sign();
                let (dxx, s) = x.add_with_sign_extra_sign(x_s, self, true);
                let dy = y*dxx;
                let yy = y.add_with_sign(false, dy, !s);

                if !yy.is_finite() || yy.is_sign_negative()
                {
                    break
                }
                y = yy;
            }
        }

        y
    }
}