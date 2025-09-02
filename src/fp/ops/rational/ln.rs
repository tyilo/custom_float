use crate::{fp::NEWTON_LN, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    pub(crate) fn ln_nonewton(self) -> Self
    {
        self.logb()*Self::ln_exp_base()
    }

    #[inline]
    pub(crate) fn ln_nonewton_extra_sign(mut self) -> (Self, bool)
    {
        //TODO: Is-less-than-one function?
        let s = !SIGN_BIT && self.abs() < Self::one();
        if s
        {
            self = self.recip()
        }
        (self.ln_nonewton(), s)
    }

    /// Returns the natural logarithm of the number.
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
    pub fn ln(self) -> Self
    {
        let mut y = self.ln_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            let one = Self::one();
            for _ in 0..NEWTON
            {
                let x = self/y.exp();
                let (dy, s) = one.add_with_sign_extra_sign(false, x, true);
                let yy = y.add_with_sign(false, dy, !s);
                if !yy.is_finite()
                {
                    break
                }
                y = yy;
            }
        }

        y
    }
}