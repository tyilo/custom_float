use num_traits::FloatConst;

use crate::{fp::NEWTON_LN, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn log10_nonewton(self) -> Self
    {
        if EXP_BASE == 10
        {
            return self.logb()
        }
        self.ln_nonewton()/Self::LN_10()
    }

    /// Returns the base 10 logarithm of the number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let ten = FpDouble::from(10.0);
    ///
    /// // log10(10) - 1 == 0
    /// let abs_difference = (ten.log10() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-4));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn log10(self) -> Self
    {
        let mut y = self.log10_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            let one = Self::one();
            let ln10 = Self::LN_10();
            for _ in 0..NEWTON
            {
                let xdyexp10 = self/y.exp10();
                let (diff, s) = one.sub_extra_sign(xdyexp10);
                let dy = diff/ln10;
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