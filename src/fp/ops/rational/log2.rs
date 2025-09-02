use num_traits::FloatConst;

use crate::{fp::NEWTON_LN, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    pub(crate) fn log2_nonewton(self) -> Self
    {
        if EXP_BASE == 2
        {
            return self.logb()
        }
        self.ln_nonewton()/Self::LN_2()
    }

    /// Returns the base 2 logarithm of the number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let two = FpDouble::from(2.0);
    ///
    /// // log2(2) - 1 == 0
    /// let abs_difference = (two.log2() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn log2(self) -> Self
    {
        let mut y = self.log2_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            let one = Self::one();
            let ln2 = Self::LN_2();
            for _ in 0..NEWTON
            {
                let xdyexp2 = self/y.exp2();
                let (diff, s) = one.sub_extra_sign(xdyexp2);
                let dy = diff/ln2;
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