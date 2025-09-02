use core::num::FpCategory;

use crate::{fp::NEWTON_RT, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Take the cubic root of a number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(8.0);
    ///
    /// // x^(1/3) - 2 == 0
    /// let abs_difference = (x.cbrt() - FpDouble::from(2.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn cbrt(self) -> Self
    {
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Infinite | FpCategory::Zero => self,
            FpCategory::Normal | FpCategory::Subnormal => {
                let third = Self::from(3u8).recip();
                let y = self.abs().powf_generic::<false>(third).copysign(self);
        
                const NEWTON: usize = NEWTON_RT;
                let two = Self::from(2u8);
                let mut y = y;
                for _ in 0..NEWTON
                {
                    let yy = third*(self/(y*y) + two*y);
                    if !yy.is_finite()
                    {
                        break
                    }
                    y = yy;
                }
                y
            }
        }
    }
}