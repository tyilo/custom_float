use core::num::FpCategory;

use num_traits::Zero;

use crate::{Fp, FpRepr};


impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns `self*EXP_BASE`.
    ///
    /// This is generally faster than using regular multiplication.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpDouble, DecDouble};
    ///
    /// let f = FpDouble::from(2.0);
    /// let d = DecDouble::from(2.0);
    ///
    /// // 2*2 - 4 == 0
    /// let abs_difference_f = (f.mulb() - FpDouble::from(4.0)).abs();
    ///
    /// // 2*10 - 20 == 0
    /// let abs_difference_d = (d.mulb() - DecDouble::from(20.0)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn mulb(self) -> Self
    {
        let s = self.is_sign_negative();
        if EXP_BASE.is_zero()
        {
            return Self::zero().with_sign(s)
        }
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Infinite | FpCategory::Zero => self,
            FpCategory::Subnormal | FpCategory::Normal => {
                let mut e = self.exp_bits();
                let mut f = self.mantissa_bits();
        
                e = e + U::one();
        
                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_exp_mantissa(e, f).with_sign(s)
            }
        }
    }
}