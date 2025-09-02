use core::num::FpCategory;

use num_traits::Zero;

use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns `self/EXP_BASE`.
    ///
    /// This is generally faster than using regular division.
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
    /// // 2/2 - 1 == 0
    /// let abs_difference_f = (f.divb() - FpDouble::from(1.0)).abs();
    ///
    /// // 2/10 - 0.2 == 0
    /// let abs_difference_d = (d.divb() - DecDouble::from(0.2)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn divb(self) -> Self
    {
        let s = self.is_sign_negative();
        match self.classify()
        {
            FpCategory::Zero if EXP_BASE.is_zero() => Self::nan().with_sign(s),
            FpCategory::Nan | FpCategory::Infinite | FpCategory::Zero => self,
            FpCategory::Subnormal | FpCategory::Normal => {
                if EXP_BASE.is_zero()
                {
                    return Self::infinity().with_sign(s)
                }
                let mut e = self.exp_bits();
                let mut f = self.mantissa_bits();
        
                if let Some(ee) = e.checked_sub(&U::one())
                {
                    e = ee
                }
                else if let Some(base) = U::from(EXP_BASE)
                {
                    f = util::rounding_div(f, base);
                }
                else
                {
                    return Self::zero().with_sign(s)
                }
                
                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_exp_mantissa(e, f).with_sign(s)
            },
        }
    }
}