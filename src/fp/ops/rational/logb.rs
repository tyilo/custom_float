use core::num::FpCategory;

use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the logarithm base `EXP_BASE` of the number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpDouble, DecDouble};
    ///
    /// let two = FpDouble::from(2.0);
    /// let ten = DecDouble::from(10.0);
    ///
    /// // log2(2) - 1 == 0
    /// let abs_difference_2 = (two.logb() - FpDouble::one()).abs();
    /// 
    /// // log10(10) - 1 == 0
    /// let abs_difference_10 = (ten.logb() - DecDouble::one()).abs();
    ///
    /// assert!(abs_difference_2 < FpDouble::from(1e-10));
    /// assert!(abs_difference_10 < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn logb(self) -> Self
    {
        match self.classify()
        {
            FpCategory::Nan => self,
            FpCategory::Zero => Self::neg_infinity(),
            FpCategory::Infinite => if self.is_sign_negative()
            {
                Self::nan()
            }
            else
            {
                self
            },
            FpCategory::Normal | FpCategory::Subnormal => {
                let mut e = self.exp_bits();
                let mut f = self.mantissa_bits();

                let bias = Self::exp_bias();
                let mut o = bias;
                
                let base = U::from(EXP_BASE).unwrap();
                if e.is_zero() && Self::IS_INT_IMPLICIT
                {
                    while (f >> (Self::MANTISSA_OP_SIZE - Self::BASE_PADDING)).is_zero()
                        && let Some(ff) = f.checked_mul(&base)
                        && util::complementary_add_sub_assign(Some(&mut o), Some(&mut e), U::one()).is_ok()
                    {
                        f = ff;
                    }
                }
                if f != Self::shift_int(U::one())
                {
                    // TODO: Avoid conversion to f64
                    let u = Self::from_exp_mantissa(bias, f);
                    let u: f64 = u.into();

                    Self::from(u.log(EXP_BASE as f64)).add_int_diff(e, o)
                }
                else
                {
                    Self::from_int_diff(e, o)
                }
            }
        }
    }
}