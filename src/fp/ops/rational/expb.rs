use core::num::FpCategory;

use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns `EXP_BASE^(self)`.
    ///
    /// This implementation is roughly based on the exp2 implementation described here: https://stackoverflow.com/questions/65554112/fast-double-exp2-function-in-c.
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
    /// // 2^2 - 4 == 0
    /// let abs_difference_f = (f.expb() - FpDouble::from(4.0)).abs();
    ///
    /// // 10^2 - 100 == 0
    /// let abs_difference_d = (d.expb() - DecDouble::from(100.0)).abs();
    ///
    /// assert!(abs_difference_f < FpDouble::from(1e-10));
    /// assert!(abs_difference_d < DecDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn expb(self) -> Self
    {
        match self.classify()
        {
            FpCategory::Nan => self,
            FpCategory::Infinite => {
                if self.is_sign_negative()
                {
                    return Self::zero()
                }
                self
            },
            FpCategory::Zero => Self::one(),
            FpCategory::Subnormal | FpCategory::Normal => {
                match self.is_sign_negative()
                {
                    true => if let Some(exp_frac) = U::from(util::exp2_ilog(FRAC_SIZE + INT_SIZE, EXP_BASE))
                        && self <= -Self::from_int((U::one() << (EXP_SIZE - 1)) + exp_frac - U::one())
                    {
                        return -Self::from_bits(U::one() << (Self::EXP_POS - EXP_SIZE/2))*self;
                    },
                    false => if self >= Self::from_int(U::one() << (EXP_SIZE - 1))
                    {
                        return Self::from_bits((Self::max_exponent_bits()) << Self::EXP_POS)*self;
                    }
                }
        
                let inv = self.is_sign_negative();
                let bias = Self::from_int(Self::exp_bias());
                let mut x = self.abs();
                let max = Self::exp_bias() >> 1usize;
                let mf = (x/max).ceil();
                let m = if !mf.is_zero()
                {
                    if let Some(m) = mf.to_int::<u32>()
                    {
                        x /= mf;
                        Some(m)
                    }
                    else
                    {
                        return if inv {Self::zero()} else {Self::infinity()}
                    }
                }
                else
                {
                    None
                };
                
                x += bias;
                let e = x.trunc();
                let f = x - e;
        
                let z = if f.is_zero()
                {
                    None
                }
                else
                {
                    // This would otherwise be solved by a LUT
                    Some(Self::from(match EXP_BASE
                    {
                        2 => <f64 as From<_>>::from(f).exp2(),
                        _ => (EXP_BASE as f64).powf(f.into())
                    }))
                };
        
                if e.is_nan()
                {
                    return e
                }
                match if e.is_sign_negative()
                {
                    Err(!inv)
                }
                else if e > Self::from_int(Self::max_exponent_bits())
                {
                    Err(inv)
                }
                else if let Some(ee) = e.to_int::<U>()
                {
                    Ok(ee)
                }
                else
                {
                    Err(inv)
                }
                {
                    Err(sgn) => if sgn {Self::zero()} else {Self::infinity()},
                    Ok(ee) => {
                        let mut bits = Self::shift_exp(ee);
                        if !Self::IS_INT_IMPLICIT
                        {
                            bits = bits + Self::shift_int(U::one())
                        };
                        let mut y = Self::from_bits(bits);
                        if let Some(z) = z
                        {
                            y *= z
                        }
                        if inv
                        {
                            y = y.recip()
                        }
                        if let Some(m) = m
                        {
                            y = y.powu(m)
                        }
                        y
                    }
                }
            }
        }
    }
}