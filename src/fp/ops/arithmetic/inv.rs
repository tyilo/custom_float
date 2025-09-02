use core::num::FpCategory;

use num_traits::Inv;

use crate::{Fp, FpRepr, fp::as_lossless};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Take the reciprocal (inverse) of a number, `1/x`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(2.0);
    /// let abs_difference = (x.recip() - (FpDouble::one()/x)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[doc(alias = "inv")]
    #[inline]
    pub fn recip(self) -> Self
    {
        as_lossless!(
            [self],
            |[x]| [x.recip()],
            |[x]| [1.0/x],
            {
                let s = self.is_sign_negative();
                match self.classify()
                {
                    FpCategory::Nan => self,
                    FpCategory::Zero => Self::infinity().with_sign(s),
                    FpCategory::Infinite => Self::zero().with_sign(s),
                    FpCategory::Normal | FpCategory::Subnormal => {
                        if self.abs().is_one()
                        {
                            return self
                        }
                
                        let bias = Self::exp_bias();
                        let mut e: U = bias + bias - self.exp_bits();
                        let f0: U = U::one() << Self::INT_POS;
                        let f1: U = self.mantissa_bits();
                        
                        let mut o = U::zero();
                        let mut f = match Self::mantissa_div(f0, f1, &mut e, &mut o)
                        {
                            Ok(f) => f,
                            Err(done) => return done.with_sign(s)
                        };
                        Self::normalize_mantissa_up(&mut e, &mut f, Some(o));
                        let mut e = match e.checked_sub(&o)
                        {
                            Some(e) => e,
                            None => return Self::zero().with_sign(s)
                        };
                
                        Self::normalize_mantissa(&mut e, &mut f, None);
                        Self::from_exp_mantissa(e, f).with_sign(s)
                    }
                }
            }
        )
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Inv for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn inv(self) -> Self::Output
    {
        self.recip()
    }
}

#[cfg(test)]
mod test
{
    use num_traits::Inv;
    use test::Bencher;

    use crate::{khronos::KhronosFp10, tests::{bench_op1, test_op1}};

    #[test]
    fn test_inv_once()
    {
        type F = KhronosFp10;

        let x = -F::one();
        let y = x.inv();

        println!("1 / {x} = {y}")
    }

    #[test]
    fn test_inv()
    {
        test_op1!("inv", Inv::inv, None, Some(0.1..10.0))
    }
    #[bench]
    fn bench_inv(bencher: &mut Bencher)
    {
        test_inv();
        bench_op1!(bencher, Inv::inv)
    }
}