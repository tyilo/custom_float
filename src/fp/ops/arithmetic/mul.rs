use core::{num::FpCategory, ops::{Mul, MulAssign}};

use crate::{util::{self, Saturation}, AnyInt, Fp, FpRepr, fp::as_lossless};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub(super) fn exponent_add(mut exp1: U, mut exp2: U, mantissa1: &mut U, mut mantissa2: Option<&mut U>) -> Result<U, Self>
    {
        let bias = Self::exp_bias();
        if let Some(e) = exp1.checked_add(&exp2)
        {
            if let Some(e) = e.checked_sub(&bias)
            {
                return Ok(e)
            }
            let mut o = bias - e;
            if let Some(m2) = &mut mantissa2
            {
                while o > U::zero()
                {
                    match Self::mantissa_pairs_div_either_base(mantissa1, *m2)
                    {
                        Ok(()) => o = o - U::one(),
                        Err(()) => return Err(Self::zero())
                    }
                }
            }
            else
            {
                while o > U::zero()
                {
                    match Self::mantissa_div_sqrt_base_or_base(mantissa1)
                    {
                        Ok(add_exp) => if let Some(oo) = o.checked_sub(&add_exp)
                        {
                            o = oo
                        }
                        else
                        {
                            return Ok(add_exp - o)
                        },
                        Err(()) => return Err(Self::zero())
                    }
                }
            }
            return Ok(U::zero())
        }
        let exp_most = if exp1 < exp2
        {
            &mut exp2
        }
        else
        {
            &mut exp1
        };
        *exp_most = *exp_most - bias;
        if let Some(e) = exp1.checked_add(&exp2)
        {
            return Ok(e)
        }
        Err(Self::infinity())
    }

    pub(super) fn integral_mul<I: AnyInt>(mut mantissa1: U, mantissa2: I, exp: &mut U) -> Result<U, Self>
    {
        let mut mul = match Self::convert_mantissa(mantissa2, exp, None)
        {
            Ok(m) => m,
            Err(Saturation::Overflow) => return Err(Self::infinity()),
            Err(Saturation::Underflow) => return Err(Self::zero())
        };
        if EXP_BASE.is_power_of_two()
        {
            let (mantissa, overflow) = util::widening_mul(mantissa1, mul);
            return Ok(Self::mantissa_from_low_high(mantissa, overflow, exp))
        }
        // TODO: This loop is slow!
        if let Some(f) = mantissa1.checked_mul(&mul)
        {
            return Ok(f)
        }
        loop
        {
            match Self::mantissa_pairs_div_either_base_mul(&mut mantissa1, &mut mul)
            {
                Ok(result) => {
                    *exp = *exp + U::one();
                    if let Some(y) = result
                    {
                        return Ok(y)
                    }
                },
                Err(()) => return Ok(Self::max_mantissa_bits())
            }
        }
    }

    pub(super) fn mantissa_mul(mut mantissa1: U, mut mantissa2: U, exp: &mut U) -> Result<U, Self>
    {
        let mut shifts = FRAC_SIZE;
        shifts -= Self::shr_mantissa_without_loss::<_, usize>(&mut mantissa1, Some(shifts), 1, None);
        shifts -= Self::shr_mantissa_without_loss::<_, usize>(&mut mantissa2, Some(shifts), 1, None);
        let mut mantissa = Self::integral_mul(mantissa1, mantissa2, exp)?;
        match Self::mantissa_shr(&mut mantissa, shifts, exp)
        {
            Ok(()) => (),
            Err(sat) => match sat
            {
                Saturation::Overflow => return Err(Self::infinity()),
                Saturation::Underflow => return Err(Self::zero())
            }
        }
        Ok(mantissa)
    }

    pub fn mul_int<I: AnyInt>(self, rhs: I) -> Self
    {
        let s = self.is_sign_negative() ^ (rhs < I::zero());
        if rhs.is_zero()
        {
            return Self::zero().with_sign(s)
        }
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Infinite | FpCategory::Zero => self,
            FpCategory::Subnormal | FpCategory::Normal => {
                let mut e = self.exp_bits();
                let mut f = self.mantissa_bits();

                if let Some(base) = I::from(EXP_BASE) && base == rhs && let Some(ee) = e.checked_add(&U::one())
                {
                    e = ee
                }
                else
                {
                    f = match Self::integral_mul(f, rhs, &mut e)
                    {
                        Ok(ff) => ff,
                        Err(done) => return done
                    }
                }
        
                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_exp_mantissa(e, f).with_sign(s)
            }
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Mul<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    fn mul(self, rhs: U) -> Self::Output
    {
        let s = self.is_sign_negative();
        if rhs.is_zero()
        {
            return Self::zero().with_sign(s)
        }
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Infinite | FpCategory::Zero => self,
            FpCategory::Subnormal | FpCategory::Normal => {
                let mut e = self.exp_bits();
                let mut f = self.mantissa_bits();

                if let Some(base) = U::from(EXP_BASE) && base == rhs && let Some(ee) = e.checked_add(&U::one())
                {
                    e = ee
                }
                else
                {
                    f = match Self::integral_mul(f, rhs, &mut e)
                    {
                        Ok(ff) => ff,
                        Err(done) => return done
                    }
                }
        
                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_exp_mantissa(e, f).with_sign(s)
            }
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Mul<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [lhs*rhs],
            {
                // Slow! Better to avoid branching.
                /*if self.to_bits() == rhs.to_bits()
                {
                    return self.squared()
                }*/
                let s = self.is_sign_negative()^rhs.is_sign_negative();
                match (self.classify(), rhs.classify())
                {
                    (FpCategory::Nan, _) | (_, FpCategory::Nan) => self.add_nan(rhs).with_sign(s),
                    (FpCategory::Infinite, FpCategory::Zero) | (FpCategory::Zero, FpCategory::Infinite) => Self::qnan().with_sign(s),
                    (FpCategory::Zero, _) | (_, FpCategory::Zero) => Self::zero().with_sign(s),
                    (FpCategory::Infinite, _) | (_, FpCategory::Infinite) => Self::infinity().with_sign(s),
                    (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => {
                        if rhs.abs().is_one()
                        {
                            return self.with_sign(s)
                        }
                        if self.abs().is_one()
                        {
                            return rhs.with_sign(s)
                        }

                        let e0: U = self.exp_bits();
                        let e1: U = rhs.exp_bits();
                
                        let mut f0: U = self.mantissa_bits();
                        let mut f1: U = rhs.mantissa_bits();
                
                        let mut e = match Self::exponent_add(e0, e1, &mut f0, Some(&mut f1))
                        {
                            Ok(e) => e,
                            Err(done) => return done.with_sign(s)
                        };
                        let mut f = match Self::mantissa_mul(f0, f1, &mut e)
                        {
                            Ok(e) => e,
                            Err(done) => return done.with_sign(s)
                        };
                
                        Self::normalize_mantissa(&mut e, &mut f, None);
                        Self::from_exp_mantissa(e, f).with_sign(s)
                    }
                }
            }
        )
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> MulAssign<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn mul_assign(&mut self, rhs: U)
    {
        *self = *self * rhs
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> MulAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn mul_assign(&mut self, rhs: Self)
    {
        *self = *self * rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Mul;

    use test::Bencher;

    use crate::{ieee754::FpHalf, tests::{bench_op2, test_op2}};
    
    #[test]
    fn test_mul_once()
    {
        type F = FpHalf;

        let a = F::one();
        let b = F::from(8388608);
        let c = a * b;
        println!("{a} * {b} = {c}");
    }
    #[test]
    fn test_mul()
    {
        test_op2!("mul", Mul::mul, Some(0.0001))
    }
    #[bench]
    fn bench_mul(bencher: &mut Bencher)
    {
        test_mul();
        bench_op2!(bencher, Mul::mul)
    }
}