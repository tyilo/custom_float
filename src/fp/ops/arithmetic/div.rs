use core::{cmp::Ordering, num::FpCategory, ops::{Div, DivAssign}};

use num_traits::Zero;

use crate::{fp::as_lossless, util::{self, Saturation}, AnyInt, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub(super) fn exponent_sub(exp1: U, exp2: U, mantissa1: &mut U, mantissa2: &mut U) -> Result<U, Self>
    {
        let bias = Self::exp_bias();
        if let Some(e) = exp1.checked_sub(&exp2)
        {
            return e.checked_add(&bias).ok_or_else(Self::infinity)
        }
        if let Some(e) = bias.checked_sub(&exp2)
        {
            return exp1.checked_add(&e).ok_or_else(Self::infinity)
        }
        let e = exp2 - bias;
        if let Some(e) = exp1.checked_sub(&e)
        {
            return Ok(e)
        }
        let mut o = e - exp1;
        if !o.is_zero() && EXP_BASE.is_zero()
        {
            return Err(match exp1.cmp(&exp2)
            {
                Ordering::Greater => Self::infinity(),
                Ordering::Equal => Self::nan(),
                Ordering::Less => Self::zero()
            })
        }
        'lp1:
        while o > U::zero()
        {
            o = o - U::one();
            let lz = mantissa2.leading_zeros() as usize;
            if let Some(base) = U::from(EXP_BASE)
            {
                if (*mantissa1 % base).is_zero()
                {
                    *mantissa1 = *mantissa1/base;
                    continue 'lp1
                }
                else if lz > Self::BASE_PADDING
                {
                    *mantissa2 = *mantissa2*base;
                    continue 'lp1
                }
            }
            if !Self::BASE_FACTORS_2.is_empty()
                && let Some((_, [b1, b2])) = Self::BASE_FACTORS_2.iter()
                    .zip(Self::BASE_FACTORS_2_PADDING)
                    .filter_map(|(&[base1, base2], &[p1, _])| if let Some(b1) = U::from(base1)
                        && let Some(b2) = U::from(base2)
                        && lz > p1
                    {
                        Some(((*mantissa1 % b1).is_zero(), [b1, b2]))
                    }
                    else
                    {
                        None
                    }).reduce(|a, b| if b.0 && !a.0
                    {
                        b
                    }
                    else
                    {
                        a
                    })
            {
                *mantissa1 = util::rounding_div(*mantissa1, b2);
                *mantissa2 = *mantissa2*b1;
                continue 'lp1
            }
            if let Some(base) = U::from(EXP_BASE)
            {
                *mantissa1 = util::rounding_div(*mantissa1, base);
                continue 'lp1
            }
            return Err(Self::zero())
        }
        Ok(U::zero())
    }

    pub(super) fn integral_div<M: AnyInt>(mut mantissa1: U, mut mantissa2: M, exp: &mut U, exp_offset: &mut U) -> Result<U, Self>
    {
        let mut div = match Self::convert_mantissa(mantissa2, exp_offset, Some(exp))
        {
            Ok(m) => m,
            Err(Saturation::Overflow) => return Err(Self::zero()),
            Err(Saturation::Underflow) => return Err(Self::infinity())
        };
        if EXP_BASE.is_power_of_two() && !mantissa2.is_zero()
            && (
                util::complementary_add_sub_assign(Some(exp_offset), Some(exp), Self::shl_mantissa_without_loss(&mut mantissa1, None, EXP_BASE.ilog2(), None)).is_err()
                || util::complementary_add_sub_assign(Some(exp_offset), Some(exp), Self::shr_mantissa_without_loss(&mut mantissa2, None, EXP_BASE.ilog2(), None)).is_err()
            )
        {
            return Err(Self::zero())
        }
        loop
        {
            if div.is_zero()
            {
                return Err(Self::infinity())
            }
            let f = util::rounding_div(mantissa1, div);
            if let Some(base) = U::from(EXP_BASE) && f.leading_zeros() as usize >= Self::BASE_PADDING
            {
                if mantissa1.leading_zeros() as usize > Self::BASE_PADDING
                {
                    mantissa1 = mantissa1*base
                }
                else
                {
                    if div % base != U::zero()
                    {
                        break Ok(f)
                    }
                    div = div/base;
                }
                util::complementary_add_sub_assign(Some(exp_offset), Some(exp), U::one())
                    .map_err(|()| Self::zero())?
            }
            else
            {
                break Ok(f)
            }
        }
    }

    pub(super) fn mantissa_div(mut mantissa1: U, mut mantissa2: U, exp: &mut U, exp_offset: &mut U) -> Result<U, Self>
    {
        let mut shifts = FRAC_SIZE;
        shifts -= Self::shl_mantissa_without_loss::<_, usize>(&mut mantissa1, Some(shifts), 1, None);
        shifts -= Self::shr_mantissa_without_loss::<_, usize>(&mut mantissa2, Some(shifts), 1, None);
        let mut mantissa = Self::integral_div(mantissa1, mantissa2, exp, exp_offset)?;
        Self::mantissa_shl(&mut mantissa, shifts, exp, Some(exp_offset))
            .map_err(|sat| match sat
            {
                Saturation::Overflow => Self::infinity(),
                Saturation::Underflow => Self::zero()
            })?;
        Ok(mantissa)
    }

    pub fn div_int<I: AnyInt>(self, rhs: I) -> Self
    {
        let s = self.is_sign_negative() ^ (rhs < I::zero());
        match self.classify()
        {
            FpCategory::Zero if rhs.is_zero() => Self::nan().with_sign(s),
            FpCategory::Nan | FpCategory::Infinite | FpCategory::Zero => self.with_sign(s),
            FpCategory::Subnormal | FpCategory::Normal => {
                let mut e = self.exp_bits();
                let mut f = self.mantissa_bits();
                let mut o = U::zero();
        
                if rhs.to_usize() == Some(EXP_BASE) && let Some(ee) = e.checked_sub(&U::one())
                {
                    e = ee
                }
                else
                {
                    f = match Self::integral_div(f, rhs, &mut e, &mut o)
                    {
                        Ok(ff) => ff,
                        Err(done) => return done.with_sign(s)
                    };
                }
                
                Self::normalize_mantissa_up(&mut e, &mut f, Some(o));
                let mut e = match e.checked_sub(&o)
                {
                    Some(e) => e,
                    None => return Self::zero().with_sign(s)
                };
                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_exp_mantissa(e, f).with_sign(s)
            },
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Div for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [lhs/rhs],
            {
                let s = self.is_sign_negative()^rhs.is_sign_negative();
                match (self.classify(), rhs.classify())
                {
                    (FpCategory::Nan, _) | (_, FpCategory::Nan) => self.add_nan(rhs).with_sign(s),
                    (FpCategory::Zero, FpCategory::Zero) | (FpCategory::Infinite, FpCategory::Infinite) => Self::qnan().with_sign(s),
                    (FpCategory::Infinite, _) | (_, FpCategory::Zero) => Self::infinity().with_sign(s),
                    (_, FpCategory::Infinite) | (FpCategory::Zero, _) => Self::zero().with_sign(s),
                    (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => {
                        if rhs.abs().is_one()
                        {
                            return self.with_sign(s)
                        }
                        if self.abs().is_one()
                        {
                            return rhs.recip().with_sign(s)
                        }
                
                        let e0: U = self.exp_bits();
                        let e1: U = rhs.exp_bits();
                
                        let mut f0: U = self.mantissa_bits();
                        let mut f1: U = rhs.mantissa_bits();
                
                        if e0 == e1 && f0 == f1
                        {
                            return Self::one().with_sign(s)
                        }
                
                        let mut e = match Self::exponent_sub(e0, e1, &mut f0, &mut f1)
                        {
                            Ok(e) => e,
                            Err(done) => return done.with_sign(s)
                        };
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
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Div<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    fn div(self, rhs: U) -> Self::Output
    {
        self.div_int(rhs)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> DivAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn div_assign(&mut self, rhs: Self)
    {
        *self = *self / rhs
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> DivAssign<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn div_assign(&mut self, rhs: U)
    {
        *self = *self / rhs
    }
}

#[cfg(test)]
mod test
{
    use std::ops::Div;

    use test::Bencher;

    use crate::{ieee754::FpDouble, tests::{bench_op2, test_op2}};

    #[test]
    fn test_div_once()
    {
        type F = FpDouble;
        
        let a = F::from(4.0);
        let b = F::from(2.0);
        let c = a / b;
        println!("{a} / {b} = {c}");
    }

    #[test]
    fn test_div()
    {
        test_op2!("div", Div::div, None)
    }
    #[bench]
    fn bench_div(bencher: &mut Bencher)
    {
        test_div();
        bench_op2!(bencher, Div::div)
    }
}