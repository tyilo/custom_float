use core::{cmp::Ordering, num::FpCategory, ops::{Add, AddAssign}};

use num_traits::Zero;

use crate::{fp::as_lossless, util::{self, Saturation}, AnyInt, Fp, FpRepr, UInt};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub(super) fn abs_add_mantissas(exp: &mut U, mut mantissa1: U, mut mantissa2: U, neg: bool) -> U
    {
        if neg
        {
            if mantissa1 < mantissa2
            {
                mantissa2 - mantissa1
            }
            else
            {
                mantissa1 - mantissa2
            }
        }
        else 
        {
            match mantissa1.checked_add(&mantissa2)
            {
                Some(sum) => sum,
                None => match Self::mantissa_pairs_div_both_base_add(&mut mantissa1, &mut mantissa2)
                {
                    Ok(sum) => {
                        *exp = *exp + U::one();
                        sum
                    },
                    _ => U::max_value()
                }
            }
        }
    }

    pub(super) fn add_signs(sign1: bool, sign2: bool, mantissa1: U, mantissa2: U) -> bool
    {
        ((sign1 != sign2) && mantissa1 < mantissa2) != sign1
    }

    pub(super) fn max_exponents(mut exp1: U, mut exp2: U, mantissa1: &mut U, mantissa2: &mut U) -> U
    {
        if EXP_BASE == 0
        {
            if exp1 == exp2
            {
                return exp1
            }
            *mantissa1 = U::zero();
            *mantissa2 = U::zero();
            return Self::max_exponent_bits()
        }
        Self::normalize_mantissa_down(&mut exp1, mantissa1, Some(exp2));
        Self::normalize_mantissa_down(&mut exp2, mantissa2, Some(exp1));
        let (exp, shr, mantissa) = match exp1.cmp(&exp2)
        {
            Ordering::Less => (
                exp2,
                exp2 - exp1,
                mantissa1
            ),
            Ordering::Equal => return exp1,
            Ordering::Greater => (
                exp1,
                exp1 - exp2,
                mantissa2
            )
        };
        *mantissa = if let Some(base) = U::from(EXP_BASE)
        {
            util::rounding_div_pow(*mantissa, base, shr)
        }
        else
        {
            U::zero()
        };
        exp
    }
    
    pub(crate) fn add_with_sign_extra_sign(self, s_lhs: bool, rhs: Self, s_rhs: bool) -> (Self, bool)
    {
        let s0 = self.is_sign_negative() ^ s_lhs;
        let s1 = rhs.is_sign_negative() ^ s_rhs;

        match (self.classify(), rhs.classify())
        {
            (FpCategory::Nan, _) | (_, FpCategory::Nan) => (self.add_nan(rhs), false),
            (FpCategory::Infinite, FpCategory::Infinite) if s0 != s1 => (Self::qnan(), false),
            (_, FpCategory::Zero) | (FpCategory::Infinite, _) => (self, s_lhs),
            (FpCategory::Zero, _) | (_, FpCategory::Infinite) => (rhs, s_rhs),
            (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => {
                let e0 = self.exp_bits();
                let e1 = rhs.exp_bits();
                let mut f0 = self.mantissa_bits();
                let mut f1 = rhs.mantissa_bits();
                
                let mut e = Self::max_exponents(e0, e1, &mut f0, &mut f1);
                let mut f = Self::abs_add_mantissas(&mut e, f0, f1, s0 != s1);
                let s = Self::add_signs(s0, s1, f0, f1);
                Self::normalize_mantissa(&mut e, &mut f, None);
                (Self::from_exp_mantissa(e, f), s)
            }
        }
    }

    pub(crate) fn add_with_sign(self, s_lhs: bool, rhs: Self, s_rhs: bool) -> Self
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [match (s_lhs, s_rhs)
            {
                (false, false) => lhs + rhs,
                (true, false) => rhs - lhs,
                (false, true) => lhs - rhs,
                (true, true) => -(lhs + rhs)
            }],
            {
                let (y, s) = self.add_with_sign_extra_sign(s_lhs, rhs, s_rhs);
                y.xor_sign(s)
            }
        )
    }

    /*pub(crate) fn add_extra_sign(self, rhs: Self) -> (Self, bool)
    {
        self.add_with_sign_extra_sign(false, rhs, false)
    }*/

    pub fn add_int_diff<I: UInt>(self, rhs_add: I, rhs_sub: I) -> Self
    {
        if rhs_add < rhs_sub
        {
            self.sub_int(rhs_sub - rhs_add)
        }
        else
        {
            self.add_int(rhs_add - rhs_sub)
        }
    }

    pub fn add_int<I: AnyInt>(self, rhs: I) -> Self
    {
        self.add_int_with_sign(false, rhs, false)
    }

    pub(crate) fn add_int_with_sign<I: AnyInt>(self, s_lhs: bool, rhs: I, s_rhs: bool) -> Self
    {
        let (y, s) = self.add_int_with_sign_extra_sign(s_lhs, rhs, s_rhs);
        y.xor_sign(s)
    }

    pub(crate) fn add_int_with_sign_extra_sign<I: AnyInt>(self, s_lhs: bool, mut rhs: I, s_rhs: bool) -> (Self, bool)
    {
        if rhs.is_zero()
        {
            return (self, false)
        }
        let s0 = self.is_sign_negative() ^ s_lhs;
        let s1 = (rhs < I::zero()) ^ s_rhs;

        match self.classify()
        {
            FpCategory::Infinite | FpCategory::Nan => (self, s_lhs),
            FpCategory::Zero => (Self::from_int(rhs), s_rhs),
            FpCategory::Normal | FpCategory::Subnormal => {
                let e0 = self.exp_bits();
                let mut e1 = Self::exp_bias();
                let mut f0 = self.mantissa_bits();
                if !FRAC_SIZE.is_zero()
                {
                    match Self::mantissa_shr(&mut rhs, FRAC_SIZE, &mut e1)
                    {
                        Ok(()) => (),
                        Err(Saturation::Overflow) => return (Self::infinity(), s_rhs),
                        Err(Saturation::Underflow) => return (self, s_lhs)
                    }
                }
                let mut f1 = loop
                {
                    if let Some(x) = if util::is_signed::<I>() && rhs == I::min_value()
                    {
                        U::from(I::max_value()).and_then(|u| u.checked_add(&U::one()))
                    }
                    else
                    {
                        U::from(util::abs(rhs))
                    }
                    {
                        break x
                    }
                    if let Some(ee) = e1.checked_add(&U::one())
                    {
                        match I::from(EXP_BASE)
                        {
                            Some(base) => {
                                rhs = util::rounding_div(rhs, base);
                                e1 = ee;
                                continue
                            },
                            None => return (self, s_lhs)
                        }
                    }
                    return (Self::infinity(), s_rhs) 
                };
                
                let mut e = Self::max_exponents(e0, e1, &mut f0, &mut f1);
                let mut f = Self::abs_add_mantissas(&mut e, f0, f1, s0 != s1);
                let s = Self::add_signs(s0, s1, f0, f1);
                Self::normalize_mantissa(&mut e, &mut f, None);
                (Self::from_exp_mantissa(e, f), s)
            }
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Add for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output
    {
        self.add_with_sign(false, rhs, false)
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Add<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: U) -> Self::Output
    {
        self.add_int(rhs)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AddAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn add_assign(&mut self, rhs: Self)
    {
        *self = *self + rhs
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> AddAssign<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn add_assign(&mut self, rhs: U)
    {
        *self = *self + rhs
    }
}

#[cfg(test)]
mod test
{
    use core::ops::Add;
    use test::Bencher;

    use crate::{ieee754::FpDouble, tests::{bench_op2, test_op2}};

    #[test]
    fn test_add_once()
    {
        type F = FpDouble;

        let a = F::from(1f32);
        let b = F::from(-16f32);
        let c = a + b;
        println!("{a} + {b} = {c}");
    }
    #[test]
    fn test_add()
    {
        test_op2!("add", Add::add, Some(0.001))
    }
    #[bench]
    fn bench_add(bencher: &mut Bencher)
    {
        test_add();
        bench_op2!(bencher, Add::add)
    }
}