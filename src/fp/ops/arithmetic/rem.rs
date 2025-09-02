use core::{cmp::Ordering, num::FpCategory, ops::{Rem, RemAssign}};

use crate::{fp::as_lossless, util::Saturation, AnyInt, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn rem_normal<I: AnyInt>(self, exp2: U, mut mantissa2: I, mut shift: usize) -> Self
    {
        let s = self.is_sign_negative();
        let mut e1 = self.exp_bits();
        let mut e2 = exp2;

        macro_rules! sat {
            ($expr:expr) => {
                match $expr
                {
                    Ok(y) => y,
                    Err(Saturation::Overflow) => return self,
                    Err(Saturation::Underflow) => return Self::nan().with_sign(s)
                }
            };
        }

        if INT_SIZE <= 1 && e1 < e2
        {
            return self
        }
        if EXP_BASE != 0 && EXP_BASE.is_power_of_two()
            && let neg = e1 < e2
            && let de = if neg {e1 - e2} else {e2 - e1}
            && let Some(s) = de.to_usize().and_then(|ss| ss.checked_mul(EXP_BASE.ilog2() as usize))
        {
            if neg
            {
                if let Some(ss) = shift.checked_sub(s)
                {
                    e2 = e1;
                    shift = ss
                }
                else if let Some(d) = U::from((s - shift)/EXP_BASE.ilog2() as usize)
                {
                    e2 = e1 + d;
                    shift = 0
                }
            }
            else
            {
                if let Some(ss) = shift.checked_add(s)
                {
                    e2 = e1;
                    shift = ss
                }
                else
                {
                    return self
                }
            }
        };
        sat!(Self::mantissa_shl(&mut mantissa2, shift, &mut e2, None));
        let mut f2 = sat!(Self::convert_mantissa(mantissa2, &mut e2, None));
        Self::normalize_mantissa_up(&mut e2, &mut f2, Some(e1));

        if e1 == e2
        {
            let mut f1 = self.frac_bits();
            match f1.cmp(&f2)
            {
                Ordering::Less => return self,
                Ordering::Equal => return Self::zero().with_sign(s),
                Ordering::Greater => {
                    f1 = f1 % f2;
                    Self::normalize_mantissa(&mut e1, &mut f1, None);
                    return Self::from_exp_mantissa(e1, f1)
                }
            }
        }

        Self::normalize_mantissa(&mut e2, &mut f2, None);

        let fa = self.abs();
        let fb = Self::from_exp_mantissa(e2, f2);

        if fa < fb
        {
            return self
        }
    
        let mut dividend = fa;
        /* normalize divisor */
        let e1 = fa.exp_bits();


        if !Self::IS_INT_IMPLICIT
        {
            f2 = f2 + (fb.int_bits() << Self::INT_POS)
        }
        let mut divisor = Self::from_bits(Self::shift_exp(e1) + Self::shift_frac(f2));
        let min_divisor = dividend.divb();
        while divisor < min_divisor
        {
            divisor = divisor.mulb();
        }
        /* compute quotient one bit at a time */
        while divisor >= fb && matches!(divisor.classify(), FpCategory::Normal | FpCategory::Subnormal)
        {
            while dividend >= divisor
            {
                dividend -= divisor;
            }
            divisor = divisor.divb();
        }
        assert!(dividend.is_finite());
        /* dividend now represents remainder */
        dividend.with_sign(s)
    }

    pub fn rem_int<I: AnyInt>(self, rhs: I) -> Self
    {
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Zero | FpCategory::Subnormal => self,
            FpCategory::Infinite => Self::nan().copysign(self),
            FpCategory::Normal => self.rem_normal(Self::exp_bias(), rhs, FRAC_SIZE)
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Rem for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output
    {
        as_lossless!(
            [self, rhs],
            |[lhs, rhs]| [lhs % rhs],
            {
                let s = self.is_sign_negative();
                match (self.classify(), rhs.classify())
                {
                    (FpCategory::Nan, _) | (_, FpCategory::Nan) => Self::add_nan(self, rhs.with_sign(s)),
                    (FpCategory::Infinite, _) | (_, FpCategory::Zero) => Self::nan(),
                    (FpCategory::Zero, _) | (_, FpCategory::Infinite) => self,
                    (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => self.rem_normal(rhs.exp_bits(), rhs.mantissa_bits(), 0)
                }
            }
        )
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Rem<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    fn rem(self, rhs: U) -> Self::Output
    {
        self.rem_int(rhs)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> RemAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn rem_assign(&mut self, rhs: Self)
    {
        *self = *self % rhs
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> RemAssign<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn rem_assign(&mut self, rhs: U)
    {
        *self = *self % rhs
    }
}

#[cfg(test)]
mod test
{
    use core::ops::Rem;
    use test::Bencher;

    use crate::tests::{bench_op2, test_op2};

    #[test]
    fn test_rem()
    {
        test_op2!("rem", Rem::rem, None)
    }
    #[bench]
    fn bench_rem(bencher: &mut Bencher)
    {
        test_rem();
        bench_op2!(bencher, Rem::rem)
    }
}