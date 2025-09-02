use core::f128::consts::*;

use num_traits::FloatConst;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[allow(non_snake_case)]
    #[inline]
    #[doc = "Return `π / 180.0`."]
    pub(crate) fn FRAC_PI_180() -> Self
    {
        Self::FRAC_PI_6()/Self::from(30.0)
    }

    pub(crate) fn ln_exp_base() -> Self
    {
        match EXP_BASE
        {
            2 => Self::LN_2(),
            10 => Self::LN_10(),
            _ if EXP_BASE.is_power_of_two() => Self::from(EXP_BASE.ilog2())*Self::LN_2(),
            _ if util::is_power_of(EXP_BASE, 10) => Self::from(EXP_BASE.ilog10())*Self::LN_10(),
            _ => Self::from((EXP_BASE as f64).ln())
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> FloatConst for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    #[doc = "Return Euler’s number."]
    fn E() -> Self
    {
        Self::from(E)
    }

    #[inline]
    #[doc = "Return `1.0 / π`."]
    fn FRAC_1_PI() -> Self
    {
        Self::from(FRAC_1_PI)
    }

    #[inline]
    #[doc = "Return `1.0 / sqrt(2.0)`."]
    fn FRAC_1_SQRT_2() -> Self
    {
        Self::from(FRAC_1_SQRT_2)
    }

    #[inline]
    #[doc = "Return `2.0 / π`."]
    fn FRAC_2_PI() -> Self
    {
        Self::from(FRAC_2_PI)
    }

    #[inline]
    #[doc = "Return `2.0 / sqrt(π)`."]
    fn FRAC_2_SQRT_PI() -> Self
    {
        Self::from(FRAC_2_SQRT_PI)
    }

    #[inline]
    #[doc = "Return `π / 2.0`."]
    fn FRAC_PI_2() -> Self
    {
        Self::from(FRAC_PI_2)
    }

    #[inline]
    #[doc = "Return `π / 3.0`."]
    fn FRAC_PI_3() -> Self
    {
        Self::from(FRAC_PI_3)
    }

    #[inline]
    #[doc = "Return `π / 4.0`."]
    fn FRAC_PI_4() -> Self
    {
        Self::from(FRAC_PI_4)
    }

    #[inline]
    #[doc = "Return `π / 6.0`."]
    fn FRAC_PI_6() -> Self
    {
        Self::from(FRAC_PI_6)
    }

    #[inline]
    #[doc = "Return `π / 8.0`."]
    fn FRAC_PI_8() -> Self
    {
        Self::from(FRAC_PI_8)
    }

    #[inline]
    #[doc = "Return `ln(10.0)`."]
    fn LN_10() -> Self
    {
        Self::from(LN_10)
    }

    #[inline]
    #[doc = "Return `ln(2.0)`."]
    fn LN_2() -> Self
    {
        Self::from(LN_2)
    }

    #[inline]
    #[doc = "Return `log10(e)`."]
    fn LOG10_E() -> Self
    {
        Self::from(LOG10_E)
    }

    #[inline]
    #[doc = "Return `log2(e)`."]
    fn LOG2_E() -> Self
    {
        Self::from(LOG2_E)
    }

    #[inline]
    #[doc = "Return Archimedes’ constant `π`."]
    fn PI() -> Self
    {
        Self::from(PI)
    }

    #[inline]
    #[doc = "Return `sqrt(2.0)`."]
    fn SQRT_2() -> Self
    {
        Self::from(SQRT_2)
    }
}