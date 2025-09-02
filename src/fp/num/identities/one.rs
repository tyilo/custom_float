use core::{cmp::Ordering, num::FpCategory};

use num_traits::One;

use crate::{fp::Fps, Fp, FpRepr};


impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub(crate) fn total_cmp_one(self) -> Ordering
    {
        if self.is_sign_negative()
        {
            return Ordering::Less
        }
        self.abs_total_cmp_one()
    }

    pub(crate) fn abs_total_cmp_one(self) -> Ordering
    {
        let Fps(x, _) = self;
        x.abs_total_cmp_one()
    }

    pub(crate) fn cmp_one(self) -> Option<Ordering>
    {
        match self.classify()
        {
            FpCategory::Nan => None,
            FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero => Some(self.total_cmp_one())
        }
    }

    pub(crate) fn abs_cmp_one(self) -> Option<Ordering>
    {
        self.0.abs_cmp_one()
    }

    pub(crate) fn abs_is_one(self) -> bool
    {
        self.0.abs_is_one()
    }
    
    /// Returns the multiplicative identity element of `Self`, `1`.
    #[must_use]
    pub fn one() -> Self
    {
        Fp::one().into()
    }
    
    /// Sets `self` to the multiplicative identity element of `Self`, `1`.
    fn set_one(&mut self)
    {
        self.0.set_one();
        self.1 = false
    }

    /// Returns `true` if `self` is equal to the multiplicative identity.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_one(self) -> bool
    {
        self.is_sign_positive()
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub(crate) fn total_cmp_one(self) -> Ordering
    {
        if self.is_sign_negative()
        {
            return Ordering::Less
        }
        self.abs_total_cmp_one()
    }

    pub(crate) fn abs_total_cmp_one(self) -> Ordering
    {
        match self.classify()
        {
            FpCategory::Infinite | FpCategory::Nan => Ordering::Greater,
            FpCategory::Subnormal | FpCategory::Zero => Ordering::Less,
            FpCategory::Normal => {
                let e = self.exp_bits();

                if Self::IS_INT_IMPLICIT
                {
                    return self.abs_total_cmp(Self::one())
                }

                match e.cmp(&Self::exp_bias())
                {
                    Ordering::Less => Ordering::Less,
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Equal => if self.frac_bits().is_zero()
                    {
                        Ordering::Equal
                    }
                    else
                    {
                        Ordering::Greater
                    }
                }
            }
        }
    }

    pub(crate) fn cmp_one(self) -> Option<Ordering>
    {
        match self.classify()
        {
            FpCategory::Nan => None,
            FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero => Some(self.total_cmp_one())
        }
    }

    pub(crate) fn abs_cmp_one(self) -> Option<Ordering>
    {
        match self.classify()
        {
            FpCategory::Nan => None,
            FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero => Some(self.abs_total_cmp_one())
        }
    }

    pub(crate) fn abs_is_one(self) -> bool
    {
        matches!(self.abs_cmp_one(), Some(Ordering::Equal))
    }
    
    /// Returns the multiplicative identity element of `Self`, `1`.
    #[must_use]
    pub fn one() -> Self
    {
        let bias = Self::exp_bias();
        let one = if !Self::IS_INT_IMPLICIT
        {
            Self::from_bits(Self::shift_exp(bias) | Self::shift_int(U::one()))
        }
        else
        {
            Self::from_bits(Self::shift_exp(bias))
        };
        assert!(one.is_one(), "1 != 1");
        one
    }
    
    /// Sets `self` to the multiplicative identity element of `Self`, `1`.
    fn set_one(&mut self)
    {
        *self = Self::one()
    }

    /// Returns `true` if `self` is equal to the multiplicative identity.
    #[must_use = "this returns the result of the operation, without modifying the original"]
    pub fn is_one(self) -> bool
    {
        matches!(self.cmp_one(), Some(Ordering::Equal))
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> One for Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn one() -> Self
    {
        Self::one()
    }

    #[inline]
    fn is_one(&self) -> bool
    {
        (*self).is_one()
    }

    #[inline]
    fn set_one(&mut self)
    {
        self.set_one()
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> One for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn one() -> Self
    {
        Self::one()
    }

    #[inline]
    fn is_one(&self) -> bool
    {
        (*self).is_one()
    }

    #[inline]
    fn set_one(&mut self)
    {
        self.set_one()
    }
}