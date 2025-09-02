use core::{cmp::Ordering, num::FpCategory};

use num_traits::float::TotalOrder;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Return the ordering between `self` and `other`.
    ///
    /// Unlike the standard partial comparison between floating point numbers,
    /// this comparison always produces an ordering in accordance to
    /// the `totalOrder` predicate as defined in the IEEE 754 (2008 revision)
    /// floating point standard. The values are ordered in the following sequence:
    ///
    /// - negative quiet NaN
    /// - negative signaling NaN
    /// - negative infinity
    /// - negative numbers
    /// - negative subnormal numbers
    /// - negative zero
    /// - positive zero
    /// - positive subnormal numbers
    /// - positive numbers
    /// - positive infinity
    /// - positive signaling NaN
    /// - positive quiet NaN.
    ///
    /// The ordering established by this function does not always agree with the
    /// [`PartialOrd`] and [`PartialEq`] implementations. For example,
    /// they consider negative and positive zero equal, while `total_cmp`
    /// doesn't.
    ///
    /// The interpretation of the signaling NaN bit follows the definition in
    /// the IEEE 754 standard, which may not match the interpretation by some of
    /// the older, non-conformant (e.g. MIPS) hardware implementations.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::{FpSingle, FpDouble};
    /// use std::cmp::Ordering;
    ///
    /// assert_eq!(FpDouble::nan().total_cmp(FpDouble::nan()), Ordering::Equal);
    /// assert_eq!(FpSingle::nan().total_cmp(FpSingle::nan()), Ordering::Equal);
    ///
    /// assert_eq!((-FpDouble::nan()).total_cmp(FpDouble::nan()), Ordering::Less);
    /// assert_eq!(FpDouble::infinity().total_cmp(FpDouble::nan()), Ordering::Less);
    /// assert_eq!((-FpDouble::zero()).total_cmp(FpDouble::zero()), Ordering::Less);
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    pub fn total_cmp(self, other: Self) -> Ordering
    {
        let s = self.is_sign_negative();
        let s1 = other.is_sign_negative();

        if s != s1
        {
            return s1.cmp(&s)
        }

        let mut cmp = self.abs_total_cmp(other);
        if s
        {
            cmp = cmp.reverse()
        }
        cmp
    }

    pub fn abs_total_cmp(self, other: Self) -> Ordering
    {
        let c = |x: Self| match x.classify()
        {
            FpCategory::Nan => 2 + x.is_snan() as u8,
            FpCategory::Infinite => 1,
            FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero => 0
        };
        
        let c0 = c(self);
        let c1 = c(other);

        if (c0, c1) != (0, 0)
        {
            return c0.cmp(&c1)
        }

        let mut e0 = self.exp_bits();
        let mut e1 = other.exp_bits();

        let mut f0 = self.mantissa_bits_raw();
        let mut f1 = other.mantissa_bits_raw();

        if !Self::IS_INT_IMPLICIT
        {
            Self::normalize_mantissa_down(&mut e0, &mut f0, Some(e1));
            Self::normalize_mantissa_down(&mut e1, &mut f1, Some(e0));
            Self::normalize_mantissa_up(&mut e0, &mut f0, Some(e1));
            Self::normalize_mantissa_up(&mut e1, &mut f1, Some(e0));
        }
        
        let (a, b) = if e0 != e1
        {
            (e0, e1)
        }
        else
        {
            (f0, f1)
        };
        a.cmp(&b)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> TotalOrder for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn total_cmp(&self, other: &Self) -> Ordering
    {
        (*self).total_cmp(*other)
    }
}