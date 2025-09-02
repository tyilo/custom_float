use core::num::FpCategory;

use crate::{fp::Fps, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn classify(self) -> FpCategory
    {
        self.0.classify()
    }
    
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_finite(self) -> bool
    {
        self.0.is_finite()
    }

    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_normal(self) -> bool
    {
        self.0.is_normal()
    }

    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_subnormal(self) -> bool
    {
        self.0.is_subnormal()
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the floating point category of the number. If only one property
    /// is going to be tested, it is generally faster to use the specific
    /// predicate instead.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use std::num::FpCategory;
    ///
    /// let num = FpDouble::from(12.4f32);
    /// let inf = FpDouble::infinity();
    ///
    /// assert_eq!(num.classify(), FpCategory::Normal);
    /// assert_eq!(inf.classify(), FpCategory::Infinite);
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn classify(self) -> FpCategory
    {
        enum ExpCategory
        {
            Zero,
            Max,
            Normal
        }

        enum MantissaCategory
        {
            Zero,
            Nonzero
        }

        let bits = self.to_bits();
        let mant_mask = Self::mantissa_mask_raw();
        match (
            if EXP_SIZE == 0
            {
                ExpCategory::Normal
            }
            else
            {
                let exp_mask = Self::exponent_mask();
                let e = bits & exp_mask;
                match (e.is_zero(), e == exp_mask)
                {
                    (false, false) => ExpCategory::Normal,
                    (true, false) => ExpCategory::Zero,
                    (false, true) => ExpCategory::Max,
                    (true, true) => ExpCategory::Normal //What to do here?
                }
            },
            if (bits & mant_mask).is_zero()
            {
                MantissaCategory::Zero
            }
            else
            {
                MantissaCategory::Nonzero
            }
        )
        {
            (ExpCategory::Zero, MantissaCategory::Zero) => FpCategory::Zero,
            (ExpCategory::Zero, MantissaCategory::Nonzero) => if Self::IS_INT_IMPLICIT
            {
                FpCategory::Subnormal
            }
            else
            {
                FpCategory::Normal
            },
            (ExpCategory::Normal, MantissaCategory::Nonzero | MantissaCategory::Zero) => FpCategory::Normal,
            (ExpCategory::Max, MantissaCategory::Zero) => FpCategory::Infinite,
            (ExpCategory::Max, MantissaCategory::Nonzero) => FpCategory::Nan
        }
    }
    
    /// Returns `true` if this number is neither infinite nor NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let f = FpSingle::from(7.0);
    /// let inf = FpSingle::infinity();
    /// let neg_inf = FpSingle::neg_infinity();
    /// let nan = FpSingle::nan();
    ///
    /// assert!(f.is_finite());
    ///
    /// assert!(!nan.is_finite());
    /// assert!(!inf.is_finite());
    /// assert!(!neg_inf.is_finite());
    /// ```
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_finite(self) -> bool
    {
        matches!(self.classify(), FpCategory::Zero | FpCategory::Subnormal | FpCategory::Normal)
    }

    /// Returns `true` if the number is neither zero, infinite,
    /// [subnormal], or NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let min = FpSingle::min_positive_value(); // 1.17549435e-38f32
    /// let max = FpSingle::max_value();
    /// let lower_than_min = FpSingle::from(1.0e-40_f32);
    /// let zero = FpSingle::zero();
    ///
    /// assert!(min.is_normal());
    /// assert!(max.is_normal());
    ///
    /// assert!(!zero.is_normal());
    /// assert!(!FpSingle::nan().is_normal());
    /// assert!(!FpSingle::infinity().is_normal());
    /// // Values between `0` and `min` are Subnormal.
    /// assert!(!lower_than_min.is_normal());
    /// ```
    /// [subnormal]: http://en.wikipedia.org/wiki/Denormal_number
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_normal(self) -> bool
    {
        matches!(self.classify(), FpCategory::Normal)
    }

    /// Returns `true` if the number is [subnormal].
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let min = FpDouble::min_positive_value(); // 2.2250738585072014e-308_f64
    /// let max = FpDouble::max_value();
    /// let lower_than_min = FpDouble::from(1.0e-308_f64);
    /// let zero = FpDouble::zero();
    ///
    /// assert!(!min.is_subnormal());
    /// assert!(!max.is_subnormal());
    ///
    /// assert!(!zero.is_subnormal());
    /// assert!(!FpDouble::nan().is_subnormal());
    /// assert!(!FpDouble::infinity().is_subnormal());
    /// // Values between `0` and `min` are Subnormal.
    /// assert!(lower_than_min.is_subnormal());
    /// ```
    /// [subnormal]: https://en.wikipedia.org/wiki/Denormal_number
    #[must_use = "this returns the result of the operation, without modifying the original"]
    #[inline]
    pub fn is_subnormal(self) -> bool
    {
        matches!(self.classify(), FpCategory::Subnormal)
    }
}