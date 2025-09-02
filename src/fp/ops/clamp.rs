use core::cmp::Ordering;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the maximum of the two numbers.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.0);
    /// let y = FpDouble::from(2.0);
    ///
    /// assert_eq!(x.max(y), y);
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    #[inline]
    pub fn max(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    self.add_nan(other)
                }
                else if s1
                {
                    other
                }
                else
                {
                    self
                }
            },
            (true, false) => other,
            (false, true) => self,
            (false, false) => match self.total_cmp(other)
            {
                Ordering::Equal => if self.is_sign_negative() && other.is_sign_positive()
                {
                    other
                }
                else
                {
                    self
                },
                Ordering::Greater => self,
                Ordering::Less => other
            }
        }
    }

    /// Returns the minimum of the two numbers.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.0);
    /// let y = FpDouble::from(2.0);
    ///
    /// assert_eq!(x.min(y), x);
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    #[inline]
    pub fn min(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    self.add_nan(other)
                }
                else if s2
                {
                    other
                }
                else
                {
                    self
                }
            },
            (true, false) => other,
            (false, true) => self,
            (false, false) => match self.total_cmp(other)
            {
                Ordering::Equal => if self.is_sign_positive() && other.is_sign_negative()
                {
                    other
                }
                else
                {
                    self
                },
                Ordering::Greater => other,
                Ordering::Less => self
            }
        }
    }

    /// Returns the maximum of the two numbers, propagating NaN.
    ///
    /// This returns NaN when *either* argument is NaN, as opposed to
    /// [`Fp::max`] which only returns NaN when *both* arguments are NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// let x = FpSingle::from(1.0);
    /// let y = FpSingle::from(2.0);
    ///
    /// assert_eq!(x.maximum(y), y);
    /// assert!(x.maximum(FpSingle::nan()).is_nan());
    /// ```
    ///
    /// If one of the arguments is NaN, then NaN is returned. Otherwise this returns the greater
    /// of the two numbers. For this operation, -0.0 is considered to be less than +0.0.
    /// Note that this follows the semantics specified in IEEE 754-2019.
    ///
    /// Also note that "propagation" of NaNs here doesn't necessarily mean that the bitpattern of a NaN
    /// operand is conserved; see [explanation of NaN as a special value](f32) for more info.
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    #[inline]
    pub fn maximum(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    self.add_nan(other)
                }
                else if s1
                {
                    other
                }
                else
                {
                    self
                }
            },
            (true, false) | (false, true) => self.add_nan(other),
            (false, false) => match self.total_cmp(other)
            {
                Ordering::Equal => if self.is_sign_negative() && other.is_sign_positive()
                {
                    other
                }
                else
                {
                    self
                },
                Ordering::Greater => self,
                Ordering::Less => other
            }
        }
    }
    
    /// Returns the minimum of the two numbers, propagating NaN.
    ///
    /// This returns NaN when *either* argument is NaN, as opposed to
    /// [`Fp::min`] which only returns NaN when *both* arguments are NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// let x = FpSingle::from(1.0);
    /// let y = FpSingle::from(2.0);
    ///
    /// assert_eq!(x.minimum(y), x);
    /// assert!(x.minimum(FpSingle::nan()).is_nan());
    /// ```
    ///
    /// If one of the arguments is NaN, then NaN is returned. Otherwise this returns the lesser
    /// of the two numbers. For this operation, -0.0 is considered to be less than +0.0.
    /// Note that this follows the semantics specified in IEEE 754-2019.
    ///
    /// Also note that "propagation" of NaNs here doesn't necessarily mean that the bitpattern of a NaN
    /// operand is conserved; see [explanation of NaN as a special value](f32) for more info.
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    #[inline]
    pub fn minimum(self, other: Self) -> Self
    {
        match (self.is_nan(), other.is_nan())
        {
            (true, true) => {
                let s1 = self.is_sign_negative();
                let s2 = other.is_sign_negative();
                if s1 == s2
                {
                    self.add_nan(other)
                }
                else if s2
                {
                    other
                }
                else
                {
                    self
                }
            },
            (true, false) | (false, true) => self.add_nan(other),
            (false, false) => match self.total_cmp(other)
            {
                Ordering::Equal => if self.is_sign_positive() && other.is_sign_negative()
                {
                    other
                }
                else
                {
                    self
                },
                Ordering::Greater => other,
                Ordering::Less => self
            }
        }
    }

    /// Restrict a value to a certain interval unless it is NaN.
    ///
    /// Returns `max` if `self` is greater than `max`, and `min` if `self` is
    /// less than `min`. Otherwise this returns `self`.
    ///
    /// Note that this function returns NaN if the initial value was NaN as
    /// well.
    ///
    /// # Panics
    ///
    /// Panics if `min > max`, `min` is NaN, or `max` is NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let min = FpSingle::from(-2.0);
    /// let max = FpSingle::one();
    ///
    /// assert_eq!(FpSingle::from(-3.0).clamp(min, max), FpSingle::from(-2.0));
    /// assert_eq!(FpSingle::zero().clamp(min, max), FpSingle::zero());
    /// assert_eq!(FpSingle::from(2.0).clamp(min, max), FpSingle::one());
    /// assert!(FpSingle::nan().clamp(min, max).is_nan());
    /// ```
    #[must_use = "this returns the result of the comparison, without modifying either input"]
    pub fn clamp(mut self, min: Self, max: Self) -> Self
    {
        assert!(min <= max, "min > max, or either was NaN. min = {min:?}, max = {max:?}");
        if self < min
        {
            self = min;
        }
        if self > max
        {
            self = max;
        }
        self
    }
}