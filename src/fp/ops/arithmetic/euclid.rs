use num_traits::Euclid;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Calculates Euclidean division, the matching method for `rem_euclid`.
    ///
    /// This computes the integer `n` such that
    /// `self = n * rhs + self.rem_euclid(rhs)`.
    /// In other words, the result is `self / rhs` rounded to the integer `n`
    /// such that `self >= n * rhs`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let a = FpDouble::from(7.0);
    /// let b = FpDouble::from(4.0);
    /// assert_eq!(a.div_euclid(b), FpDouble::from(1.0)); // 7.0 > 4.0 * 1.0
    /// assert_eq!((-a).div_euclid(b), FpDouble::from(-2.0)); // -7.0 >= 4.0 * -2.0
    /// assert_eq!(a.div_euclid(-b), FpDouble::from(-1.0)); // 7.0 >= -4.0 * -1.0
    /// assert_eq!((-a).div_euclid(-b), FpDouble::from(2.0)); // -7.0 >= -4.0 * 2.0
    /// ```
    #[doc(alias = "modulo", alias = "mod")]
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn div_euclid(self, rhs: Self) -> Self
    {
        // TODO: do this better
        let q = (self / rhs).trunc();
        if (self % rhs).is_nonzero_negative()
        {
            let s = rhs.is_nonzero_negative();
            return (q.xor_sign(s) - U::one()).xor_sign(s)
        }
        q
    }

    /// Calculates the least nonnegative remainder of `self (mod rhs)`.
    ///
    /// In particular, the return value `r` satisfies `0.0 <= r < rhs.abs()` in
    /// most cases. However, due to a floating point round-off error it can
    /// result in `r == rhs.abs()`, violating the mathematical definition, if
    /// `self` is much smaller than `rhs.abs()` in magnitude and `self < 0.0`.
    /// This result is not an element of the function's codomain, but it is the
    /// closest floating point number in the real numbers and thus fulfills the
    /// property `self == self.div_euclid(rhs) * rhs + self.rem_euclid(rhs)`
    /// approximately.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let a = FpDouble::from(7.0);
    /// let b = FpDouble::from(4.0);
    /// assert_eq!(a.rem_euclid(b), FpDouble::from(3.0));
    /// assert_eq!((-a).rem_euclid(b), FpDouble::from(1.0));
    /// assert_eq!(a.rem_euclid(-b), FpDouble::from(3.0));
    /// assert_eq!((-a).rem_euclid(-b), FpDouble::from(1.0));
    /// // limitation due to round-off error
    /// assert!((-FpDouble::epsilon()).rem_euclid(FpDouble::from(3.0)) != FpDouble::from(0.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn rem_euclid(self, rhs: Self) -> Self
    {
        // TODO: do this better
        let r = self % rhs;
        if r.is_nonzero_negative() {r + rhs.abs()} else {r}
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Euclid for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn div_euclid(&self, rhs: &Self) -> Self
    {
        (*self).div_euclid(*rhs)
    }

    #[inline]
    fn rem_euclid(&self, rhs: &Self) -> Self
    {
        (*self).rem_euclid(*rhs)
    }
}