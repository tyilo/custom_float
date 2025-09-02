use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Calculates the middle point of `self` and `rhs`.
    ///
    /// This returns NaN when *either* argument is NaN or if a combination of
    /// +inf and -inf is provided as arguments.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// assert_eq!(FpSingle::from(1.0).midpoint(FpSingle::from(4.0)), FpSingle::from(2.5));
    /// assert_eq!(FpSingle::from(-5.5).midpoint(FpSingle::from(8.0)), FpSingle::from(1.25));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn midpoint(self, rhs: Self) -> Self
    {
        let two = U::from(2).unwrap();

        let lo = Self::min_positive_value()*two;
        let hi = Self::max_value()/two;

        let (a, b) = (self, rhs);
        let abs_a = a.abs();
        let abs_b = b.abs();

        if abs_a <= hi && abs_b <= hi
        {
            // Overflow is impossible
            (a + b)/two
        }
        else if abs_a < lo
        {
            // Not safe to halve a
            a + b/two
        }
        else if abs_b < lo
        {
            // Not safe to halve b
            a/two + b
        }
        else
        {
            // Not safe to halve a and b
            a/two + b/two
        }
    }
}