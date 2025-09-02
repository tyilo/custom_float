use core::num::FpCategory;

use num_traits::NumCast;

use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Compute the distance between the origin and a point (`x`, `y`) on the
    /// Euclidean plane. Equivalently, compute the length of the hypotenuse of a
    /// right-angle triangle with other sides having length `x.abs()` and
    /// `y.abs()`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(2.0);
    /// let y = FpDouble::from(3.0);
    ///
    /// // sqrt(x^2 + y^2)
    /// let abs_difference = (x.hypot(y) - (x.powi(2) + y.powi(2)).sqrt()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn hypot(self, other: Self) -> Self
    {
        let mut x = self.abs();
        let mut y = other.abs();

        match (x.classify(), y.classify())
        {
            (FpCategory::Nan, _) | (_, FpCategory::Nan) => self.add_nan(other),
            (FpCategory::Infinite, _) | (_, FpCategory::Zero) => x,
            (_, FpCategory::Infinite) | (FpCategory::Zero, _) => y,
            (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => {
                /* arrange |x| >= |y| */
                if x < y
                {
                    core::mem::swap(&mut x, &mut y)
                }
            
                /* special cases */
                let ex = x.exp_bits();
                let ey = y.exp_bits();
                /* note: hypot(x,y) ~= x + y*y/x/2 with inexact for small y/x */
                if match <usize as NumCast>::from(ex - ey)
                {
                    Some(de) => de > util::exp2_ilog(FRAC_SIZE + 1, EXP_BASE),
                    None => true
                }
                {
                    return x + y;
                }
        
                let bias = Self::exp_bias();
                let hi = Self::from(bias.to_f64().unwrap()*0.7).expb();
                let lo = Self::from(-bias.to_f64().unwrap()*0.7).expb();
            
                /* precise sqrt argument in nearest rounding mode without overflow */
                /* xh*xh must not overflow and xl*xl must not underflow in sq */
                let bias_half = bias/U::from(2).unwrap();
                let mut z = Self::one();
                if ex > bias + bias_half
                {
                    z = hi;
                    x *= lo;
                    y *= lo;
                }
                else if ey < bias - bias_half
                {
                    z = lo;
                    x *= hi;
                    y *= hi;
                }
                z * (x*x + y*y).sqrt()
            }
        }
    }
}