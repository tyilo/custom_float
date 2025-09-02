use core::{cmp::Ordering, num::FpCategory};

use num_traits::FloatConst;

use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn sin_extra_sign(self) -> (Self, bool)
    {
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Zero => (self, false),
            FpCategory::Infinite => (Self::nan(), false),
            FpCategory::Normal | FpCategory::Subnormal => {
                if matches!(self.abs_partial_cmp(Self::from(0.000244140625)), None | Some(Ordering::Less))
                {
                    return (self, false)
                }
        
                const N: usize = 6;
                const C: [f64; N] = [
                    1.276278962,
                    -0.285261569,
                    0.009118016,
                    -0.000136587,
                    0.000001185,
                    -0.000000007
                ];
        
                static mut P: Option<[f64; N]> = None;
                let mut p = unsafe {
                    #[allow(static_mut_refs)]
                    *P.get_or_insert_with(|| util::chebychev_approximation(C))
                };
        
                let one = U::one();
                let two = one + one;
                let three = two + one;
                let four = two + two;

                let mut w = self*Self::FRAC_2_PI();
                w %= four;
                for i in 0..4
                {
                    let r = i % 2 == 0 && matches!(w.abs_cmp_one(), Some(Ordering::Less));
                    w += three;
                    if r
                    {
                        w = w.add_int_with_sign(true, four, false);
                    }
                    w %= four;
                }
                let mut s = false;
                if let Some(b) = match w.abs_cmp_one()
                {
                    None => return (Self::nan(), false),
                    Some(Ordering::Greater) => Some(false),
                    Some(Ordering::Less) => Some(true),
                    Some(Ordering::Equal) => None
                }
                {
                    (w, s) = w.add_int_with_sign_extra_sign(true, two, b ^ w.is_sign_negative());
                }
                
                if matches!(w.abs_partial_cmp(Self::from(1.5542474911317903883680055016847e-4)), None | Some(Ordering::Less))
                {
                    return (w*Self::FRAC_PI_2(), false)
                }
        
                let ww2 = w*w*two;
                let y = {
                    let (mut z, mut z_s) = ww2.sub_int_extra_sign(one);
        
                    (z, z_s) = z.polynomial(&p, z_s);
                    s ^= z_s;

                    w*z
                };
        
                (y, s)
            }
        }
    }

    /// Computes the sine of a number (in radians).
    ///
    /// This implementation is based on Harvey M. Wagner's [Polynomial approximations to elementary functions](https://www.ams.org/journals/mcom/1954-08-047/S0025-5718-1954-0063487-2/S0025-5718-1954-0063487-2.pdf).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::FRAC_PI_2();
    ///
    /// let abs_difference = (x.sin() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn sin(self) -> Self
    {
        let (mut y, s) = self.sin_extra_sign();
        if s
        {
            y = -y
        }
        y
    }
}