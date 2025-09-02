use core::{cmp::Ordering, num::FpCategory};

use num_traits::FloatConst;

use crate::{util, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn cos_extra_sign(self) -> (Self, bool)
    {
        match self.classify()
        {
            FpCategory::Nan => (self, false),
            FpCategory::Infinite => (Self::nan(), false),
            FpCategory::Zero | FpCategory::Normal | FpCategory::Subnormal => {
                const N: usize = 6;
                const C: [f64; N] = [
                    0.472001216,
                    -0.499403258,
                    0.027992080,
                    -0.000596695,
                    0.000006704,
                    -0.000000047
                ];
        
                static mut P: Option<[f64; N]> = None;
                let p = unsafe {
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
                w = w.abs();
                let mut s = matches!(w.abs_cmp_one(), Some(Ordering::Greater));
                if s
                {
                    let b;
                    (w, b) = w.add_int_with_sign_extra_sign(true, two, false);
                    s ^= b;
                }
        
                let ww2 = w*w*two;
                let y = {
                    let (mut z, mut z_s) = ww2.sub_int_extra_sign(one);
        
                    (z, z_s) = z.polynomial(&p, z_s);
                    s ^= z_s;
        
                    z
                };
        
                (y, s)
            }
        }
    }

    /// Computes the cosine of a number (in radians).
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
    /// let x = FpDouble::TAU();
    ///
    /// let abs_difference = (x.cos() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn cos(self) -> Self
    {
        let (y, s) = self.cos_extra_sign();
        y.xor_sign(s)
    }
}