use core::num::FpCategory;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the least number greater than `self`.
    ///
    /// Let `TINY` be the smallest representable positive value. Then,
    ///  - if `self.is_nan()`, this returns `self`;
    ///  - if `self` is [`NEG_INFINITY`], this returns [`MIN`];
    ///  - if `self` is `-TINY`, this returns -0.0;
    ///  - if `self` is -0.0 or +0.0, this returns `TINY`;
    ///  - if `self` is [`MAX`] or [`INFINITY`], this returns [`INFINITY`];
    ///  - otherwise the unique least value greater than `self` is returned.
    ///
    /// 
    /// The identity `x.next_up() == -(-x).next_down()` holds for all non-NaN `x`.
    /// When `x` is finite and the radix is 2, `x == x.next_up().next_down()` also holds.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    /// 
    /// // epsilon is the difference between 1.0 and the next number up.
    /// assert_eq!(FpSingle::one().next_up(), FpSingle::one() + FpSingle::epsilon());
    /// // But not for most numbers.
    /// assert!(FpSingle::from(0.1).next_up() < FpSingle::from(0.1) + FpSingle::epsilon());
    /// assert_eq!(FpSingle::from(16777216.0).next_up(), FpSingle::from(16777218.0));
    /// ```
    ///
    /// [`NEG_INFINITY`]: Self::neg_infinity
    /// [`INFINITY`]: Self::infinity
    /// [`MIN`]: Self::min_value
    /// [`MAX`]: Self::max_value
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn next_up(self) -> Self
    {
        self.next(false)
    }
    
    /// Returns the greatest number less than `self`.
    ///
    /// Let `TINY` be the smallest representable positive value. Then,
    ///  - if `self.is_nan()`, this returns `self`;
    ///  - if `self` is [`INFINITY`], this returns [`MAX`];
    ///  - if `self` is `TINY`, this returns 0.0;
    ///  - if `self` is -0.0 or +0.0, this returns `-TINY`;
    ///  - if `self` is [`MIN`] or [`NEG_INFINITY`], this returns [`NEG_INFINITY`];
    ///  - otherwise the unique greatest value less than `self` is returned.
    ///
    /// The identity `x.next_down() == -(-x).next_up()` holds for all non-NaN `x`. When `x`
    /// is finite `x == x.next_down().next_up()` also holds.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let x = FpSingle::one();
    /// // Clamp value into range [0, 1).
    /// let clamped = x.clamp(FpSingle::zero(), FpSingle::one().next_down());
    /// assert!(clamped < FpSingle::one());
    /// assert_eq!(clamped.next_up(), FpSingle::one());
    /// ```
    ///
    /// [`NEG_INFINITY`]: Self::neg_infinity
    /// [`INFINITY`]: Self::infinity
    /// [`MIN`]: Self::min_value
    /// [`MAX`]: Self::max_value
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn next_down(self) -> Self
    {
        self.next(true)
    }

    fn next(self, sign: bool) -> Self
    {
        let mut s = self.is_sign_negative();
        match self.classify()
        {
            FpCategory::Nan => self,
            FpCategory::Infinite => if s == sign
            {
                self
            }
            else
            {
                Self::max_value().with_sign(s)
            },
            FpCategory::Zero | FpCategory::Subnormal | FpCategory::Normal => {
                let mut e = self.exp_bits();
        
                if Self::MANTISSA_DIGITS == 0
                {
                    if s ^ sign
                    {
                        if e.is_zero()
                        {
                            return -self
                        }
                        else
                        {
                            e = e - U::one()
                        }
                    }
                    else
                    {
                        e = e + U::one()
                    };
                    let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
                    Self::from_bits(s_bits + Self::shift_exp(e))
                }
                else if Self::IS_INT_IMPLICIT
                {
                    let mut f = self.frac_bits();
        
                    if s ^ sign
                    {
                        if f.is_zero()
                        {
                            if e.is_zero()
                            {
                                s = !s;
                            }
                            else
                            {
                                e = e - U::one();
                                f = Self::max_frac_bits()
                            }
                        }
                        else
                        {
                            f = f - U::one()
                        }
                    }
                    else if f == Self::max_frac_bits()
                    {
                        e = e + U::one();
                        f = U::zero()
                    }
                    else
                    {
                        f = f + U::one()
                    }
        
                    let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
                    Self::from_bits(s_bits + Self::shift_exp(e) + Self::shift_frac(f))
                }
                else
                {
                    let mut f = self.mantissa_bits_raw();
        
                    let base = U::from(EXP_BASE);
        
                    if s ^ sign
                    {
                        if f.is_zero()
                        {
                            s = !s;
                            e = U::zero()
                        }
                        else if let Some(base) = base && let Some(mut ff) = f.checked_mul(&base)
                        {
                            let mut o = U::one();
                            while !e.is_zero() && let oo = ff.saturating_sub(Self::max_mantissa_bits()) && oo < base
                            {
                                o = oo.max(U::one());
                                e = e - U::one();
                                f = ff;
                                if let Some(fff) = f.checked_mul(&base)
                                {
                                    ff = fff
                                }
                                else
                                {
                                    break
                                }
                            }
                            f = f - o
                        }
                        else if !e.is_zero() && f == U::one()
                        {
                            e = e - U::one();
                            f = Self::max_mantissa_bits()
                        }
                        else
                        {
                            f = f - U::one()
                        }
                    }
                    else
                    {
                        if let Some(base) = base && let Some(mut ff) = f.checked_mul(&base)
                        {
                            while !e.is_zero() && ff <= Self::max_mantissa_bits()
                            {
                                e = e - U::one();
                                f = ff;
                                if let Some(fff) = f.checked_mul(&base)
                                {
                                    ff = fff
                                }
                                else
                                {
                                    break
                                }
                            }
                        }
                        if f == Self::max_mantissa_bits()
                        {
                            e = e + U::one();
                            if e == Self::max_exponent_bits()
                            {
                                let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
                                return Self::from_bits(s_bits + Self::shift_exp(e))
                            }
                            if let Some(base) = base
                            {
                                f = f/base + U::one();
                            }
                            else
                            {
                                f = U::one()
                            }
                        }
                        else
                        {
                            f = f + U::one()
                        }
                    }
        
                    let s_bits = if s {Self::shift_sign(U::one())} else {U::zero()};
                    Self::from_bits(s_bits + Self::shift_exp(e) + Self::shift_frac(f))
                }
            }
        }
    }
}