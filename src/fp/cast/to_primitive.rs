use core::{cmp::Ordering, num::FpCategory, ops::Neg};

use num_traits::ToPrimitive;

use crate::{util::{self, Saturation}, AnyInt, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn checked_to_int<I: AnyInt>(self, convert: impl FnOnce(U, &mut bool) -> Result<I, Saturation>) -> Result<I, Option<Saturation>>
    {
        let mut s = self.is_sign_negative();
        let apply_sign = |sgn| if sgn {Neg::neg} else {util::do_nothing};
        match self.classify()
        {
            FpCategory::Zero => Ok(I::zero()),
            FpCategory::Nan => Err(None),
            FpCategory::Infinite => Err(Some(apply_sign(s)(Saturation::Overflow))),
            FpCategory::Normal | FpCategory::Subnormal => if !util::is_signed::<I>() && s
            {
                Err(Some(Saturation::Underflow))
            }
            else
            {
                let mut e = self.exp_bits();
                let mut f = self.mantissa_bits();
                let mut n;
                let bias = Self::exp_bias();

                let check = |result| match result
                {
                    Ok(()) => None,
                    Err(sat) => Some(
                        match sat
                        {
                            Saturation::Overflow => Err(Some(apply_sign(s)(Saturation::Overflow))),
                            Saturation::Underflow => Ok(I::zero())
                        }
                    )
                };
        
                if util::bitsize_of::<I>() - util::is_signed::<I>() as usize > util::bitsize_of::<U>()
                {
                    n = I::from(f).unwrap();
                    if s
                    {
                        assert!(util::is_signed::<I>());
                        n = util::neg(n);
                    }
        
                    if let Some(done) = check(
                        Self::mantissa_shr(&mut n, FRAC_SIZE, &mut e)
                            .and_then(|()| Self::force_exp_to(&mut e, &mut n, bias))
                    )
                    {
                        return done
                    }
                }
                else
                {
                    if let Some(done) = check(
                        Self::mantissa_shr(&mut f, FRAC_SIZE, &mut e)
                            .and_then(|()| Self::force_exp_to(&mut e, &mut f, bias))
                    )
                    {
                        return done
                    }

                    if s && f == U::one() << (util::bitsize_of::<I>() - 1)
                    {
                        assert!(util::is_signed::<I>());
                        return Ok(I::min_value())
                    }
        
                    n = convert(f, &mut s).map_err(apply_sign(s))?;
                    if s
                    {
                        assert!(util::is_signed::<I>());
                        n = util::neg(n);
                    }
                }
        
                Ok(n)
            }
        }
    }
    
    /// Converts a custom floating-point type into an integer.
    ///
    /// Returns [None](None) if out of bounds.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn to_int<I: AnyInt>(self) -> Option<I>
    {
        self.checked_to_int(|x, _| I::from(x).ok_or(Saturation::Overflow)).ok()
    }
    
    /// Converts a custom floating-point type into an integer.
    ///
    /// Wraps if out of bounds.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn to_int_wrapping<I: AnyInt>(mut self) -> I
    {
        let s = !util::is_signed::<I>() && self.is_sign_negative();
        if !util::is_signed::<I>()
        {
            self = self.abs()
        }
        let mut y = match self.checked_to_int(|mut x, s| {
            let mut iwrap = U::zero();
            if core::mem::size_of::<U>() > core::mem::size_of::<I>()
            {
                iwrap = U::one() << (util::bitsize_of::<I>());
                x = x % iwrap;
            }
            if util::is_signed::<I>() && core::mem::size_of::<U>() >= core::mem::size_of::<I>()
            {
                while match x.cmp(&(U::one() << (util::bitsize_of::<I>() - 1)))
                {
                    Ordering::Less => false,
                    Ordering::Equal => !*s,
                    Ordering::Greater => true
                }
                {
                    x = iwrap.wrapping_sub(&x);
                    *s = !*s;
                }
            }
            I::from(x).ok_or(Saturation::Overflow)
        })
        {
            Ok(y) => y,
            Err(err) => match err
            {
                Some(sat) => match sat
                {
                    Saturation::Overflow => I::max_value(),
                    Saturation::Underflow => I::min_value(),
                },
                None => I::zero(),
            },
        };
        if s
        {
            assert!(!util::is_signed::<I>());
            y = (I::max_value() - y).wrapping_add(&I::one());
        }
        y
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> ToPrimitive for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn to_i8(&self) -> Option<i8>
    {
        self.to_int()
    }
    #[inline]
    fn to_u8(&self) -> Option<u8>
    {
        self.to_int()
    }

    #[inline]
    fn to_i16(&self) -> Option<i16>
    {
        self.to_int()
    }
    #[inline]
    fn to_u16(&self) -> Option<u16>
    {
        self.to_int()
    }

    #[inline]
    fn to_i32(&self) -> Option<i32>
    {
        self.to_int()
    }
    #[inline]
    fn to_u32(&self) -> Option<u32>
    {
        self.to_int()
    }

    #[inline]
    fn to_isize(&self) -> Option<isize>
    {
        self.to_int()
    }
    #[inline]
    fn to_usize(&self) -> Option<usize>
    {
        self.to_int()
    }

    #[inline]
    fn to_i64(&self) -> Option<i64>
    {
        self.to_int()
    }
    #[inline]
    fn to_u64(&self) -> Option<u64>
    {
        self.to_int()
    }

    #[inline]
    fn to_i128(&self) -> Option<i128>
    {
        self.to_int()
    }
    #[inline]
    fn to_u128(&self) -> Option<u128>
    {
        self.to_int()
    }

    #[inline]
    fn to_f32(&self) -> Option<f32>
    {
        Some(Into::<f32>::into(*self))
    }
    
    #[inline]
    fn to_f64(&self) -> Option<f64>
    {
        Some(Into::<f64>::into(*self))
    }
}