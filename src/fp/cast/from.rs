use core::{cmp::Ordering, num::FpCategory};

use num_traits::{NumCast, Zero};

use crate::{fp::Fps, ieee754::{FpDouble, FpHalf, FpQuadruple, FpSingle}, util::{self, Saturation}, AnyInt, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Converts to this type from the input type.
    #[inline]
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from<T>(from: T) -> Self
    where
        Self: From<T>
    {
        <Self as From<T>>::from(from)
    }

    /// Converts from one custom floating-point number to another.
    /// Rounding errors may occurr.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_fp<V, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize>(fp: Fp<V, S, E, I, F, B>) -> Self
    where
        V: FpRepr<S, E, I, F, B>
    {
        Fps::from_fp(fp).collect()
    }
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub(crate) fn from_fps<V, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize>(fps: Fps<V, S, E, I, F, B>) -> Self
    where
        V: FpRepr<S, E, I, F, B>
    {
        Fps::from_fps(fps).collect()
    }

    pub(crate) fn from_int_diff<I: AnyInt>(lhs: I, rhs: I) -> Self
    {
        Fps::from_int_diff(lhs, rhs).collect()
    }
    
    /// Converts an integer into a custom floating-point type.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_int<I: AnyInt>(from: I) -> Self
    {
        Fps::from_int(from).collect()
    }

    #[must_use = "method returns a new number and does not mutate the original value"]
    fn abs_from_fp<V, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize>(mut fp: Fp<V, S, E, I, F, B>) -> Self
    where
        V: FpRepr<S, E, I, F, B>
    {
        if EXP_SIZE == E && INT_SIZE == I && EXP_BASE == B && (I == 0 || FRAC_SIZE >= F)
        {
            fp = fp.abs();
            if let Some(b) = if FRAC_SIZE == F
            {
                <U as NumCast>::from(fp.to_bits())
            }
            else if util::bitsize_of::<U>() >= util::bitsize_of::<V>()
            {
                <U as NumCast>::from(fp.to_bits())
                    .map(|b| if FRAC_SIZE >= F
                    {
                        b << (FRAC_SIZE - F)
                    }
                    else
                    {
                        util::rounding_div_2_pow(b, F - FRAC_SIZE)
                    })
            }
            else
            {
                let b = if FRAC_SIZE >= F
                {
                    fp.to_bits() << (FRAC_SIZE - F)
                }
                else
                {
                    util::rounding_div_2_pow(fp.to_bits(), F - FRAC_SIZE)
                };
                <U as NumCast>::from(b)
            }
            {
                return Self::from_bits(b)
            }
        }

        match fp.classify()
        {
            FpCategory::Nan => {
                if fp.is_snan()
                {
                    Self::snan()
                }
                else
                {
                    Self::qnan()
                }
            },
            FpCategory::Infinite => Self::infinity(),
            FpCategory::Zero => Self::zero(),
            FpCategory::Subnormal | FpCategory::Normal => {
                let mut e1 = fp.exp_bits();
                let mut f = fp.mantissa_bits();

                let df = FRAC_SIZE as isize - F as isize;

                let base1 = V::from(B);

                let mut f = loop
                {
                    match if df >= 0
                    {
                        U::from(f).and_then(|f| if f.leading_zeros() as usize >= df as usize
                        {
                            f.checked_shl(df as u32)
                        }
                        else
                        {
                            None
                        })
                    }
                    else
                    {
                        U::from(util::rounding_div_2_pow(f, (-df) as usize))
                    }
                    {
                        Some(f) => break f,
                        None => {
                            e1 = e1 + V::one();
                            f = if let Some(base1) = base1
                            {
                                if B.is_zero()
                                {
                                    return Self::infinity()
                                }
                                util::rounding_div(f, base1)
                            }
                            else
                            {
                                V::zero()
                            }
                        }
                    }
                };
                
                let bias1 = Fp::<V, S, E, I, F, B>::exp_bias();
                let bias2 = Self::exp_bias();

                let mut e = {
                    match if EXP_BASE == B
                    {
                        if util::bitsize_of::<U>() < util::bitsize_of::<V>() && let Some(bias) = V::from(bias2)
                        {
                            match bias1.cmp(&bias)
                            {
                                Ordering::Greater => e1.checked_sub(&(bias1 - bias)),
                                Ordering::Equal => Some(e1),
                                Ordering::Less => e1.checked_add(&(bias - bias1))
                            }.and_then(U::from)
                        }
                        else if let Some(bias) = U::from(bias1) && let Some(e) = U::from(e1)
                        {
                            match bias.cmp(&bias2)
                            {
                                Ordering::Greater => e.checked_sub(&(bias - bias2)),
                                Ordering::Equal => Some(e),
                                Ordering::Less => e.checked_add(&(bias2 - bias))
                            }
                        }
                        else
                        {
                            None
                        }
                    }
                    else
                    {
                        None
                    }
                    {
                        Some(e) => e,
                        None => {
                            let base1 = U::from(B);
                            let base2 = U::from(EXP_BASE);
        
                            let mut e = bias2;
                            if let Some(base1) = base1
                            {
                                while e1 > bias1
                                {
                                    e1 = e1 - V::one();
                                    loop
                                    {
                                        if let Some(ff) = f.checked_mul(&base1)
                                        {
                                            f = ff;
                                            break
                                        }
                                        else if let Some(ee) = e.checked_add(&U::one())
                                        {
                                            e = ee;
                                            if let Some(base2) = base2
                                            {
                                                if EXP_BASE.is_zero()
                                                {
                                                    return Self::infinity()
                                                }
                                                f = util::rounding_div(f, base2)
                                            }
                                            else
                                            {
                                                f = U::zero();
                                                break
                                            }
                                        }
                                        else
                                        {
                                            return Self::infinity()
                                        }
                                    }
                                }
                            }
                            if let Some(base2) = base2
                            {
                                while e1 < bias1
                                {
                                    e1 = e1 + V::one();
                                    while let Some(ee) = e.checked_sub(&U::one())
                                        && let Some(ff) = f.checked_mul(&base2)
                                    {
                                        e = ee;
                                        f = ff;
                                    }
                                    f = if let Some(base1) = base1
                                    {
                                        if B.is_zero()
                                        {
                                            return Self::infinity()
                                        }
                                        util::rounding_div(f, base1)
                                    }
                                    else
                                    {
                                        U::zero()
                                    };
                                }
                            }
                            e
                        }
                    }
                };

                Self::normalize_mantissa(&mut e, &mut f, None);
                Self::from_exp_mantissa(e, f)
            }
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Converts to this type from the input type.
    #[inline]
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from<T>(from: T) -> Self
    where
        Self: From<T>
    {
        <Self as From<T>>::from(from)
    }

    /// Converts from one custom floating-point number to another.
    /// Rounding errors may occurr.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_fp<V, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize>(fp: Fp<V, S, E, I, F, B>) -> Self
    where
        V: FpRepr<S, E, I, F, B>
    {
        if EXP_SIZE == E && INT_SIZE == I && EXP_BASE == B && SIGN_BIT == S && FRAC_SIZE == F
        {
            if let Some(b) = <U as NumCast>::from(fp.to_bits())
            {
                return Fp::from_bits(b).extra_sign(false)
            }
        }

        let s = fp.is_sign_negative();
        Fp::abs_from_fp(fp).extra_sign(s)
    }
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_fps<V, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize>(fps: Fps<V, S, E, I, F, B>) -> Self
    where
        V: FpRepr<S, E, I, F, B>
    {
        let Fps(fp, s_x) = fps;
        Self::from_fp(fp).xor_sign(s_x)
    }

    pub(crate) fn from_int_diff<I: AnyInt>(lhs: I, rhs: I) -> Self
    {
        let s = lhs < rhs;
        let abs_diff = if s
        {
            rhs - lhs
        }
        else
        {
            lhs - rhs
        };
        Self::from_int(abs_diff).xor_sign(s)
    }
    
    /// Converts an integer into a custom floating-point type.
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn from_int<I: AnyInt>(mut from: I) -> Self
    {
        let s = from < I::zero();

        macro_rules! sat {
            ($expr:expr) => {
                match $expr
                {
                    Ok(y) => y,
                    Err(Saturation::Overflow) => return Fp::infinity().extra_sign(s),
                    Err(Saturation::Underflow) => return Fp::zero().extra_sign(s)
                }
            };
        }

        let overflow = s && from == I::min_value();
        from = if overflow
        {
            I::max_value()
        }
        else
        {
            util::abs(from)
        };
        let mut e = Fp::exp_bias();
        let f = if util::bitsize_of::<I>() - 1 > util::bitsize_of::<U>()
        {
            let mut f = from;
            if overflow && let Some(base) = I::from(EXP_BASE)
            {
                f = f/base + I::one();
                e = e + U::one();
            }
            sat!(Fp::mantissa_shl(&mut f, FRAC_SIZE, &mut e, None));
            Fp::normalize_mantissa(&mut e, &mut f, None);
            sat!(Fp::convert_mantissa(f, &mut e, None))
        }
        else
        {
            let mut f = if overflow
            {
                U::one() << (util::bitsize_of::<I>() - 1)
            }
            else
            {
                <U as NumCast>::from(from).unwrap()
            };
            sat!(Fp::mantissa_shl(&mut f, FRAC_SIZE, &mut e, None));
            Fp::normalize_mantissa(&mut e, &mut f, None);
            f
        };

        Fp::from_exp_mantissa(e, f).extra_sign(s)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> for Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn from(value: Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>) -> Self
    {
        value.extra_sign(false)
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn from(value: Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>) -> Self
    {
        value.collect()
    }
}

macro_rules! impl_from_int {
    ($($i:ty),*) => {
        $(
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$i> for Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                fn from(n: $i) -> Self
                {
                    Self::from_int(n)
                }
            }
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$i> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                fn from(n: $i) -> Self
                {
                    Self::from_int(n)
                }
            }
        )*
    };
}

impl_from_int!(
    u8, u16, u32, usize, u64, u128,
    i8, i16, i32, isize, i64, i128
);

/*#[cfg(feature = "ethnum")]
impl_from_int!(ethnum::U256, ethnum::I256);*/

macro_rules! impl_from_float {
    ($($f:ty: $b:ty, $fp:ty),*) => {
        $(
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$f> for Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                fn from(f: $f) -> Self
                {
                    Self::from_fp::<$b, _, _, _, _, 2>(<$fp>::from_bits(f.to_bits()))
                }
            }
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<$f> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                fn from(f: $f) -> Self
                {
                    Self::from_fp::<$b, _, _, _, _, 2>(<$fp>::from_bits(f.to_bits()))
                }
            }
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> for $f
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                fn from(value: Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>) -> Self
                {
                    <$f>::from_bits(<$fp>::from_fp(value).to_bits())
                }
            }
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> From<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> for $f
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                fn from(value: Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>) -> Self
                {
                    <$f>::from_bits(<$fp>::from_fp(value).to_bits())
                }
            }
        )*
    };
}
impl_from_float!(f16: u16, FpHalf, f32: u32, FpSingle, f64: u64, FpDouble, f128: u128, FpQuadruple);