use core::cmp::Ordering;
use core::num::FpCategory;
use core::ops::Neg;

use num_traits::{Float, Signed, ConstZero, FloatConst, FromBytes, NumCast, ToBytes, Zero};

use crate::util::Saturation;
use crate::{util, AnyInt, FpRepr, Int, UInt};

moddef::moddef!(
    flat(pub) mod {
        bytemuck for cfg(feature = "bytemuck"),
        cast,
        cmp,
        fmt,
        num,
        ops,
    },
    flat(pub) mod {
        default,
        parse
    }
);

const NO_NEWTON: bool = false;
const NEWTON_EXP: usize = if NO_NEWTON {0} else {3};
const NEWTON_LN: usize = if NO_NEWTON {0} else {2};
const NEWTON_RT: usize = if NO_NEWTON {0} else {4};
const NEWTON_TRIG: usize = if NO_NEWTON {0} else {3};

macro_rules! as_lossless {
    ($value:expr, $fn_as_lossless:expr, $fn:block) => {
        {
            #[cfg(any(
                feature = "use_std_float",
                all(debug_assertions, test)
            ))]
            let _as_lossless = crate::Fp::as_lossless(
                $value,
                $fn_as_lossless,
                $fn_as_lossless,
                $fn_as_lossless,
                $fn_as_lossless
            );
            #[cfg(all(not(test), feature = "use_std_float"))]
            if let Some([as_lossless]) = _as_lossless
            {
                return as_lossless
            }
            #[allow(clippy::redundant_closure_call)]
            let y = (|| $fn)();
            #[cfg(all(debug_assertions, test))]
            if let Some([as_lossless]) = _as_lossless
            {
                if !y.approx_eq(as_lossless)
                {
                    debug_assert_eq!(y, as_lossless, "Error is too big!")
                }
            }
            y
        }
    };
    ($value:expr, $fn_as_lossless:expr, $fn_as_lossless_alt:expr, $fn:block) => {
        {
            #[cfg(any(
                feature = "use_std_float",
                all(debug_assertions, test)
            ))]
            let _as_lossless = crate::Fp::as_lossless(
                $value,
                $fn_as_lossless,
                $fn_as_lossless,
                $fn_as_lossless,
                $fn_as_lossless
            );
            #[cfg(all(not(test), feature = "use_std_float"))]
            if let Some([as_lossless]) = _as_lossless
            {
                return as_lossless
            }
            #[allow(clippy::redundant_closure_call)]
            let y = (|| $fn)();
            #[cfg(all(debug_assertions, test))]
            if let Some([as_lossless]) = _as_lossless
            {
                if !y.approx_eq(as_lossless) && crate::Fp::as_lossless(
                    $value,
                    $fn_as_lossless_alt,
                    $fn_as_lossless_alt,
                    $fn_as_lossless_alt,
                    $fn_as_lossless_alt
                ).is_none_or(|[alt]| {
                    !y.approx_eq(alt) && y > as_lossless.max(alt) && y < as_lossless.min(alt)
                })
                {
                    debug_assert_eq!(y, as_lossless, "Error is too big!")
                }
            }
            y
        }
    };
}
use as_lossless as as_lossless;

#[derive(Clone, Copy)]
#[cfg_attr(feature = "bytemuck", derive(::bytemuck::Zeroable))]
struct Fps<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>(
    pub Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>,
    pub bool
)
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>;

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    fn collect(self) -> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
    {
        let Self (y, s) = self;
        y.xor_sign(s)
    }
}

/// A custom floating point type, where the bit size of the exponent and mantissa can be set separately.
/// 
/// `U` is the underlying unsigned integer type which is used to represent the number.
///
/// `SIGN_BIT` is wether or not the number has a sign bit.
/// 
/// `EXP_SIZE` is the size of the exponent in bits.
/// 
/// `INT_SIZE` is the size of the integer part of the mantissa in bits. If zero, then the integer bit is implicit.
/// 
/// `FRAC_SIZE` is the size of the fractional part of the mantissa in bits.
/// 
/// `EXP_BASE` is the base of the exponent.
/// 
/// The total bit size of `U` must be greater or equal to `SIGN_BIT` + `EXP_SIZE` + `INT_SIZE` + `FRAC_SIZE` to contain the entire number.
///
/// The bit layout is as follows:
/// ```txt
/// No data: | Sign:      | Exponent:  | Integer:   | Fractional: |
/// <  ..  > | <SIGN_BIT> | <EXP_SIZE> | <INT_SIZE> | <FRAC_SIZE> |
/// ```
/// 
/// The value of a real floating-point number is the following:
/// ```txt
/// x = (-1)**sign*EXP_BASE**(exponent - bias)*mantissa
/// ```
/// 
/// where the bias equals
/// ```txt
/// bias = 2**(EXP_SIZE - 1) - 1
/// ```
///
/// If the exponent has the maximum value, the number is either infinity or NaN.
/// 
/// The number then automatically implements `num::Float`, and supports all ordinary floating point operations.
/// 
/// This allows simple implementation of special floating point types, such as TensorFloat, IEEE754 Quadruple/binary128, Fp80, and BFloat16.
/// 
/// The accuracy of all of the floating point operations are not perfect, but work well enough to be usable. Various plots showing the accuracy of basic functions are shown in the [plots](https://github.com/sigurd4/custom_float/tree/master/plots) subfolder.
/// 
/// All floats can be converted into each other painlessly, though the conversion may produce rounding errors or unbounded outputs when converting to a float with lesser resolution.
///
/// # Examples
///
/// ```rust
/// #![feature(generic_const_exprs)]
/// 
/// use custom_float::Fp;
///
/// type FpSingle = Fp<u32, true, 8, 0, 23, 2>;
///
/// let two = FpSingle::from(2);
/// let four = FpSingle::from(4);
/// 
/// assert_eq!(two + two, four);
/// ```
#[derive(Clone, Copy)]
#[cfg_attr(feature = "bytemuck", derive(::bytemuck::Pod, ::bytemuck::Zeroable))]
#[repr(transparent)]
pub struct Fp<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize>(U)
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>;

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Size of floating-point number in bits
    pub const BIT_SIZE: usize = Self::SIGN_SIZE + Self::EXP_SIZE + Self::INT_SIZE + Self::FRAC_SIZE;
    /// Size of the sign bit
    pub const SIGN_SIZE: usize = Self::SIGN_BIT as usize;

    /// Wether or not the number has a sign bit. If not, it can only be positive.
    pub const SIGN_BIT: bool = SIGN_BIT;
    /// The size of the exponent part in bits
    pub const EXP_SIZE: usize = EXP_SIZE;
    /// The size of the integer part in bits
    pub const INT_SIZE: usize = INT_SIZE;
    /// The size of the fractional part in bits
    pub const FRAC_SIZE: usize = FRAC_SIZE;
    /// The base for the exponent
    pub const EXP_BASE: usize = EXP_BASE;

    /// Position of the sign bit
    pub const SIGN_POS: usize = Self::EXP_SIZE + Self::INT_SIZE + Self::FRAC_SIZE;
    /// Position of the first exponent bit
    pub const EXP_POS: usize = Self::INT_SIZE + Self::FRAC_SIZE;
    /// Position of the first integer bit
    pub const INT_POS: usize = Self::FRAC_SIZE;
    /// Position of the first fractional bit
    pub const FRAC_POS: usize = 0;

    /// Number of significant digits in base 2.
    pub const MANTISSA_DIGITS: usize = Self::INT_SIZE + Self::FRAC_SIZE;

    /// `true` if the number contains an implicit integer bit
    pub const IS_INT_IMPLICIT: bool = Self::INT_SIZE == 0;

    const BASE_FACTORS_2: &[[usize; 2]] = util::factorize(EXP_BASE);
    const BASE_FACTORS_3: &[[usize; 3]] = util::factorize(EXP_BASE);

    const MANTISSA_OP_SIZE: usize = FRAC_SIZE + INT_SIZE + Self::IS_INT_IMPLICIT as usize;
    const BASE_PADDING: usize = util::base_padding(EXP_BASE);
    const BASE_FACTORS_2_PADDING: &[[usize; 2]] = util::base_factor_paddings(Self::BASE_FACTORS_2);
    //const BASE_FACTORS_3_PADDING: &[[usize; 3]] = util::base_factor_paddings(Self::BASE_FACTORS_3);
    const HALF_BASE_PADDING: usize = util::base_padding(EXP_BASE/2);

    #[cfg(any(
        feature = "use_std_float",
        all(debug_assertions, test)
    ))]
    fn as_lossless<const N: usize, const M: usize>(
        value: [Self; N],
        as_f16: impl FnOnce([f16; N]) -> [f16; M],
        as_f32: impl FnOnce([f32; N]) -> [f32; M],
        as_f64: impl FnOnce([f64; N]) -> [f64; M],
        as_f128: impl FnOnce([f128; N]) -> [f128; M]
    ) -> Option<[Self; M]>
    {
        macro_rules! as_float {
            ($($t:ty => $fn:expr),*) => {
                $(
                    if util::is_float_conversion_mutually_lossless::<Self, $t>()
                    {
                        return Some($fn(value.map(<$t as From<Self>>::from)).map(Self::from))
                    }
                )*
            };
        }

        as_float!(
            f16 => as_f16,
            f32 => as_f32,
            f64 => as_f64,
            f128 => as_f128
        );

        None
    }

    const fn extra_sign(self, s: bool) -> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
    {
        Fps(self, s)
    }

    pub fn sign_eq(self, rhs: Self) -> bool
    {
        if !SIGN_BIT
        {
            return true
        }
        let mask = Self::sign_mask();
        (self.to_bits() & mask) == (rhs.to_bits() & mask)
    }

    pub fn approx_eq(mut self, mut rhs: Self) -> bool
    {
        let (c1, c2) = (self.classify(), rhs.classify());
        let sign_eq = self.sign_eq(rhs) || matches!(c1, FpCategory::Zero) || matches!(c2, FpCategory::Zero);
        match (c1, c2)
        {
            (FpCategory::Nan, FpCategory::Nan)
                | (FpCategory::Zero | FpCategory::Subnormal, FpCategory::Zero | FpCategory::Subnormal) => true,
            (FpCategory::Infinite, FpCategory::Infinite) => sign_eq,
            (FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero, FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero) => {
                if matches!(c1, FpCategory::Infinite)
                {
                    self = Self::max_value().copysign(self)
                }
                if matches!(c2, FpCategory::Infinite)
                {
                    rhs = Self::max_value().copysign(rhs)
                }
                if !Self::IS_INT_IMPLICIT
                {
                    self.normalize_up();
                    rhs.normalize_up();
                }
                
                if let Some(e0) = self.exp_bits().to_usize() && let Some(f0) = <u128 as NumCast>::from(self.mantissa_bits())
                    && let Some(e1) = rhs.exp_bits().to_usize() && let Some(f1) = <u128 as NumCast>::from(rhs.mantissa_bits())
                {
                    const PRECISION: usize = 2;
                    let div = Self::FRAC_SIZE.saturating_sub(PRECISION + 1) + 1;
        
                    let quantize = |mut e: usize, f: u128| if f != 0
                    {
                        if EXP_BASE != 2
                        {
                            e = util::pow_ilog2(e, EXP_BASE);
                        }
                        e += f.ilog2() as usize;
                        e/div
                    }
                    else
                    {
                        0
                    };
        
                    let a = quantize(e0, f0);
                    let b = quantize(e1, f1);

                    if if sign_eq
                    {
                        a.max(b) - a.min(b) <= 1
                    }
                    else
                    {
                        a + b <= 1
                    }
                    {
                        return true
                    }
                }

                false
            },
            _ => false
        }
    }

    #[inline]
    fn sign_mask() -> U
    {
        Self::shift_sign(U::one())
    }
    #[inline]
    fn max_exponent_bits() -> U
    {
        U::max_value() >> (util::bitsize_of::<U>() - EXP_SIZE)
    }
    #[inline]
    fn exponent_mask() -> U
    {
        Self::shift_exp(Self::max_exponent_bits())
    }
    #[inline]
    fn max_int_bits() -> U
    {
        U::max_value() >> (util::bitsize_of::<U>() - INT_SIZE)
    }
    #[inline]
    fn int_mask() -> U
    {
        Self::shift_int(Self::max_int_bits())
    }
    #[inline]
    fn max_frac_bits() -> U
    {
        U::max_value() >> (util::bitsize_of::<U>() - FRAC_SIZE)
    }
    #[inline]
    fn frac_mask() -> U
    {
        Self::shift_frac(Self::max_frac_bits())
    }
    #[inline]
    fn max_mantissa_bits() -> U
    {
        U::max_value() >> (util::bitsize_of::<U>() - Self::MANTISSA_DIGITS - Self::IS_INT_IMPLICIT as usize)
    }
    /*#[inline]
    fn mantissa_mask() -> U
    {
        Self::shift_frac(Self::max_mantissa_bits())
    }*/
    #[inline]
    fn max_mantissa_bits_raw() -> U
    {
        U::max_value() >> (util::bitsize_of::<U>() - Self::MANTISSA_DIGITS)
    }
    #[inline]
    fn mantissa_mask_raw() -> U
    {
        Self::shift_frac(Self::max_mantissa_bits_raw())
    }

    fn shift_sign(s: U) -> U
    {
        match SIGN_BIT
        {
            true => s << Self::SIGN_POS,
            false => Zero::zero()
        }
    }
    fn shift_exp(e: U) -> U
    {
        match EXP_SIZE != 0
        {
            true => e << Self::EXP_POS,
            false => Zero::zero()
        }
    }
    fn shift_int(e: U) -> U
    {
        match INT_SIZE != 0
        {
            true => e << Self::INT_POS,
            false => Zero::zero()
        }
    }
    fn shift_frac(f: U) -> U
    {
        match FRAC_SIZE != 0
        {
            true => f << Self::FRAC_POS,
            false => Zero::zero()
        }
    }
    fn shift_mantissa(f: U) -> U
    {
        match FRAC_SIZE != 0 || INT_SIZE != 0
        {
            true => f << Self::FRAC_POS,
            false => Zero::zero()
        }
    }

    fn mantissa_bits(&self) -> U
    {
        let mut f = self.mantissa_bits_raw();
        if Self::IS_INT_IMPLICIT
        {
            Self::mantissa_from_raw(&mut f, self.is_subnormal());
        }
        f
    }
    fn mantissa_bits_raw(&self) -> U
    {
        (self.to_bits() & Self::mantissa_mask_raw()) >> Self::FRAC_POS
    }

    fn mantissa_to_raw(mantissa: &mut U, subnormal: bool)
    {
        if Self::IS_INT_IMPLICIT
        {
            let implicit_one = U::one() << Self::INT_POS;
            if subnormal
            {
                *mantissa = util::rounding_div_2(*mantissa)
            }
            else
            {
                debug_assert!(*mantissa >= implicit_one, "Mantissa is less than one even though it's supposed to be normal.");
                *mantissa = *mantissa - implicit_one;
            }
            debug_assert!(*mantissa < implicit_one, "Mantissa is impossibly large!")
        }
    }
    fn mantissa_from_raw(mantissa_raw: &mut U, subnormal: bool)
    {
        if Self::IS_INT_IMPLICIT
        {
            if subnormal
            {
                *mantissa_raw = *mantissa_raw << 1usize
            }
            else
            {
                *mantissa_raw = *mantissa_raw | (U::one() << Self::INT_POS)
            }
        }
    }

    fn from_exp_mantissa(exp: U, mut mantissa: U) -> Self
    {
        if mantissa.is_zero()
        {
            return Self::zero()
        }
        if exp >= Self::max_exponent_bits()
        {
            return Self::infinity()
        }
        let subnormal = Self::IS_INT_IMPLICIT && exp.is_zero();
        Self::mantissa_to_raw(&mut mantissa, subnormal);
        Self::from_bits(Self::shift_exp(exp) | Self::shift_mantissa(mantissa))
    }

    fn add_nan(self, rhs: Self) -> Self
    {
        if rhs.is_snan() || !self.is_nan()
        {
            rhs
        }
        else
        {
            self
        }
    }

    fn convert_mantissa<I: AnyInt>(mut f: I, e: &mut U, mut e_minus: Option<&mut U>) -> Result<U, Saturation>
    {
        loop
        {
            if (
                f >= I::zero()
                || (f != I::min_value() && {
                    f = util::neg(f);
                    true
                })
            ) && let Some(ff) = NumCast::from(f)
            {
                break Ok(ff)
            }
            if let Some(base) = I::from(EXP_BASE)
            {
                if !EXP_BASE.is_zero() && util::complementary_add_sub_assign(Some(e), e_minus.as_deref_mut(), U::one()).is_ok()
                {
                    f = util::rounding_div(f, base);
                    continue
                }
                return Err(Saturation::Overflow)
            }
            else
            {
                return Err(Saturation::Underflow)
            }
        }
    }

    fn normalize_up(&mut self)
    {
        let s = self.is_sign_negative();
        let mut e = self.exp_bits();
        let mut m = self.mantissa_bits();
        Self::normalize_mantissa_up(&mut e, &mut m, None);
        *self = Self::from_exp_mantissa(e, m).with_sign(s)
    }
    /*fn normalize_down(&mut self)
    {
        let s = self.is_sign_negative();
        let mut e = self.exp_bits();
        let mut m = self.mantissa_bits();
        Self::normalize_mantissa_down(&mut e, &mut m, None);
        *self = Self::from_exp_mantissa(e, m).with_sign(s)
    }
    fn normalize(&mut self)
    {
        let s = self.is_sign_negative();
        let mut e = self.exp_bits();
        let mut m = self.mantissa_bits();
        Self::normalize_mantissa(&mut e, &mut m, None);
        *self = Self::from_exp_mantissa(e, m).with_sign(s)
    }*/

    fn normalize_mantissa_down<M: AnyInt>(exp: &mut U, mantissa: &mut M, target_exp: Option<U>)
    {
        if exp.is_zero()
        {
            return
        }
        if mantissa.is_zero()
        {
            *exp = U::zero();
            return
        }
        let min_exp = target_exp.unwrap_or_else(U::zero);
        if EXP_BASE.is_power_of_two()
            && util::bitsize_of::<M>() + Self::BASE_PADDING > Self::MANTISSA_OP_SIZE
            && let Some(min_lz) = <u32 as NumCast>::from(util::bitsize_of::<M>() + Self::BASE_PADDING - Self::MANTISSA_OP_SIZE - 1)
            && let Some(shifts) = <usize as NumCast>::from(exp.saturating_sub(min_exp))
        {
            *exp = *exp - Self::shl_mantissa_without_loss(mantissa, Some(shifts), EXP_BASE.ilog2(), Some(min_lz));
            return
        }
        let base = M::from(EXP_BASE);
        if let Some(base) = base
        {
            while *exp > min_exp && (*mantissa >> (Self::MANTISSA_OP_SIZE - Self::BASE_PADDING)).is_zero()
            {
                *exp = *exp - U::one();
                *mantissa = *mantissa*base;
            }
        }
    }
    fn normalize_mantissa_up<M: AnyInt>(exp: &mut U, mantissa: &mut M, target_exp: Option<U>)
    {
        if mantissa.is_zero()
        {
            *exp = U::zero();
            return
        }
        let max_exp = target_exp.unwrap_or_else(|| Self::max_exponent_bits());
        if EXP_BASE.is_power_of_two()
            && util::bitsize_of::<M>() > Self::MANTISSA_OP_SIZE
            && let Some(max_lz) = <u32 as NumCast>::from(util::bitsize_of::<M>() - Self::MANTISSA_OP_SIZE - 1)
            && let Some(shifts) = <usize as NumCast>::from(max_exp.saturating_sub(*exp))
        {
            *exp = *exp + Self::shr_mantissa_without_loss(mantissa, Some(shifts), EXP_BASE.ilog2(), Some(max_lz));
        }
        if let Some(base) = M::from(EXP_BASE)
        {
            while *exp < max_exp && !(*mantissa >> Self::MANTISSA_OP_SIZE).is_zero()
            {
                if EXP_BASE == 0
                {
                    *exp = Self::max_exponent_bits();
                    *mantissa = M::zero();
                    return
                }
                *exp = *exp + U::one();
                *mantissa = util::rounding_div(*mantissa, base)
            }
        }
        else if *exp < max_exp && !(*mantissa >> Self::MANTISSA_OP_SIZE).is_zero()
        {
            *exp = U::zero();
            *mantissa = M::zero()
        }
    }

    fn normalize_mantissa<M: AnyInt>(exp: &mut U, mantissa: &mut M, target_exp: Option<U>)
    {
        Self::normalize_mantissa_down(exp, mantissa, target_exp);
        Self::normalize_mantissa_up(exp, mantissa, target_exp);
    }

    fn mantissa_shl<M: AnyInt>(mantissa: &mut M, mut shifts: usize, exp: &mut U, mut exp_offset: Option<&mut U>) -> Result<(), Saturation>
    {
        if mantissa.is_zero()
        {
            return Err(Saturation::Underflow)
        }
        if shifts == 0
        {
            return Ok(())
        }
        shifts -= Self::shl_mantissa_without_loss::<_, usize>(mantissa, Some(shifts), 1, None);
        if EXP_BASE.is_power_of_two() && let Some(mut change) = U::from(shifts/EXP_BASE.ilog2() as usize)
        {
            if let Some(ofs) = exp_offset
            {
                if let Some(diff) = ofs.checked_sub(&change)
                {
                    *ofs = diff;
                    return Ok(())
                }
                change = change - *ofs;
                *ofs = U::zero();
            }
            if let Some(diff) = exp.checked_add(&change)
            {
                *exp = diff;
                return Ok(())
            }
            Err(Saturation::Overflow)
        }
        else if EXP_BASE.is_multiple_of(2) && let Some(half_base) = M::from(EXP_BASE/2)
        {
            let mut step = 1usize;
            loop
            {
                if !EXP_BASE.is_multiple_of(2 << step)
                {
                    break
                }
                step += 1;
            }
            loop
            {
                if mantissa.is_zero()
                {
                    return Err(Saturation::Underflow)
                }

                let lz = mantissa.leading_zeros();
                if lz > 0
                {
                    let s = shifts.min(lz as usize);
                    shifts -= s;
                    *mantissa = *mantissa << s;
                    if shifts == 0
                    {
                        break
                    }
                }

                if let Some(ofs) = &mut exp_offset && **ofs > U::zero()
                {
                    **ofs = **ofs - U::one();
                }
                else if *exp < U::max_value()
                {
                    *exp = *exp + U::one();
                }
                else
                {
                    return Err(Saturation::Overflow)
                }
                let s = shifts.min(step);
                shifts -= s;
                let div = half_base >> (s - 1);
                if div.is_zero()
                {
                    return Err(Saturation::Overflow)
                }
                *mantissa = util::rounding_div(*mantissa, div);
                if shifts == 0
                {
                    break
                }
            }
            Ok(())
        }
        else if let Some(base) = M::from(EXP_BASE)
        {
            loop
            {
                if mantissa.is_zero()
                {
                    return Err(Saturation::Underflow)
                }
                
                let lz = mantissa.leading_zeros();
                if lz > 0
                {
                    let s = shifts.min(lz as usize);
                    shifts -= s;
                    *mantissa = *mantissa << s;
                    if shifts == 0
                    {
                        break
                    }
                }

                if let Some(ofs) = &mut exp_offset && **ofs > U::zero()
                {
                    **ofs = **ofs - U::one();
                }
                else if *exp < U::max_value()
                {
                    *exp = *exp + U::one();
                }
                else
                {
                    return Err(Saturation::Overflow)
                }
                if EXP_BASE.is_zero()
                {
                    return Err(Saturation::Overflow)
                }
                *mantissa = util::rounding_div(*mantissa, base);
            }
            Ok(())
        }
        else
        {
            Err(Saturation::Overflow)
        }
    }

    fn mantissa_shr<M: AnyInt>(mantissa: &mut M, mut shifts: usize, exp: &mut U) -> Result<(), Saturation>
    {
        shifts -= Self::shr_mantissa_without_loss::<_, usize>(mantissa, Some(shifts), 1, None);
        if shifts == 0
        {
            return Ok(())
        }
        if EXP_BASE.is_power_of_two() && let Some(mut change) = U::from(shifts/EXP_BASE.ilog2() as usize)
        {
            if let Some(diff) = exp.checked_sub(&change)
            {
                *exp = diff;
                return Ok(())
            }
            change = change - *exp;
            *exp = U::zero();
            if let Some(unshift) = change.to_usize()
            {
                shifts -= unshift*EXP_BASE.ilog2() as usize
            }
            else
            {
                return Err(Saturation::Overflow)
            }
            *mantissa = util::rounding_div_2_pow(*mantissa, shifts);
            return Ok(())
        }
        if EXP_BASE.is_multiple_of(2) && let Some(half_base) = M::from(EXP_BASE / 2)
        {
            while !shifts.is_zero()
            {
                if mantissa.leading_zeros() as usize > Self::HALF_BASE_PADDING && *exp > U::zero()
                {
                    *exp = *exp - U::one();
                    *mantissa = *mantissa*half_base;
                }
                else
                {
                    *mantissa = util::rounding_div_2(*mantissa);
                }
                shifts -= 1;
            }
            return Ok(())
        }
        if let Some(base) = M::from(EXP_BASE)
        {
            while !shifts.is_zero()
            {
                while mantissa.leading_zeros() as usize > Self::BASE_PADDING
                {
                    if *exp > U::zero()
                    {
                        *exp = *exp - U::one();
                    }
                    else
                    {
                        break
                    }
                    *mantissa = *mantissa*base;
                }
                *mantissa = util::rounding_div_2(*mantissa);
                shifts -= 1;
            }
            return Ok(())
        }
        Err(Saturation::Underflow)
    }

    fn force_exp_to<M: AnyInt>(exp: &mut U, mantissa: &mut M, target_exp: U) -> Result<(), Saturation>
    {
        if mantissa.is_zero()
        {
            return Err(Saturation::Underflow)
        }
        if *exp > target_exp
        {
            let base = match M::from(EXP_BASE)
            {
                Some(b) => b,
                None => return Err(Saturation::Overflow)
            };
            loop
            {
                *mantissa = match mantissa.checked_mul(&base)
                {
                    Some(m) => m,
                    None => return Err(Saturation::Overflow)
                };
                *exp = *exp - U::one();
                if *exp <= target_exp
                {
                    break
                }
            }
        }
        else if *exp < target_exp
        {
            let base = match M::from(EXP_BASE)
            {
                Some(b) => b,
                None => return Err(Saturation::Underflow)
            };
            if EXP_BASE.is_zero()
            {
                return Err(Saturation::Overflow)
            }
            loop
            {
                *mantissa = util::rounding_div(*mantissa, base);
                if mantissa.is_zero()
                {
                    return Err(Saturation::Underflow)
                }
                *exp = *exp + U::one();
                if *exp >= target_exp
                {
                    break
                }
            }
        }
        assert_eq!(*exp, target_exp);
        Ok(())
    }

    fn shr_mantissa_without_loss<M: AnyInt, I: UInt>(mantissa: &mut M, shifts: Option<usize>, quanta: u32, max_lz: Option<u32>) -> I
    {
        if !shifts.is_some_and(|s| s.is_zero()) && !mantissa.is_zero()
        {
            let mut i: u32 = mantissa.trailing_zeros();
            if let Some(mlz) = max_lz
            {
                let lz = mantissa.leading_zeros();
                i = i.min(mlz.saturating_sub(lz))
            }
            i -= i % quanta;
            if let Some(s) = shifts.and_then(<u32 as NumCast>::from)
            {
                i = i.min(s*quanta)
            }
            if i != 0 && let Some(o) = I::from(i/quanta)
            {
                *mantissa = *mantissa >> i;
                return o
            }
        }
        I::zero()
    }

    fn shl_mantissa_without_loss<M: AnyInt, I: UInt>(mantissa: &mut M, shifts: Option<usize>, quanta: u32, min_lz: Option<u32>) -> I
    {
        if !mantissa.is_zero()
        {
            let mut i: u32 = mantissa.leading_zeros().saturating_sub(min_lz.unwrap_or(util::is_signed::<M>() as u32));
            i = i - i % quanta;
            if let Some(s) = shifts.and_then(<u32 as NumCast>::from)
            {
                i = i.min(s*quanta)
            }
            if i != 0 && let Some(o) = I::from(i/quanta)
            {
                *mantissa = *mantissa << i;
                return o
            }
        }
        I::zero()
    }

    fn mantissa_from_low_high(mut mantissa: U, overflow: U, exp: &mut U) -> U
    {
        assert!(EXP_BASE.is_power_of_two());
        if overflow.is_zero()
        {
            return mantissa
        }
        let quanta = EXP_BASE.ilog2() as usize;
        let mut lz = overflow.leading_zeros() as usize;
        let mut shift = util::bitsize_of::<U>() - lz;
        let extra = (quanta - shift % quanta) % quanta;
        shift += extra;
        lz -= extra;
        if let Some(offs) = U::from(shift/quanta)
        {
            *exp = *exp + offs;
            mantissa = util::rounding_div_2_pow(mantissa, shift);
            return mantissa | (overflow << lz)
        }
        *exp = U::zero();
        U::zero()
    }

    fn mantissa_div_sqrt_base_or_base(mantissa: &mut U) -> Result<U, ()>
    {
        assert!(!EXP_BASE.is_zero());
        if let Some(base) = U::from(EXP_BASE) && (*mantissa % base).is_zero()
        {
            *mantissa = *mantissa/base;
            Ok(U::one() << 1usize)
        }
        else if EXP_BASE.isqrt()*EXP_BASE.isqrt() == EXP_BASE && let Some(base_sqrt) = U::from(EXP_BASE.isqrt())
        {
            *mantissa = util::rounding_div(*mantissa, base_sqrt);
            Ok(U::one())
        }
        else if let Some(base) = U::from(EXP_BASE)
        {
            *mantissa = util::rounding_div(*mantissa, base);
            Ok(U::one() << 1usize)
        }
        else
        {
            Err(())
        }
    }

    fn mantissa_pairs_div_either_base(mantissa1: &mut U, mantissa2: &mut U) -> Result<(), ()>
    {
        assert!(!EXP_BASE.is_zero());
        fn sort<'a, U: UInt>(mantissa1: &'a mut U, mantissa2: &'a mut U) -> [&'a mut U; 2]
        {
            if mantissa1 < mantissa2
            {
                [mantissa2, mantissa1]
            }
            else
            {
                [mantissa1, mantissa2]
            }
        }

        let mut mantissas = sort(mantissa1, mantissa2);
        let [&mut m_most, &mut m_least] = mantissas;

        if let Some(base) = U::from(EXP_BASE)
        {
            for m in mantissas.iter_mut()
            {
                if (**m % base).is_zero()
                {
                    **m = **m/base;

                    return Ok(())
                }
            }
        }

        if !Self::BASE_FACTORS_2.is_empty()
        {
            for &[base1, base2] in Self::BASE_FACTORS_2
            {
                if let Some(b1) = U::from(base1) && let Some(b2) = U::from(base2)
                    && (m_most % b1).is_zero() && (m_least % b2).is_zero()
                {
                    for (m, b) in mantissas.iter_mut()
                        .zip([b1, b2])
                    {
                        **m = **m/b
                    }
                    
                    return Ok(())
                }
            }
        }

        if !Self::BASE_FACTORS_2.is_empty()
        {
            let (_, [b1, b2]) = Self::BASE_FACTORS_2.iter()
                .filter_map(|&[base1, base2]| {
                    if let Some(b1) = U::from(base1) && let Some(b2) = U::from(base2)
                        && let k = (m_most % b1).is_zero() as u8 + (m_least % b2).is_zero() as u8
                        && k > 0
                    {
                        Some((k, [b1, b2]))
                    }
                    else
                    {
                        None
                    }
                }).reduce(|a, b| {
                    if a.0 < b.0
                    {
                        b
                    }
                    else
                    {
                        a
                    }
                }).unwrap();

            for (m, b) in mantissas.iter_mut()
                .zip([b1, b2])
            {
                **m = util::rounding_div(**m, b)
            }
            return Ok(())
        }

        if let Some(base) = U::from(EXP_BASE)
        {
            let [m_most, _] = &mut mantissas;
            **m_most = util::rounding_div(**m_most, base);
            return Ok(())
        }

        Err(())
    }

    fn mantissa_pairs_div_either_base_mul(mantissa1: &mut U, mantissa2: &mut U) -> Result<Option<U>, ()>
    {
        fn sort<'a, U: UInt>(mantissa1: &'a mut U, mantissa2: &'a mut U) -> [&'a mut U; 2]
        {
            if mantissa1 < mantissa2
            {
                [mantissa2, mantissa1]
            }
            else
            {
                [mantissa1, mantissa2]
            }
        }

        let mut mantissas = sort(mantissa1, mantissa2);
        let [&mut m_most, &mut m_least] = mantissas;

        if let Some(base) = U::from(EXP_BASE)
        {
            for m in mantissas.iter_mut()
            {
                if (**m % base).is_zero()
                {
                    **m = **m/base;

                    return Ok(mantissa1.checked_mul(mantissa2))
                }
            }
        }

        if !Self::BASE_FACTORS_2.is_empty()
        {
            for &[base1, base2] in Self::BASE_FACTORS_2
            {
                if let Some(b1) = U::from(base1) && let Some(b2) = U::from(base2)
                    && (m_most % b1).is_zero() && (m_least % b2).is_zero()
                {
                    for (m, b) in mantissas.iter_mut()
                        .zip([b1, b2])
                    {
                        **m = **m/b
                    }
                    
                    return Ok(mantissa1.checked_mul(mantissa2))
                }
            }
        }

        if !Self::BASE_FACTORS_3.is_empty()
        {
            for &[base_final, base1, base2] in Self::BASE_FACTORS_3
            {
                if let Some(bf) = U::from(base_final)
                    && let Some(b1) = U::from(base1) && let Some(b2) = U::from(base2)
                    && (m_most % b1).is_zero() && (m_least % b2).is_zero()
                {
                    for (m, b) in mantissas.iter_mut()
                        .zip([b1, b2])
                    {
                        **m = **m/b
                    }
                    return Ok(mantissa1.checked_mul(mantissa2).map(|y| util::rounding_div(y, bf)))
                }
            }
            let (_, [bf, b1, b2]) = Self::BASE_FACTORS_3.iter()
                .filter_map(|&[base_final, base1, base2]| {
                    if let Some(bf) = U::from(base_final)
                        && let Some(b1) = U::from(base1) && let Some(b2) = U::from(base2)
                        && let k = (m_most % b1).is_zero() as u8 + (m_least % b2).is_zero() as u8
                        && k > 0
                    {
                        Some((k, [bf, b1, b2]))
                    }
                    else
                    {
                        None
                    }
                }).reduce(|a, b| {
                    if a.0 < b.0
                    {
                        b
                    }
                    else
                    {
                        a
                    }
                }).unwrap();

            for (m, b) in mantissas.iter_mut()
                .zip([b1, b2])
            {
                **m = util::rounding_div(**m, b)
            }
            return Ok(mantissa1.checked_mul(mantissa2).map(|y| util::rounding_div(y, bf)))
        }

        if !Self::BASE_FACTORS_2.is_empty()
        {
            let (_, [b1, b2]) = Self::BASE_FACTORS_2.iter()
                .filter_map(|&[base1, base2]| {
                    if let Some(b1) = U::from(base1) && let Some(b2) = U::from(base2)
                        && let k = (m_most % b1).is_zero() as u8 + (m_least % b2).is_zero() as u8
                        && k > 0
                    {
                        Some((k, [b1, b2]))
                    }
                    else
                    {
                        None
                    }
                }).reduce(|a, b| {
                    if a.0 < b.0
                    {
                        b
                    }
                    else
                    {
                        a
                    }
                }).unwrap();

            for (m, b) in mantissas.iter_mut()
                .zip([b1, b2])
            {
                **m = util::rounding_div(**m, b)
            }
            return Ok(mantissa1.checked_mul(mantissa2))
        }

        if let Some(base) = U::from(EXP_BASE)
        {
            let [m_most, _] = &mut mantissas;
            **m_most = util::rounding_div(**m_most, base);
            return Ok(mantissa1.checked_mul(mantissa2))
        }

        Err(())
    }

    fn mantissa_pairs_div_both_base_add(mantissa1: &mut U, mantissa2: &mut U) -> Result<U, ()>
    {
        assert!(!EXP_BASE.is_zero());
        if let Some(b) = U::from(EXP_BASE)
        {
            if (*mantissa1 % b).is_zero() && (*mantissa2 % b).is_zero()
            {
                *mantissa1 = *mantissa1/b;
                *mantissa2 = *mantissa2/b;
                return Ok(*mantissa1 + *mantissa2)
            }
        }

        if !Self::BASE_FACTORS_2.is_empty()
        {
            for &[base_final, base] in Self::BASE_FACTORS_2
            {
                if let Some(bf) = U::from(base_final) && let Some(b) = U::from(base)
                    && (*mantissa1 % b).is_zero() && (*mantissa2 % b).is_zero() 
                {
                    *mantissa1 = *mantissa1/b;
                    *mantissa2 = *mantissa2/b;
                    return Ok(util::rounding_div(*mantissa1 + *mantissa2, bf))
                }
            }
        }

        if !Self::BASE_FACTORS_2.is_empty()
        {   
            let (_, [bf, b]) = Self::BASE_FACTORS_2.iter()
                .filter_map(|&[base_final, base]| {
                    if let Some(bf) = U::from(base_final) && let Some(b) = U::from(base)
                        && let k = (*mantissa1 % b).is_zero() as u8 + (*mantissa2 % b).is_zero() as u8
                        && k > 0
                    {
                        Some((k, [bf, b]))
                    }
                    else
                    {
                        None
                    }
                }).reduce(|a, b| {
                    if a.0 < b.0
                    {
                        b
                    }
                    else
                    {
                        a
                    }
                }).unwrap();

            *mantissa1 = util::rounding_div(*mantissa1, b);
            *mantissa2 = util::rounding_div(*mantissa2, b);
            return Ok(util::rounding_div(*mantissa1 + *mantissa2, bf))
        }

        if let Some(b) = U::from(EXP_BASE)
        {
            *mantissa1 = util::rounding_div(*mantissa1, b);
            *mantissa2 = util::rounding_div(*mantissa2, b);
            return Ok(*mantissa1 + *mantissa2)
        }

        Err(())
    }

    /// Computes the arcsine of a number. Return value is in radians in
    /// the range [-pi/2, pi/2] or NaN if the number is outside the range
    /// [-1, 1].
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
    /// let f = FpDouble::FRAC_PI_2();
    ///
    /// // asin(sin(pi/2))
    /// let abs_difference = (f.sin().asin() - f).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-4));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn asin(self) -> Self
    {
        if self.is_nan()
        {
            return self
        }
        let xabs = self.abs();
        if xabs.is_one()
        {
            return Self::FRAC_PI_2().copysign(self)
        }
        if xabs > Self::one()
        {
            return Self::nan()
        }
        if xabs < Self::from(0.000000007450580596923828125)
        {
            return self
        }
        
        const N: usize = 10;
        const C: [f64; N] = [
            1.051231959,
            0.054946487,
            0.004080631,
            0.000407890,
            0.000046985,
            0.000005881,
            0.000000777,
            0.000000107,
            0.000000015,
            0.000000002
        ];

        static mut P: Option<[f64; N]> = None;
        let p = unsafe {
            #[allow(static_mut_refs)]
            P.get_or_insert_with(|| util::chebychev_approximation(C))
        };

        let one = Self::one();
        let w = if xabs <= Self::FRAC_1_SQRT_2()
        {
            self.abs()
        }
        else
        {
            (one - self*self).sqrt()
        };
        
        let ww2 = Self::from_int(2u8)*w*w;
        
        let mut y = {
            let z = Fps::from(ww2) - Fps::from(one);

            z.polynomial(p)
        }*w;
        if xabs > Self::FRAC_1_SQRT_2()
        {
            y = Self::FRAC_PI_2() - y
        }
        if self.is_sign_negative()
        {
            y = -y
        }

        /*const PIO2_HI: f64 = 1.57079637050628662109375;
        const PIO2_LO: f64 = -4.37113900018624283e-8;
        const PIO4_HI: f64 = 0.785398185253143310546875;

        const P0: f64 = 1.666675248e-1;
        const P1: f64 = 7.495297643e-2;
        const P2: f64 = 4.547037598e-2;
        const P3: f64 = 2.417951451e-2;
        const P4: f64 = 4.216630880e-2;

        let pio2_hi = Self::from(PIO2_HI);
        let pio2_lo = Self::from(PIO2_LO);
        let pio4_hi = Self::from(PIO4_HI);

        let mut y = if xabs.is_one()
        {
            self*pio2_hi + self*pio2_lo
        }
        else if xabs > Self::one()
        {
            return (self - self)/(self - self)
        }
        else if xabs < Self::from(0.5)
        {
            if xabs < Self::from(0.000000007450580596923828125)
            {
                self
            }
            else
            {
                let t = self*self;
                let w = t*(Self::from(P0) + t*(Self::from(P1) + t*(Self::from(P2) + t*(Self::from(P3) + t*Self::from(P4)))));
                self + self*w
            }
        }
        else
        {
            let mut w = Self::one() - xabs;
            let mut t = w*Self::from(0.5);
            let mut p = t*(Self::from(P0) + t*(Self::from(P1) + t*(Self::from(P2) + t*(Self::from(P3) + t*Self::from(P4)))));
            let s = t.sqrt();

            let two = Self::from_uint(2u8);

            if xabs >= Self::from(0.975)
            {
                t = pio2_hi - (two*(s + s*p) - pio2_lo);
            }
            else
            {
                let e_s = s.exp_bits();
                let mut f_s = s.frac_bits();
                if !Self::IS_INT_IMPLICIT
                {
                    f_s = f_s + (s.int_bits() << Self::INT_POS)
                }
                f_s = f_s & ((U::max_value() >> FRAC_SIZE/2) << FRAC_SIZE/2);
                w = Self::from_bits((e_s << Self::EXP_POS) + (f_s << Self::FRAC_POS));

	            let c = (t - w*w)/(s + w);
                let r = p;
                p = two*s*r - (pio2_lo - two*c);
                let q  = pio4_hi - two*w;
                t = pio4_hi - (p - q);
            }

            t.copysign(self)
        };*/

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let sin = y.sin();
                let cos = y.cos();
                let dsx = sin - self;
                let dy = dsx/cos;
                let yy = y - dy;
                if !yy.is_finite()
                {
                    break
                }
                y = yy;
            }
        }

        y
    }

    /// Computes the arccosine of a number. Return value is in radians in
    /// the range [0, pi] or NaN if the number is outside the range
    /// [-1, 1].
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
    /// let f = FpDouble::FRAC_PI_4();
    ///
    /// // acos(cos(pi/4))
    /// let abs_difference = (f.cos().acos() - f).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn acos(self) -> Self
    {
        const N: usize = 10;
        const C: [f64; N] = [
            1.051231959,
            0.054946487,
            0.004080631,
            0.000407890,
            0.000046985,
            0.000005881,
            0.000000777,
            0.000000107,
            0.000000015,
            0.000000002
        ];

        match self.classify()
        {
            FpCategory::Nan => self,
            FpCategory::Infinite => Self::nan().copysign(self),
            FpCategory::Zero => Self::FRAC_PI_2(),
            FpCategory::Normal | FpCategory::Subnormal => {
                let one = Self::one();
                let xabs = self.abs();
                if xabs == one
                {
                    return if self.is_sign_negative()
                    {
                        Self::PI()
                    }
                    else
                    {
                        Self::zero()
                    }
                }
                if xabs > one
                {
                    return Self::nan()
                }
                if xabs < Self::from(0.00000001490116119384765625)
                {
                    return Self::FRAC_PI_2() - self
                }
        
                static mut P: Option<[f64; N]> = None;
                let p = unsafe {
                    #[allow(static_mut_refs)]
                    P.get_or_insert_with(|| util::chebychev_approximation(C))
                }.map(Self::from);
        
                let w = if xabs <= Self::FRAC_1_SQRT_2()
                {
                    xabs
                }
                else
                {
                    (one - xabs.squared()).sqrt()
                };
                
                let ww = w.squared();
                let ww2 = ww + ww;
                let z = ww2 - one;
        
                let mut y = util::polynomial(&p, z, true, false)*w;
                if xabs <= Self::FRAC_1_SQRT_2()
                {
                    y = Self::FRAC_PI_2() - y
                }
                if self.is_sign_negative()
                {
                    y = Self::PI() - y
                }
                
                if !y.is_finite()
                {
                    let xx = xabs.squared();
                    y = if xx > Self::from(0.5)
                    {
                        ((Self::one() - xx).sqrt()/self).atan()
                    }
                    else
                    {
                        Self::FRAC_PI_2() - (self/(Self::one() - xx).sqrt()).atan()
                    };
                }
        
                /*const PIO2_HI: f64 = 1.5707962513e+00;
                const PIO2_LO: f64 = 7.5497894159e-08;
        
                const PS0: f64 = 1.6666667163e-01;
                const PS1: f64 = -3.2556581497e-01;
                const PS2: f64 = 2.0121252537e-01;
                const PS3: f64 = -4.0055535734e-02;
                const PS4: f64 = 7.9153501429e-04;
                const PS5: f64 = 3.4793309169e-05;
                const QS1: f64 = -2.4033949375e+00;
                const QS2: f64 = 2.0209457874e+00;
                const QS3: f64 = -6.8828397989e-01;
                const QS4: f64 = 7.7038154006e-02;
        
                let mut y = if xabs < Self::from(0.5)
                {
                    if xabs < Self::from(0.00000001490116119384765625)
                    {
                        Self::from(PIO2_HI) + Self::from(PIO2_LO)
                    }
                    else
                    {
                        let z = self*self;
                        let p = z*(Self::from(PS0) + z*(Self::from(PS1) + z*(Self::from(PS2) + z*(Self::from(PS3) + z*(Self::from(PS4) + z*Self::from(PS5))))));
                        let q = Self::one() + z*(Self::from(QS1) + z*(Self::from(QS2) + z*(Self::from(QS3) + z*Self::from(QS4))));
                        let r = p/q;
                        Self::from(PIO2_HI) - (self - (Self::from(PIO2_LO) - self*r))
                    }
                }
                else if self.is_sign_negative()
                {
                    let z = (Self::one() + self)*Self::from(0.5);
                    let p = z*(Self::from(PS0) + z*(Self::from(PS1) + z*(Self::from(PS2) + z*(Self::from(PS3) + z*(Self::from(PS4) + z*Self::from(PS5))))));
                    let q = Self::one() + z*(Self::from(QS1) + z*(Self::from(QS2) + z*(Self::from(QS3) + z*Self::from(QS4))));
                    let r = p/q;
                    let s = z.sqrt();
                    let w = r*s - Self::from(PIO2_LO);
                    Self::PI() - Self::from_uint(2u8)*(s + w)
                }
                else
                {
                    let z = (Self::one() - self)*Self::from(0.5);
                    let s = z.sqrt();
        
                    let e_s = s.exp_bits();
                    let mut f_s = s.frac_bits();
                    if !Self::IS_INT_IMPLICIT
                    {
                        f_s = f_s + (s.int_bits() << Self::INT_POS)
                    }
                    f_s = f_s & ((U::max_value() >> FRAC_SIZE/2) << FRAC_SIZE/2);
                    let df = Self::from_bits((e_s << Self::EXP_POS) + (f_s << Self::FRAC_POS));
                    let c = (z - df*df)/(s + df);
        
                    let p = z*(Self::from(PS0) + z*(Self::from(PS1) + z*(Self::from(PS2) + z*(Self::from(PS3) + z*(Self::from(PS4) + z*Self::from(PS5))))));
                    let q = Self::one() + z*(Self::from(QS1) + z*(Self::from(QS2) + z*(Self::from(QS3) + z*Self::from(QS4))));
                    let r = p/q;
                    let w = r*s + c;
                    Self::from_uint(2u8)*(df + w)
                };*/
        
                if y.is_finite()
                {
                    const NEWTON: usize = NEWTON_TRIG;
        
                    for _ in 0..NEWTON
                    {
                        let (cos, cos_s) = y.cos_extra_sign();
                        let (sin, sin_s) = y.sin_extra_sign();
                        let (dcx, dcx_s) = cos.add_with_sign_extra_sign(cos_s, self, true);
                        let s = !dcx_s^sin_s;
                        let dy = dcx/sin;
                        let yy = y.add_with_sign(false, dy, !s);
                        if !yy.is_finite()
                        {
                            break
                        }
                        y = yy;
                    }
                }
        
                y
            }
        }
    }

    /// Computes the arctangent of a number. Return value is in radians in the
    /// range [-pi/2, pi/2];
    ///
    /// This implementation is based on Harvey M. Wagner's [Polynomial approximations to elementary functions](https://www.ams.org/journals/mcom/1954-08-047/S0025-5718-1954-0063487-2/S0025-5718-1954-0063487-2.pdf).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::one();
    ///
    /// // atan(tan(1))
    /// let abs_difference = (f.tan().atan() - f).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-5));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn atan(self) -> Self
    {
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Zero => self,
            FpCategory::Infinite => Self::FRAC_PI_2().copysign(self),
            FpCategory::Normal | FpCategory::Subnormal => {
                const N: usize = 11;
                const C: [f64; N] = [
                    0.881373587,
                    -0.105892925,
                    0.011135843,
                    -0.001381195,
                    0.000185743,
                    -0.000026215,
                    0.000003821,
                    -0.00000057,
                    0.000000086,
                    -0.000000013,
                    0.000000002
                ];
        
                static mut P: Option<[f64; N]> = None;
                let p = unsafe {
                    #[allow(static_mut_refs)]
                    P.get_or_insert_with(|| util::chebychev_approximation(C))
                };
        
                let one = Self::one();
                let xabs = self.abs();
                let inv = xabs > one;
                let w = if !inv
                {
                    self
                }
                else
                {
                    self.recip()
                };
                
                let ww = w.squared();
                let ww2 = ww + ww;
                let mut s = false;
                let mut y = if !SIGN_BIT && ww2 < one
                {
                    let z = one - ww2;
        
                    let mut y = util::polynomial(p, z, SIGN_BIT, true);
                    if y.is_nan()
                    {
                        s = !s;
                        y = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, true);
                    }
                    y
                }
                else
                {
                    let z = ww2 - one;
        
                    let mut y = util::polynomial(p, z, SIGN_BIT, false);
                    if !SIGN_BIT && y.is_nan()
                    {
                        s = !s;
                        y = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
                    }
                    y
                }*w;

                let pi = Self::PI();
                let half_pi = Self::FRAC_PI_2();

                if inv
                {
                    if s
                    {
                        y = half_pi + y
                    }
                    else if !SIGN_BIT && y > half_pi
                    {
                        y %= pi;
                        y = pi + half_pi - y
                    }
                    else
                    {
                        y = half_pi - y
                    }
                }
                else if s
                {
                    y %= pi;
                    y = pi - y
                }
        
                /*if y.is_nan()
                {
                    const TAYLOR: usize = 8;
                    y = if self.abs() < Self::one()
                    {
                        let mut s = false;
                        let mut z = self;
                        let mut y = z;
        
                        for k in 1..TAYLOR
                        {
                            z *= self*self;
                            s = !s;
                            let dy = z/Self::from_uint(1 + 2*k);
                            if !s
                            {
                                y += dy
                            }
                            else
                            {
                                y -= dy
                            }
                        }
        
                        y
                    }
                    else
                    {
                        let mut s = false;
                        let mut z = Self::one()/self;
                        let mut y = Self::FRAC_PI_2() - z;
        
                        for k in 1..TAYLOR
                        {
                            z /= self*self;
                            s = !s;
                            let dy = z/Self::from_uint(1 + 2*k);
                            if !s
                            {
                                y -= dy
                            }
                            else
                            {
                                y += dy
                            }
                        }
        
                        y
                    }
                }*/
        
                y %= pi;
                // TODO: Is "if" enough?
                while y.abs() > half_pi
                {
                    y -= pi.copysign(y)
                }
                
                if y.is_finite()
                {
                    const NEWTON: usize = NEWTON_TRIG;
        
                    for _ in 0..NEWTON
                    {
                        let (sin, sin_s) = y.sin_extra_sign();
                        let (cos, cos_s) = y.cos_extra_sign();
                        let xcos = self*cos;
                        let (dsxc, dsxc_s) = sin.add_with_sign_extra_sign(sin_s, xcos, !cos_s);
                        let dy = dsxc*cos;
                        let s = dsxc_s^cos_s;
                        let yy = y.add_with_sign(false, dy, !s);
                        if !yy.is_finite()
                        {
                            break
                        }
                        y = yy;
                    }
                }
                
                if y.abs() > half_pi
                {
                    return half_pi.copysign(y)
                }
        
                y
            }
        }
    }

    /// Computes the four quadrant arctangent of `self` (`y`) and `other` (`x`).
    ///
    /// * `x = 0`, `y = 0`: `0`
    /// * `x >= 0`: `arctan(y/x)` -> `[-pi/2, pi/2]`
    /// * `y >= 0`: `arctan(y/x) + pi` -> `(pi/2, pi]`
    /// * `y < 0`: `arctan(y/x) - pi` -> `(-pi, -pi/2)`
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// // All angles from horizontal right (+x)
    /// // 45 deg counter-clockwise
    /// let x1 = FpDouble::from(3.0);
    /// let y1 = FpDouble::from(-3.0);
    ///
    /// // 135 deg clockwise
    /// let x2 = FpDouble::from(-3.0);
    /// let y2 = FpDouble::from(3.0);
    ///
    /// let abs_difference_1 = (y1.atan2(x1) - (-FpDouble::FRAC_PI_4())).abs();
    /// let abs_difference_2 = (y2.atan2(x2) - (FpDouble::PI() - FpDouble::FRAC_PI_4())).abs();
    ///
    /// assert!(abs_difference_1 < FpDouble::from(1e-5));
    /// assert!(abs_difference_2 < FpDouble::from(1e-5));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn atan2(self, other: Self) -> Self
    {
        if other.is_zero()
        {
            return if self.is_zero()
            {
                self
            }
            else
            {
                Self::FRAC_PI_2().copysign(self)
            }
        }
        let atan = (self/other).atan();
        if other.is_sign_positive()
        {
            return atan
        }
        atan + Self::PI().copysign(self)
    }

    /// Simultaneously computes the sine and cosine of the number, `x`. Returns `(sin(x), cos(x))`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::FRAC_PI_4();
    /// let f = x.sin_cos();
    ///
    /// let abs_difference_0 = (f.0 - x.sin()).abs();
    /// let abs_difference_1 = (f.1 - x.cos()).abs();
    ///
    /// assert!(abs_difference_0 < FpDouble::from(1e-10));
    /// assert!(abs_difference_1 < FpDouble::from(1e-10));
    /// ```
    #[must_use]
    pub fn sin_cos(self) -> (Self, Self)
    {
        // TODO: Can i do more before and after?
        (self.sin(), self.cos())
    }

    /// Returns `e^(self) - 1`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(7.0);
    ///
    /// // e^(ln(7)) - 1
    /// let abs_difference = (x.ln().exp_m1() - FpDouble::from(6.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-8));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn exp_m1(self) -> Self
    {
        if self.is_zero()
        {
            return Self::zero()
        }
        if self.is_sign_positive()
        {
            return self.exp() - Self::one()
        }

        let mut y = self.exp_nonewton() - Self::one();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_EXP;

            for _ in 0..NEWTON
            {
                let yp1 = y + Self::one();
                let x = y.ln_1p();
                let (diff, s) = x.sub_extra_sign(self);
                let dy = yp1*diff;
                let yy = y.add_with_sign(false, dy, !s);
                if !yy.is_finite()
                {
                    break
                }
                y = yy;
            }
        }

        y
    }

    /// Returns `ln(1+n)` (natural logarithm) more accurately than if
    /// the operations were performed separately.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let x = FpDouble::E() - FpDouble::one();
    ///
    /// // ln(1 + (e - 1)) == ln(e) == 1
    /// let abs_difference = (x.ln_1p() - FpDouble::one()).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-9));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn ln_1p(self) -> Self
    {
        let xp1 = self + Self::one();
        let mut y = xp1.ln_nonewton();

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_LN;

            let one = Self::one();
            for _ in 0..NEWTON
            {
                let xp1dx = xp1/y.exp();
                let (dy, s) = one.sub_extra_sign(xp1dx);
                let yy = y.add_with_sign(false, dy, !s);
                if !yy.is_finite()
                {
                    break
                }
                y = yy;
            }
        }

        y
    }

    /// Hyperbolic sine function.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let x = FpDouble::one();
    ///
    /// let f = x.sinh();
    /// // Solving sinh() at 1 gives `(e^2-1)/(2e)`
    /// let g = (e*e - FpDouble::one())/(FpDouble::from(2.0)*e);
    /// let abs_difference = (f - g).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn sinh(self) -> Self
    {
        if !self.is_finite()
        {
            return self
        }

        let emx = (-self.abs()).exp();
    
        let mut y = ((Self::one() - emx*emx)/emx*Self::from(0.5)).copysign(self);

        if y.is_nan()
        {
            let ex = (self.abs()).exp();
        
            y = ((ex*ex - Self::one())/ex*Self::from(0.5)).copysign(self);
        }

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let x = y.asinh();
                let (diff, s) = x.sub_extra_sign(self);
                let dy = diff*(y*y + Self::one()).sqrt();
                let yy = y.add_with_sign(false, dy, !s);
                if !yy.is_finite()
                {
                    break
                }
                y = yy;
            }
        }

        if y.is_nan()
        {
            return Self::infinity().copysign(self)
        }

        y
    }

    /// Hyperbolic cosine function.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let x = FpDouble::one();
    /// 
    /// let f = x.cosh();
    /// // Solving cosh() at 1 gives this result
    /// let g = (e*e + FpDouble::one())/(FpDouble::from(2.0)*e);
    /// let abs_difference = (f - g).abs();
    ///
    /// // Same result
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn cosh(self) -> Self
    {
        let xabs = self.abs();
        if !self.is_finite()
        {
            return xabs
        }

        let ex = xabs.exp();
        let mut y = ex.midpoint(ex.recip());

        if y.is_finite()
        {
            const NEWTON: usize = NEWTON_TRIG;

            for _ in 0..NEWTON
            {
                let x = y.acosh();
                let (diff, s) = x.sub_extra_sign(xabs);
                let dy = diff*(y.squared() - Self::one()).sqrt();
                let yy = y.add_with_sign(false, dy, !s);
                if !yy.is_finite()
                {
                    break
                }
                y = yy;
            }
        }

        y
    }
    
    /// Hyperbolic tangent function.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let x = FpDouble::one();
    ///
    /// let f = x.tanh();
    /// // Solving tanh() at 1 gives `(1 - e^(-2))/(1 + e^(-2))`
    /// let g = (FpDouble::one() - e.powi(-2))/(FpDouble::one() + e.powi(-2));
    /// let abs_difference = (f - g).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn tanh(self) -> Self
    {
        match self.classify()
        {
            FpCategory::Nan | FpCategory::Zero => self,
            FpCategory::Infinite => Self::one().copysign(self),
            FpCategory::Normal | FpCategory::Subnormal => {
                let one = Self::one();
                let xabs = self.abs();
                if xabs < Self::from(2.7755575615628913510590791702271e-17)
                {
                    return self*(one + self)
                }
        
                let mut y = if SIGN_BIT
                {
                    let ex = (-xabs).exp();
                    let ex2 = ex.squared();
                    let ex2p1 = one + ex2;
                    let ex2m1 = one - ex2;
                    
                    (ex2m1/ex2p1).copysign(self)
                }
                else
                {
                    let ex = xabs.exp().recip();
                    let ex2 = ex.squared();
                    let ex2p1 = one + ex2;
                    let ex2m1 = one - ex2;
                    
                    ex2m1/ex2p1
                };
        
                if y.is_finite()
                {
                    const NEWTON: usize = NEWTON_TRIG;
        
                    for _ in 0..NEWTON
                    {
                        let x = y.atanh();
                        if !x.is_finite()
                        {
                            break
                        }
                        let (diff, s) = x.sub_extra_sign(self);
                        let dy = diff*(one - y)*(one + y);
                        let yy = y.add_with_sign(false, dy, !s);
                        if !yy.is_finite()
                        {
                            break
                        }
                        y = yy;
                    }
                }
        
                // TODO: "abs larger than one"-function
                if y.abs() >= Self::one()
                {
                    return Self::one().copysign(y)
                }
                if y.is_nan()
                {
                    return Self::one().copysign(self)
                }
        
                y
            }
        }
    }

    /// Inverse hyperbolic sine function.
    ///
    /// This implementation is based on the glibc implementation of asinhf.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::one();
    /// let f = x.sinh().asinh();
    ///
    /// let abs_difference = (f - x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn asinh(self) -> Self
    {
        match self.classify()
        {
            FpCategory::Infinite | FpCategory::Nan | FpCategory::Zero => self,
            FpCategory::Subnormal | FpCategory::Normal => {
                //(self + (self*self + Self::one()).sqrt()).ln()

                let xabs = self.abs();
                if xabs < Self::from(0.00006103515625)
                {
                    return self
                }
                let w = if xabs > Self::from(16384.0)
                {
                    xabs.ln() + Self::LN_2()
                }
                else
                {
                    let one = Self::one();
                    let two = Self::from_uint(2u8);
                    let t = xabs*xabs;
                    if xabs > two
                    {
                        (xabs + xabs + ((t + one).sqrt() + xabs).recip()).ln()
                    }
                    else
                    {
                        (xabs + t/(one + (one + t).sqrt())).ln_1p()
                    }
                };
                w.copysign(self)
            }
        }
    }

    /// Inverse hyperbolic cosine function.
    /// 
    /// This implementation is based on the glibc implementation of acoshf.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::one();
    /// let f = x.cosh().acosh();
    ///
    /// let abs_difference = (f - x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn acosh(self) -> Self
    {
        match self.classify()
        {
            FpCategory::Nan => self,
            FpCategory::Infinite => if self.is_sign_negative()
            {
                Self::nan()
            }
            else
            {
                self
            },
            FpCategory::Zero => Self::nan(),
            FpCategory::Normal | FpCategory::Subnormal => if self < Self::one()
            {
                Self::snan()
            }
            else
            {
                //(self + (self*self - Self::one()).sqrt()).ln()
        
                if self > Self::from(268435456.0)
                {
                    self.ln() + Self::LN_2()
                }
                else if self.is_one()
                {
                    Self::zero()
                }
                else
                {
                    let two = Self::from_uint(2u8);
                    let one = Self::one();
                    if self > two
                    {
                        let t = self*self;
                        (two*self - (self + (t - one).sqrt()).recip()).ln()
                    }
                    else
                    {
                        let t = self - one;
                        (t + (two*t + t*t).sqrt()).ln_1p()
                    }
                }
            }
        }
    }

    /// Inverse hyperbolic tangent function.
    ///
    /// This implementation is based on the glibc implementation of atanhf.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// use num::traits::FloatConst;
    ///
    /// let e = FpDouble::E();
    /// let f = e.tanh().atanh();
    ///
    /// let abs_difference = (f - e).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1.0e-2));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn atanh(self) -> Self
    {
        let y = match self.classify()
        {
            FpCategory::Nan | FpCategory::Zero => return self,
            FpCategory::Infinite => Self::nan(),
            FpCategory::Normal | FpCategory::Subnormal => {
                let xabs = self.abs();
                let one = Self::one();
                match xabs.total_cmp(one)
                {
                    Ordering::Greater => Self::nan(),
                    Ordering::Equal => Self::infinity(),
                    Ordering::Less => {
                        //Self::from(0.5)*((Self::one() + self.abs())/(Self::one() - self.abs())).ln().copysign(self)
                        
                        let half = Self::from(0.5);
                        let t = xabs + xabs;
                        if t < one
                        {
                            if xabs < Self::from(0.0000000037252902984619140625)
                            {
                                return self
                            }
                
                            half*(t + t*xabs/(one - xabs)).ln_1p()
                        }
                        else
                        {
                            half*(t/(one - xabs)).ln_1p()
                        }
                    }
                }
            }
        };
        y.copysign(self)
    }
    
    fn ln_gamma_lanczos(self) -> Self
    {
        const LANCZOS_CHEB_7: [f64; 9] = [
            0.99999999999980993227684700473478,
            676.520368121885098567009190444019,
            -1259.13921672240287047156078755283,
            771.3234287776530788486528258894,
            -176.61502916214059906584551354,
            12.507343278686904814458936853,
            -0.13857109526572011689554707,
            9.984369578019570859563e-6,
            1.50563273514931155834e-7
        ];

        const LOGROOT2PI: f64 = 0.9189385332046727417803297364056176;

        let mut sum = Self::from(LANCZOS_CHEB_7[0]);
        for (k, c) in LANCZOS_CHEB_7.iter()
            .skip(1)
            .enumerate()
        {
            if SIGN_BIT || c.is_sign_positive()
            {
                sum += Self::from(*c)/(self + Self::from_uint(k))
            }
        }
        if !SIGN_BIT
        {
            for (k, c) in LANCZOS_CHEB_7.iter()
                .skip(1)
                .enumerate()
            {
                if c.is_sign_negative()
                {
                    sum -= Self::from(-c)/(self + Self::from_uint(k))
                }
            }
        }

        let term1 = (self - Self::from(0.5))
            *((self + Self::from(6.5))/Self::E()).ln();
        let term2 = Self::from(LOGROOT2PI) + sum.ln();

        let seven = Self::from(7u8);
        if term2 < seven
        {
            return term1 + term2 - seven
        }
        term1 + (term2 - seven)
    }

    /// Natural logarithm of the absolute value of the gamma function
    ///
    /// The integer part of the tuple indicates the sign of the gamma function.
    ///
    /// This implementation is based on [the libstdc++ implementation of the gamma function](https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.4/a01203.html).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let x = FpDouble::from(2.0);
    ///
    /// let abs_difference = (x.ln_gamma().0 - FpDouble::zero()).abs();
    ///
    /// assert!(abs_difference <= FpDouble::from(1e-2));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn ln_gamma(self) -> (Self, i32)
    {
        match self.classify()
        {
            FpCategory::Nan => (self.abs(), 0),
            FpCategory::Infinite => (self.abs(), if self.is_sign_negative() {0} else {1}),
            FpCategory::Zero => (Self::infinity(), if self.is_sign_negative() {-1} else {1}),
            FpCategory::Normal | FpCategory::Subnormal => {
                if self >= Self::from(0.5)
                {
                    return (self.ln_gamma_lanczos(), 1)
                }
        
                let sin_fact = (Self::PI()*self).sin();
        
                if sin_fact.is_zero()
                {
                    return (Self::infinity(), 0)
                }
        
                let sin_fact_abs = sin_fact.abs();
                let ln_pi = Self::PI().ln();
        
                (
                    if !SIGN_BIT && sin_fact_abs < Self::one()
                    {
                        ln_pi + sin_fact_abs.recip().ln()
                    }
                    else
                    {
                        ln_pi - sin_fact_abs.ln()
                    } - (Self::one() - self).ln_gamma_lanczos(),
                    if sin_fact.is_sign_negative() {-1} else {1}
                )
            }
        }
    }

    /// Gamma function.
    /// 
    /// This implementation is based on [the libstdc++ implementation of the gamma function](https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.4/a01203.html).
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    /// 
    /// let x = FpDouble::from(5.0f32);
    ///
    /// let abs_difference = (x.gamma() - FpDouble::from(24.0)).abs();
    ///
    /// assert!(abs_difference <= FpDouble::from(1e-1));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn gamma(self) -> Self
    {
        let (lgamma, sign) = self.ln_gamma();
        if sign == 0
        {
            return Self::nan()
        }
        lgamma.exp().with_sign(sign < 0)
    }

    /// Returns a number composed of the magnitude of `self` and the sign of
    /// `sign`.
    ///
    /// Equal to `self` if the sign of `self` and `sign` are the same, otherwise
    /// equal to `-self`. If `self` is a NaN, then a NaN with the sign bit of
    /// `sign` is returned. Note, however, that conserving the sign bit on NaN
    /// across arithmetical operations is not generally guaranteed.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.5);
    /// let s = FpDouble::from(0.42);
    ///
    /// assert_eq!(f.copysign(s), f);
    /// assert_eq!(f.copysign(-s), -f);
    /// assert_eq!((-f).copysign(s), f);
    /// assert_eq!((-f).copysign(-s), -f);
    ///
    /// assert!(FpDouble::nan().copysign(FpDouble::one()).is_nan());
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn copysign(self, sign: Self) -> Self
    {
        if !SIGN_BIT
        {
            return self
        }

        let mask = Self::sign_mask();
        Self::from_bits((self.to_bits() & (!mask)) | (sign.to_bits() & mask))
    }

    fn with_sign(self, sign: bool) -> Self
    {
        if !SIGN_BIT
        {
            if !sign || self.is_zero()
            {
                return self
            }
            return Self::qnan()
        }
        
        let mask = Self::sign_mask();
        let mut bits = self.to_bits();
        if sign
        {
            bits = bits | mask;
        }
        else
        {
            bits = bits & !mask
        }
        
        Self::from_bits(bits)
    }

    fn xor_sign(self, sign: bool) -> Self
    {
        if !sign
        {
            return self
        }
        if !SIGN_BIT
        {
            if self.is_zero()
            {
                return self
            }
            return Self::qnan()
        }
        
        let mask = Self::sign_mask();
        let mut bits = self.to_bits();
        bits = bits ^ mask;
        Self::from_bits(bits)
    }

    fn erfc1(self) -> Self
    {
        const ERX: f64 = 8.45062911510467529297e-01; /* 0x3FEB0AC1, 0x60000000 */
        
        /*
        * Coefficients for approximation to  erf  in [0.84375,1.25]
        */
        const PA: [f64; 7] = [
            -2.36211856075265944077e-03, /* 0xBF6359B8, 0xBEF77538 */
            4.14856118683748331666e-01, /* 0x3FDA8D00, 0xAD92B34D */
            -3.72207876035701323847e-01, /* 0xBFD7D240, 0xFBB8C3F1 */
            3.18346619901161753674e-01, /* 0x3FD45FCA, 0x805120E4 */
            -1.10894694282396677476e-01, /* 0xBFBC6398, 0x3D3E28EC */
            3.54783043256182359371e-02, /* 0x3FA22A36, 0x599795EB */
            -2.16637559486879084300e-03 /* 0xBF61BF38, 0x0A96073F */
        ];
        const QA: [f64; 6] = [
            1.06420880400844228286e-01, /* 0x3FBB3E66, 0x18EEE323 */
            5.40397917702171048937e-01, /* 0x3FE14AF0, 0x92EB6F33 */
            7.18286544141962662868e-02, /* 0x3FB2635C, 0xD99FE9A7 */
            1.26171219808761642112e-01, /* 0x3FC02660, 0xE763351F */
            1.36370839120290507362e-02, /* 0x3F8BEDC2, 0x6B51DD1C */
            1.19844998467991074170e-02 /* 0x3F888B54, 0x5735151D */
        ];
        let one = Self::one();
        let xabs = self.abs();
        let (s, b) = xabs.sub_extra_sign(one);

        // TODO: CLEAN THIS UP!
        let mut r = false;
        let mut p = util::polynomial(&PA, s, SIGN_BIT, b);
        if !SIGN_BIT && p.is_nan()
        {
            p = util::polynomial(&PA.map(Neg::neg), s, SIGN_BIT, b);
            r = !r;
        }
        let q = if b
        {
            let mut q = s*util::polynomial(&QA, s, SIGN_BIT, b);
            if !SIGN_BIT && q.is_nan()
            {
                q = s*util::polynomial(&QA.map(Neg::neg), s, SIGN_BIT, b);
                one + q
            }
            else if SIGN_BIT || one >= q
            {
                one - q
            }
            else
            {
                r = !r;
                q - one
            }
        }
        else
        {
            let mut q = s*util::polynomial(&QA, s, SIGN_BIT, b);
            if SIGN_BIT || !q.is_nan()
            {
                one + q
            }
            else
            {
                q = s*util::polynomial(&QA.map(Neg::neg), s, SIGN_BIT, b);
                if one >= q
                {
                    one - q
                }
                else
                {
                    r = !r;
                    q - one
                }
            }
        };
    
        (one - Self::from(ERX)).add_with_sign(false, p/q, !r)
    }
    
    fn erfc2(mut self) -> Self
    {
        /*
         * Coefficients for approximation to  erfc in [1.25,1/0.35]
         */
        const RA: [f64; 8] = [
            -9.86494403484714822705e-03, /* 0xBF843412, 0x600D6435 */
            -6.93858572707181764372e-01, /* 0xBFE63416, 0xE4BA7360 */
            -1.05586262253232909814e+01, /* 0xC0251E04, 0x41B0E726 */
            -6.23753324503260060396e+01, /* 0xC04F300A, 0xE4CBA38D */
            -1.62396669462573470355e+02, /* 0xC0644CB1, 0x84282266 */
            -1.84605092906711035994e+02, /* 0xC067135C, 0xEBCCABB2 */
            -8.12874355063065934246e+01, /* 0xC0545265, 0x57E4D2F2 */
            -9.81432934416914548592e+00 /* 0xC023A0EF, 0xC69AC25C */
        ];
        const SA: [f64; 8] = [
            1.96512716674392571292e+01, /* 0x4033A6B9, 0xBD707687 */
            1.37657754143519042600e+02, /* 0x4061350C, 0x526AE721 */
            4.34565877475229228821e+02, /* 0x407B290D, 0xD58A1A71 */
            6.45387271733267880336e+02, /* 0x40842B19, 0x21EC2868 */
            4.29008140027567833386e+02, /* 0x407AD021, 0x57700314 */
            1.08635005541779435134e+02, /* 0x405B28A3, 0xEE48AE2C */
            6.57024977031928170135e+00, /* 0x401A47EF, 0x8E484A93 */
            -6.04244152148580987438e-02 /* 0xBFAEEFF2, 0xEE749A62 */
        ];
        /*
         * Coefficients for approximation to  erfc in [1/.35,28]
         */
        const RB: [f64; 7] = [
            -9.86494292470009928597e-03, /* 0xBF843412, 0x39E86F4A */
            -7.99283237680523006574e-01, /* 0xBFE993BA, 0x70C285DE */
            -1.77579549177547519889e+01, /* 0xC031C209, 0x555F995A */
            -1.60636384855821916062e+02, /* 0xC064145D, 0x43C5ED98 */
            -6.37566443368389627722e+02, /* 0xC083EC88, 0x1375F228 */
            -1.02509513161107724954e+03, /* 0xC0900461, 0x6A2E5992 */
            -4.83519191608651397019e+02 /* 0xC07E384E, 0x9BDC383F */
        ];
        const SB: [f64; 7] = [
            3.03380607434824582924e+01, /* 0x403E568B, 0x261D5190 */
            3.25792512996573918826e+02, /* 0x40745CAE, 0x221B9F0A */
            1.53672958608443695994e+03, /* 0x409802EB, 0x189D5118 */
            3.19985821950859553908e+03, /* 0x40A8FFB7, 0x688C246A */
            2.55305040643316442583e+03, /* 0x40A3F219, 0xCEDF3BE6 */
            4.74528541206955367215e+02, /* 0x407DA874, 0xE79FE763 */
            -2.24409524465858183362e+01 /* 0xC03670E2, 0x42712D62 */
        ];

        if self.abs() < Self::from(1.25)
        {
            /* |x| < 1.25 */
            return self.erfc1();
        }

        self = self.abs();
        let mut q_s = false;
        let s = (self*self).recip();
        let mut r;
        let mut big_s;
        let one = Self::one();
        if self < Self::from(1.0/0.35)
        {
            /* |x| < 1/.35 ~ 2.85714 */
            r = util::polynomial(&RA, s, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                r = util::polynomial(&RA.map(Neg::neg), s, SIGN_BIT, false);
                q_s = !q_s;
            }
            big_s = one + s*util::polynomial(&SA, s, SIGN_BIT, false);
            if !SIGN_BIT && big_s.is_nan()
            {
                big_s = s*util::polynomial(&SA.map(Neg::neg), s, SIGN_BIT, false) - one;
                q_s = !q_s;
            }
        }
        else
        {
            /* |x| > 1/.35 */
            r = util::polynomial(&RB, s, SIGN_BIT, false);
            if !SIGN_BIT && r.is_nan()
            {
                r = util::polynomial(&RB.map(Neg::neg), s, SIGN_BIT, false);
                q_s = !q_s;
            }
            big_s = one + s*util::polynomial(&SB, s, SIGN_BIT, false);
            if !SIGN_BIT && big_s.is_nan()
            {
                big_s = s*util::polynomial(&SB.map(Neg::neg), s, SIGN_BIT, false) - one;
                q_s = !q_s;
            }
        }
        let z = Self::from_bits((self.to_bits() >> (FRAC_SIZE/2)) << (FRAC_SIZE/2));

        let ofs = Self::from(0.5625);
        let div = (z * z + ofs).exp();
        let (p_diff, p_s) = z.sub_extra_sign(self);
        let p = p_diff*(z + self);
        let q = r / big_s;
        let (diff, s) = p.add_with_sign_extra_sign(p_s, q, q_s);
        let mut diff_exp = diff.exp();
        if s
        {
            diff_exp = diff_exp.recip()
        }
        diff_exp/div/self
    }
        
    /// Error function (f64)
    ///
    /// Calculates an approximation to the “error function”, which estimates
    /// the probability that an observation will fall within x standard
    /// deviations of the mean (assuming a normal distribution)
    pub fn erf(self) -> Self
    {
        /*
        * Coefficients for approximation to  erf on [0,0.84375]
        */
        const EFX8: f64 = 1.02703333676410069053e+00; /* 0x3FF06EBA, 0x8214DB69 */
        const PP: [f64; 5] = [
            1.28379167095512558561e-01, /* 0x3FC06EBA, 0x8214DB68 */
            -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
            -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
            -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
            -2.37630166566501626084e-05 /* 0xBEF8EAD6, 0x120016AC */
        ];
        const QQ: [f64; 5] = [
            3.97917223959155352819e-01, /* 0x3FD97779, 0xCDDADC09 */
            6.50222499887672944485e-02, /* 0x3FB0A54C, 0x5536CEBA */
            5.08130628187576562776e-03, /* 0x3F74D022, 0xC4D36B0F */
            1.32494738004321644526e-04, /* 0x3F215DC9, 0x221C1A10 */
            -3.96022827877536812320e-06 /* 0xBED09C43, 0x42A26120 */
        ];

        match self.classify()
        {
            FpCategory::Nan | FpCategory::Zero => self,
            FpCategory::Infinite => Self::one().copysign(self),
            FpCategory::Subnormal | FpCategory::Normal => {
                let xabs = self.abs();
                let one = Self::one();
                
                if xabs < Self::from(0.84375)
                {
                    /* |x| < 0.84375 */
                    if xabs < Self::from(0.0000000037252902984619140625)
                    {
                        /* |x| < 2**-28 */
                        /* avoid underflow */
                        return Self::from(0.125)*(Self::from(8u8)*self + Self::from(EFX8)*self);
                    }
                    let z = self*self;
                    let mut b = false;
                    let mut numer = util::polynomial(&PP, z, SIGN_BIT, false);
                    let mut denom = one + z*util::polynomial(&QQ, z, SIGN_BIT, false);
                    if !SIGN_BIT && numer.is_nan()
                    {
                        numer = util::polynomial(&PP.map(Neg::neg), z, SIGN_BIT, false);
                        b = !b;
                    }
                    if !SIGN_BIT && denom.is_nan()
                    {
                        denom = z*util::polynomial(&QQ.map(Neg::neg), z, SIGN_BIT, false) - one;
                        b = !b;
                    }
                    let y = numer / denom;
                    return self.add_with_sign(false, self*y, b)
                }

                let s = self.is_sign_negative();

                let y = one - if xabs < Self::from_uint(6u8)
                {
                    /* 0.84375 <= |x| < 6 */
                    self.erfc2()
                }
                else
                {
                    Self::from(f64::from_bits(0x0010000000000000))
                };
            
                y.with_sign(s)
            }
        }
    }

    /// Complementary error function
    ///
    /// Calculates the complementary probability.
    /// Is `1 - erf(x)`. Is computed directly, so that you can use it to avoid
    /// the loss of precision that would result from subtracting
    /// large probabilities (on large `x`) from 1.
    pub fn erfc(self) -> Self
    {
        /*
        * Coefficients for approximation to  erf on [0,0.84375]
        */
        const PP: [f64; 5] = [
            1.28379167095512558561e-01, /* 0x3FC06EBA, 0x8214DB68 */
            -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
            -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
            -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
            -2.37630166566501626084e-05 /* 0xBEF8EAD6, 0x120016AC */
        ];
        const QQ: [f64; 5] = [
            3.97917223959155352819e-01, /* 0x3FD97779, 0xCDDADC09 */
            6.50222499887672944485e-02, /* 0x3FB0A54C, 0x5536CEBA */
            5.08130628187576562776e-03, /* 0x3F74D022, 0xC4D36B0F */
            1.32494738004321644526e-04, /* 0x3F215DC9, 0x221C1A10 */
            -3.96022827877536812320e-06 /* 0xBED09C43, 0x42A26120 */
        ];

        let sign = self.is_sign_negative();
        let one = Self::one();
        let two = Self::from_uint(2u8); //TODO try "one + one"

        match self.classify()
        {
            FpCategory::Nan => self,
            FpCategory::Infinite => if sign
            {
                two
            }
            else
            {
                Self::zero()
            },
            FpCategory::Zero => one,
            FpCategory::Normal | FpCategory::Subnormal => {
                let xabs = self.abs();
                if xabs < Self::from(0.84375)
                {
                    /* |x| < 0.84375 */
                    if xabs < Self::from(1.3877787807814456755295395851135e-17)
                    {
                        /* |x| < 2**-56 */
                        return one - self;
                    }
                    let z = self*self;
                    let mut b = false;

                    // TODO: Make cleaner
                    let mut numer = util::polynomial(&PP, z, SIGN_BIT, false);
                    let mut denom = one + z*util::polynomial(&QQ, z, SIGN_BIT, false);
                    if !SIGN_BIT && numer.is_nan()
                    {
                        numer = util::polynomial(&PP.map(Neg::neg), z, SIGN_BIT, false);
                        b = !b;
                    }
                    if !SIGN_BIT && denom.is_nan()
                    {
                        denom = z*util::polynomial(&QQ.map(Neg::neg), z, SIGN_BIT, false) - one;
                        b = !b;
                    }

                    let y = numer / denom;
                    let xy = self*y;
                    let half = Self::from(0.5);
                    if sign || xabs < Self::from(1.0/4.0) || (!SIGN_BIT && self < half)
                    {
                        /* x < 1/4 */
                        let (diff, s_diff) = self.add_with_sign_extra_sign(false, xy, b);
                        return one.add_with_sign(false, diff, !s_diff)
                    }
                    let (diff_xhalf, s_xhalf) = self.sub_extra_sign(half);
                    let (diff, s_diff) = diff_xhalf.add_with_sign_extra_sign(s_xhalf, xy, b);
                    return half.add_with_sign(false, diff, !s_diff)
                }
                if xabs < Self::from_uint(28u8)
                {
                    /* 0.84375 <= |x| < 28 */
                    let y = self.erfc2();
                    if sign
                    {
                        return two - y;
                    }
                    return y;
                }
        
                let x1p_1022 = Self::from(f64::from_bits(0x0010000000000000));
                if sign
                {
                    two - x1p_1022
                }
                else
                {
                    x1p_1022 * x1p_1022
                }
            }
        }
    }
    
    fn bessel0_p(self) -> (Self, bool)
    {
        /* The asymptotic expansions of pzero is
        *      1 - 9/128 s^2 + 11025/98304 s^4 - ...,  where s = 1/x.
        * For x >= 2, We approximate pzero by
        *      pzero(x) = 1 + (R/S)
        * where  R = pR0 + pR1*s^2 + pR2*s^4 + ... + pR5*s^10
        *        S = 1 + pS0*s^2 + ... + pS4*s^10
        * and
        *      | pzero(x)-1-R/S | <= 2  ** ( -60.26)
        */
        const PR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00,  /* 0x00000000, 0x00000000 */
            -7.03124999999900357484e-02, /* 0xBFB1FFFF, 0xFFFFFD32 */
            -8.08167041275349795626e+00, /* 0xC02029D0, 0xB44FA779 */
            -2.57063105679704847262e+02, /* 0xC0701102, 0x7B19E863 */
            -2.48521641009428822144e+03, /* 0xC0A36A6E, 0xCD4DCAFC */
            -5.25304380490729545272e+03, /* 0xC0B4850B, 0x36CC643D */
        ];
        const PS8: [f64; 5] = [
            1.16534364619668181717e+02, /* 0x405D2233, 0x07A96751 */
            3.83374475364121826715e+03, /* 0x40ADF37D, 0x50596938 */
            4.05978572648472545552e+04, /* 0x40E3D2BB, 0x6EB6B05F */
            1.16752972564375915681e+05, /* 0x40FC810F, 0x8F9FA9BD */
            4.76277284146730962675e+04, /* 0x40E74177, 0x4F2C49DC */
        ];

        const PR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            -1.14125464691894502584e-11, /* 0xBDA918B1, 0x47E495CC */
            -7.03124940873599280078e-02, /* 0xBFB1FFFF, 0xE69AFBC6 */
            -4.15961064470587782438e+00, /* 0xC010A370, 0xF90C6BBF */
            -6.76747652265167261021e+01, /* 0xC050EB2F, 0x5A7D1783 */
            -3.31231299649172967747e+02, /* 0xC074B3B3, 0x6742CC63 */
            -3.46433388365604912451e+02, /* 0xC075A6EF, 0x28A38BD7 */
        ];
        const PS5: [f64; 5] = [
            6.07539382692300335975e+01, /* 0x404E6081, 0x0C98C5DE */
            1.05125230595704579173e+03, /* 0x40906D02, 0x5C7E2864 */
            5.97897094333855784498e+03, /* 0x40B75AF8, 0x8FBE1D60 */
            9.62544514357774460223e+03, /* 0x40C2CCB8, 0xFA76FA38 */
            2.40605815922939109441e+03, /* 0x40A2CC1D, 0xC70BE864 */
        ];

        const PR3: [f64; 6] = [
            /* for x in [4.547,2.8571]=1/[0.2199,0.35001] */
            -2.54704601771951915620e-09, /* 0xBE25E103, 0x6FE1AA86 */
            -7.03119616381481654654e-02, /* 0xBFB1FFF6, 0xF7C0E24B */
            -2.40903221549529611423e+00, /* 0xC00345B2, 0xAEA48074 */
            -2.19659774734883086467e+01, /* 0xC035F74A, 0x4CB94E14 */
            -5.80791704701737572236e+01, /* 0xC04D0A22, 0x420A1A45 */
            -3.14479470594888503854e+01, /* 0xC03F72AC, 0xA892D80F */
        ];
        const PS3: [f64; 5] = [
            3.58560338055209726349e+01, /* 0x4041ED92, 0x84077DD3 */
            3.61513983050303863820e+02, /* 0x40769839, 0x464A7C0E */
            1.19360783792111533330e+03, /* 0x4092A66E, 0x6D1061D6 */
            1.12799679856907414432e+03, /* 0x40919FFC, 0xB8C39B7E */
            1.73580930813335754692e+02, /* 0x4065B296, 0xFC379081 */
        ];

        const PR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            -8.87534333032526411254e-08, /* 0xBE77D316, 0xE927026D */
            -7.03030995483624743247e-02, /* 0xBFB1FF62, 0x495E1E42 */
            -1.45073846780952986357e+00, /* 0xBFF73639, 0x8A24A843 */
            -7.63569613823527770791e+00, /* 0xC01E8AF3, 0xEDAFA7F3 */
            -1.11931668860356747786e+01, /* 0xC02662E6, 0xC5246303 */
            -3.23364579351335335033e+00, /* 0xC009DE81, 0xAF8FE70F */
        ];
        const PS2: [f64; 5] = [
            2.22202997532088808441e+01, /* 0x40363865, 0x908B5959 */
            1.36206794218215208048e+02, /* 0x4061069E, 0x0EE8878F */
            2.70470278658083486789e+02, /* 0x4070E786, 0x42EA079B */
            1.53875394208320329881e+02, /* 0x40633C03, 0x3AB6FAFF */
            1.46576176948256193810e+01, /* 0x402D50B3, 0x44391809 */
        ];

        let p: &[f64; 6];
        let q: &[f64; 5];

        let xabs = self.abs();
        if xabs >= Self::from_uint(8u8)
        {
            p = &PR8;
            q = &PS8;
        }
        else if xabs >= Self::from(4.5454)
        {
            p = &PR5;
            q = &PS5;
        }
        else if xabs >= Self::from(2.857)
        {
            p = &PR3;
            q = &PS3;
        }
        else
        {
            p = &PR2;
            q = &PS2;
        }
        let one = Self::one();
        let z = (self*self).recip();
        let mut b = false;
        let mut r = util::polynomial(p, z, SIGN_BIT, false);
        let mut s = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && r.is_nan()
        {
            b = !b;
            r = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && s.is_nan()
        {
            s = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false);
            if s >= one
            {
                b = !b;
                s -= one
            }
            else
            {
                s = one - s
            }
        }
        
        let rs = r/s;
        if b
        {
            if one < rs
            {
                return (rs - one, true)
            }
            return (one - rs, false)
        }
        (one + rs, false)
    }
    
    fn bessel0_q(self) -> (Self, bool)
    {
        /* For x >= 8, the asymptotic expansions of qzero is
        *      -1/8 s + 75/1024 s^3 - ..., where s = 1/x.
        * We approximate pzero by
        *      qzero(x) = s*(-1.25 + (R/S))
        * where  R = qR0 + qR1*s^2 + qR2*s^4 + ... + qR5*s^10
        *        S = 1 + qS0*s^2 + ... + qS5*s^12
        * and
        *      | qzero(x)/s +1.25-R/S | <= 2  ** ( -61.22)
        */
        const QR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00, /* 0x00000000, 0x00000000 */
            7.32421874999935051953e-02, /* 0x3FB2BFFF, 0xFFFFFE2C */
            1.17682064682252693899e+01, /* 0x40278952, 0x5BB334D6 */
            5.57673380256401856059e+02, /* 0x40816D63, 0x15301825 */
            8.85919720756468632317e+03, /* 0x40C14D99, 0x3E18F46D */
            3.70146267776887834771e+04, /* 0x40E212D4, 0x0E901566 */
        ];
        const QS8: [f64; 6] = [
            1.63776026895689824414e+02,  /* 0x406478D5, 0x365B39BC */
            8.09834494656449805916e+03,  /* 0x40BFA258, 0x4E6B0563 */
            1.42538291419120476348e+05,  /* 0x41016652, 0x54D38C3F */
            8.03309257119514397345e+05,  /* 0x412883DA, 0x83A52B43 */
            8.40501579819060512818e+05,  /* 0x4129A66B, 0x28DE0B3D */
            -3.43899293537866615225e+05, /* 0xC114FD6D, 0x2C9530C5 */
        ];

        const QR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            1.84085963594515531381e-11, /* 0x3DB43D8F, 0x29CC8CD9 */
            7.32421766612684765896e-02, /* 0x3FB2BFFF, 0xD172B04C */
            5.83563508962056953777e+00, /* 0x401757B0, 0xB9953DD3 */
            1.35111577286449829671e+02, /* 0x4060E392, 0x0A8788E9 */
            1.02724376596164097464e+03, /* 0x40900CF9, 0x9DC8C481 */
            1.98997785864605384631e+03, /* 0x409F17E9, 0x53C6E3A6 */
        ];
        const QS5: [f64; 6] = [
            8.27766102236537761883e+01,  /* 0x4054B1B3, 0xFB5E1543 */
            2.07781416421392987104e+03,  /* 0x40A03BA0, 0xDA21C0CE */
            1.88472887785718085070e+04,  /* 0x40D267D2, 0x7B591E6D */
            5.67511122894947329769e+04,  /* 0x40EBB5E3, 0x97E02372 */
            3.59767538425114471465e+04,  /* 0x40E19118, 0x1F7A54A0 */
            -5.35434275601944773371e+03, /* 0xC0B4EA57, 0xBEDBC609 */
        ];

        const QR3: [f64; 6] = [
            /* for x in [4.547,2.8571]=1/[0.2199,0.35001] */
            4.37741014089738620906e-09, /* 0x3E32CD03, 0x6ADECB82 */
            7.32411180042911447163e-02, /* 0x3FB2BFEE, 0x0E8D0842 */
            3.34423137516170720929e+00, /* 0x400AC0FC, 0x61149CF5 */
            4.26218440745412650017e+01, /* 0x40454F98, 0x962DAEDD */
            1.70808091340565596283e+02, /* 0x406559DB, 0xE25EFD1F */
            1.66733948696651168575e+02, /* 0x4064D77C, 0x81FA21E0 */
        ];
        const QS3: [f64; 6] = [
            4.87588729724587182091e+01,  /* 0x40486122, 0xBFE343A6 */
            7.09689221056606015736e+02,  /* 0x40862D83, 0x86544EB3 */
            3.70414822620111362994e+03,  /* 0x40ACF04B, 0xE44DFC63 */
            6.46042516752568917582e+03,  /* 0x40B93C6C, 0xD7C76A28 */
            2.51633368920368957333e+03,  /* 0x40A3A8AA, 0xD94FB1C0 */
            -1.49247451836156386662e+02, /* 0xC062A7EB, 0x201CF40F */
        ];

        const QR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            1.50444444886983272379e-07, /* 0x3E84313B, 0x54F76BDB */
            7.32234265963079278272e-02, /* 0x3FB2BEC5, 0x3E883E34 */
            1.99819174093815998816e+00, /* 0x3FFFF897, 0xE727779C */
            1.44956029347885735348e+01, /* 0x402CFDBF, 0xAAF96FE5 */
            3.16662317504781540833e+01, /* 0x403FAA8E, 0x29FBDC4A */
            1.62527075710929267416e+01, /* 0x403040B1, 0x71814BB4 */
        ];
        const QS2: [f64; 6] = [
            3.03655848355219184498e+01,  /* 0x403E5D96, 0xF7C07AED */
            2.69348118608049844624e+02,  /* 0x4070D591, 0xE4D14B40 */
            8.44783757595320139444e+02,  /* 0x408A6645, 0x22B3BF22 */
            8.82935845112488550512e+02,  /* 0x408B977C, 0x9C5CC214 */
            2.12666388511798828631e+02,  /* 0x406A9553, 0x0E001365 */
            -5.31095493882666946917e+00, /* 0xC0153E6A, 0xF8B32931 */
        ];

        let xabs = self.abs();
        let (p, q) = if xabs >= Self::from_uint(8u8)
        {
            (&QR8, &QS8)
        }
        else if xabs >= Self::from(4.5454)
        {
            (&QR5, &QS5)
        }
        else if xabs >= Self::from(2.857)
        {
            (&QR3, &QS3)
        }
        else
        {
            (&QR2, &QS2)
        };

        let one = Self::one();
        let z = (self*self).recip();
        let mut b = false;

        // TODO: Make cleaner
        let mut r = util::polynomial(p, z, SIGN_BIT, false);
        let mut s = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && r.is_nan()
        {
            b = !b;
            r = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && s.is_nan()
        {
            b = !b;
            s = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false) - one;
        }

        let rs = r/s;
        let ofs = Self::from(0.125);

        let (diff, s) = rs.add_with_sign_extra_sign(b, ofs, true);
        let y = diff/self;
        (y, s)
    }

    fn bessel0_common(self, y0: bool) -> (Self, bool)
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01;

        /*
         * j0(x) = sqrt(2/(pi*x))*(p0(x)*cos(x-pi/4)-q0(x)*sin(x-pi/4))
         * y0(x) = sqrt(2/(pi*x))*(p0(x)*sin(x-pi/4)+q0(x)*cos(x-pi/4))
         *
         * sin(x-pi/4) = (sin(x) - cos(x))/sqrt(2)
         * cos(x-pi/4) = (sin(x) + cos(x))/sqrt(2)
         * sin(x) +- cos(x) = -cos(2x)/(sin(x) -+ cos(x))
         */
        let (s, s_s) = self.sin_extra_sign();
        let (c, mut c_s) = self.cos_extra_sign();
        c_s ^= y0;
        
        let (mut cc, mut cc_s) = s.add_with_sign_extra_sign(s_s, c, c_s);

        let two = Self::from_uint(2u8);
        let xabs = self.abs();
        /* avoid overflow in 2*x, big ulp error when x>=0x1p1023 */
        if xabs < Self::max_value()/two
        {
            let (mut ss, mut ss_s) = s.add_with_sign_extra_sign(s_s, c, !c_s);
            let (z, mut z_s) = (two*self).cos_extra_sign();
            z_s = !z_s;

            if (s*c).is_sign_negative() ^ s_s ^ c_s
            {
                cc_s = z_s^ss_s;
                cc = z / ss;
            }
            else
            {
                ss_s = z_s^cc_s;
                ss = z / cc;
            }

            if xabs < Self::from((EXP_BASE as f64).powf(17.0/127.0*Self::exp_bias().to_f64().unwrap()))
            {
                ss_s ^= y0;
                
                let (mut p, mut p_s) = self.bessel0_p();
                p *= cc;
                p_s ^= cc_s;

                let (mut q, mut q_s) = self.bessel0_q();
                q *= ss;
                q_s ^= ss_s;

                (cc, cc_s) = p.add_with_sign_extra_sign(p_s, q, !q_s);
            }
        }
        
        (Self::from(INVSQRTPI)*cc / self.sqrt(), cc_s)
    }

    fn j0_extra_sign(mut self) -> (Self, bool)
    {
        match self.classify()
        {
            // TODO: Add .squared()
            FpCategory::Nan => ((self*self).recip(), false),
            FpCategory::Infinite => (Self::zero(), false),
            FpCategory::Zero => (Self::one(), false),
            FpCategory::Normal | FpCategory::Subnormal => {
                self = self.abs();
                let two = Self::from_uint(2u8);
                if self >= two
                {
                    /* |x| >= 2 */
                    /* large ulp error near zeros: 2.4, 5.52, 8.6537,.. */
                    return self.bessel0_common(false)
                }
                let one = Self::one();

                /* 1 - x*x/4 + x*x*R(x^2)/S(x^2) */
                if self >= Self::from(0.0001220703125)
                {
                    /* R0/S0 on [0, 2.00] */
                    const R0: [f64; 4] = [
                        1.56249999999999947958e-02, /* 0x3F8FFFFF, 0xFFFFFFFD */
                        -1.89979294238854721751e-04, /* 0xBF28E6A5, 0xB61AC6E9 */
                        1.82954049532700665670e-06, /* 0x3EBEB1D1, 0x0C503919 */
                        -4.61832688532103189199e-09 /* 0xBE33D5E7, 0x73D63FCE */
                    ];
                    const S0: [f64; 4] = [
                        1.56191029464890010492e-02, /* 0x3F8FFCE8, 0x82C8C2A4 */
                        1.16926784663337450260e-04, /* 0x3F1EA6D2, 0xDD57DBF4 */
                        5.13546550207318111446e-07, /* 0x3EA13B54, 0xCE84D5A9 */
                        1.16614003333790000205e-09 /* 0x3E1408BC, 0xF4745D8F */
                    ];

                    /* |x| >= 2**-13 */
                    /* up to 4ulp error close to 2 */
                    // TODO: add .squared()
                    let z = self*self;
                    let mut b = false;
                    let mut numer = z*util::polynomial(&R0, z, SIGN_BIT, false);
                    let mut denom = one + z*util::polynomial(&S0, z, SIGN_BIT, false);
                    if !SIGN_BIT && numer.is_nan()
                    {
                        b = !b;
                        numer = z*util::polynomial(&R0.map(Neg::neg), z, SIGN_BIT, false)
                    }
                    if !SIGN_BIT && denom.is_nan()
                    {
                        denom = z*util::polynomial(&S0.map(Neg::neg), z, SIGN_BIT, false);
                        if denom > one
                        {
                            b = !b;
                            denom -= one
                        }
                        else
                        {
                            denom = one - denom
                        }
                    }

                    let rs = numer/denom;
                    let zrs = z*rs;
                    let d = (one + self/two)*(one - self/two);

                    return d.add_with_sign_extra_sign(false, zrs, b)
                }

                /* 1 - x*x/4 */
                /* prevent underflow */
                /* inexact should be raised when x!=0, this is not done correctly */
                if self >= Self::from((EXP_BASE as f64).powf(-13.0/127.0*Self::exp_bias().to_f64().unwrap()))
                {
                    /* |x| >= 2**-127 */
                    self = Self::from(0.25)*self*self;
                }

                one.sub_extra_sign(self)
            }
        }
    }
    
    /// Bessel function of the first kind with α = 0.
    pub fn j0(self) -> Self
    {
        let (y, s) = self.j0_extra_sign();
        y.xor_sign(s)
    }

    fn y0_extra_sign(self) -> (Self, bool)
    {
        const TPI: f64 = 6.36619772367581382433e-01; /* 0x3FE45F30, 0x6DC9C883 */
        
        const U0: [f64; 7] = [
            -7.38042951086872317523e-02, /* 0xBFB2E4D6, 0x99CBD01F */
            1.76666452509181115538e-01, /* 0x3FC69D01, 0x9DE9E3FC */
            -1.38185671945596898896e-02, /* 0xBF8C4CE8, 0xB16CFA97 */
            3.47453432093683650238e-04, /* 0x3F36C54D, 0x20B29B6B */
            -3.81407053724364161125e-06, /* 0xBECFFEA7, 0x73D25CAD */
            1.95590137035022920206e-08, /* 0x3E550057, 0x3B4EABD4 */
            -3.98205194132103398453e-11 /* 0xBDC5E43D, 0x693FB3C8 */
        ];
        const V0: [f64; 4] = [
            1.27304834834123699328e-02, /* 0x3F8A1270, 0x91C9C71A */
            7.60068627350353253702e-05, /* 0x3F13ECBB, 0xF578C6C1 */
            2.59150851840457805467e-07, /* 0x3E91642D, 0x7FF202FD */
            4.41110311332675467403e-10 /* 0x3DFE5018, 0x3BD6D9EF */
        ];

        /* y0(nan)=nan, y0(<0)=nan, y0(0)=-inf, y0(inf)=0 */
        match self.classify()
        {
            FpCategory::Zero => (Self::infinity(), true),
            FpCategory::Nan => (self, false),
            FpCategory::Infinite => if self.is_sign_negative()
            {
                (Self::nan(), false)
            }
            else
            {
                (Self::zero(), false)
            },
            FpCategory::Normal | FpCategory::Subnormal => if self.is_sign_negative()
            {
                (Self::nan(), false)
            }
            else if self >= Self::from_uint(2u8)
            {
                /* x >= 2 */
                /* large ulp errors near zeros: 3.958, 7.086,.. */
                self.bessel0_common(true)
            }
            else if self >= Self::from((EXP_BASE as f64).powf(-1.625*EXP_SIZE as f64))
            {
                /* U(x^2)/V(x^2) + (2/pi)*j0(x)*log(x) */
                /* large ulp error near the first zero, x ~= 0.89 */
                let z = self*self;
    
                // TODO: Make cleaner
                let mut b = false;
                let mut u = util::polynomial(&U0, z, SIGN_BIT, false);
                let mut v = Self::one() + z*util::polynomial(&V0, z, SIGN_BIT, false);
                if !SIGN_BIT && u.is_nan()
                {
                    b = !b;
                    u = util::polynomial(&U0.map(Neg::neg), z, SIGN_BIT, false);
                }
                let one = Self::one();
                if !SIGN_BIT && v.is_nan()
                {
                    v = z*util::polynomial(&V0.map(Neg::neg), z, SIGN_BIT, false);
                    if v > one
                    {
                        b = !b;
                        v -= one;
                    }
                    else
                    {
                        v = one - v;
                    }
                }
    
                let uv = u/v;
                let (xj0, xj0_s) = self.j0_extra_sign();
                let mut s_l = !SIGN_BIT && self < one;
                let ll = if s_l
                {
                    self.recip()
                }
                else
                {
                    self
                };
                s_l ^= xj0_s;
                let l = Self::from(TPI)*(xj0*ll.ln());

                uv.add_with_sign_extra_sign(b, l, s_l)
            }
            else
            {
                let u00 = Self::from(U0[0]);
                let s_l = !SIGN_BIT;
                let ll = if s_l
                {
                    self.recip()
                }
                else
                {
                    self
                };
                let l = Self::from(TPI)*ll.ln();
                u00.add_with_sign_extra_sign(false, l, s_l)
            }
        }
    }

    /// Bessel function of the second kind with α = 0.
    pub fn y0(self) -> Self
    {
        let (y, s) = self.y0_extra_sign();
        y.xor_sign(s)
    }
    
    fn bessel1_p(self) -> (Self, bool)
    {
        /* For x >= 8, the asymptotic expansions of pone is
        *      1 + 15/128 s^2 - 4725/2^15 s^4 - ...,   where s = 1/x.
        * We approximate pone by
        *      pone(x) = 1 + (R/S)
        * where  R = pr0 + pr1*s^2 + pr2*s^4 + ... + pr5*s^10
        *        S = 1 + ps0*s^2 + ... + ps4*s^10
        * and
        *      | pone(x)-1-R/S | <= 2  ** ( -60.06)
        */

        const PR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00, /* 0x00000000, 0x00000000 */
            1.17187499999988647970e-01, /* 0x3FBDFFFF, 0xFFFFFCCE */
            1.32394806593073575129e+01, /* 0x402A7A9D, 0x357F7FCE */
            4.12051854307378562225e+02, /* 0x4079C0D4, 0x652EA590 */
            3.87474538913960532227e+03, /* 0x40AE457D, 0xA3A532CC */
            7.91447954031891731574e+03, /* 0x40BEEA7A, 0xC32782DD */
        ];
        const PS8: [f64; 5] = [
            1.14207370375678408436e+02, /* 0x405C8D45, 0x8E656CAC */
            3.65093083420853463394e+03, /* 0x40AC85DC, 0x964D274F */
            3.69562060269033463555e+04, /* 0x40E20B86, 0x97C5BB7F */
            9.76027935934950801311e+04, /* 0x40F7D42C, 0xB28F17BB */
            3.08042720627888811578e+04, /* 0x40DE1511, 0x697A0B2D */
        ];

        const PR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            1.31990519556243522749e-11, /* 0x3DAD0667, 0xDAE1CA7D */
            1.17187493190614097638e-01, /* 0x3FBDFFFF, 0xE2C10043 */
            6.80275127868432871736e+00, /* 0x401B3604, 0x6E6315E3 */
            1.08308182990189109773e+02, /* 0x405B13B9, 0x452602ED */
            5.17636139533199752805e+02, /* 0x40802D16, 0xD052D649 */
            5.28715201363337541807e+02, /* 0x408085B8, 0xBB7E0CB7 */
        ];
        const PS5: [f64; 5] = [
            5.92805987221131331921e+01, /* 0x404DA3EA, 0xA8AF633D */
            9.91401418733614377743e+02, /* 0x408EFB36, 0x1B066701 */
            5.35326695291487976647e+03, /* 0x40B4E944, 0x5706B6FB */
            7.84469031749551231769e+03, /* 0x40BEA4B0, 0xB8A5BB15 */
            1.50404688810361062679e+03, /* 0x40978030, 0x036F5E51 */
        ];

        const PR3: [f64; 6] = [
            3.02503916137373618024e-09, /* 0x3E29FC21, 0xA7AD9EDD */
            1.17186865567253592491e-01, /* 0x3FBDFFF5, 0x5B21D17B */
            3.93297750033315640650e+00, /* 0x400F76BC, 0xE85EAD8A */
            3.51194035591636932736e+01, /* 0x40418F48, 0x9DA6D129 */
            9.10550110750781271918e+01, /* 0x4056C385, 0x4D2C1837 */
            4.85590685197364919645e+01, /* 0x4048478F, 0x8EA83EE5 */
        ];
        const PS3: [f64; 5] = [
            3.47913095001251519989e+01, /* 0x40416549, 0xA134069C */
            3.36762458747825746741e+02, /* 0x40750C33, 0x07F1A75F */
            1.04687139975775130551e+03, /* 0x40905B7C, 0x5037D523 */
            8.90811346398256432622e+02, /* 0x408BD67D, 0xA32E31E9 */
            1.03787932439639277504e+02, /* 0x4059F26D, 0x7C2EED53 */
        ];

        const PR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            1.07710830106873743082e-07, /* 0x3E7CE9D4, 0xF65544F4 */
            1.17176219462683348094e-01, /* 0x3FBDFF42, 0xBE760D83 */
            2.36851496667608785174e+00, /* 0x4002F2B7, 0xF98FAEC0 */
            1.22426109148261232917e+01, /* 0x40287C37, 0x7F71A964 */
            1.76939711271687727390e+01, /* 0x4031B1A8, 0x177F8EE2 */
            5.07352312588818499250e+00, /* 0x40144B49, 0xA574C1FE */
        ];
        const PS2: [f64; 5] = [
            2.14364859363821409488e+01, /* 0x40356FBD, 0x8AD5ECDC */
            1.25290227168402751090e+02, /* 0x405F5293, 0x14F92CD5 */
            2.32276469057162813669e+02, /* 0x406D08D8, 0xD5A2DBD9 */
            1.17679373287147100768e+02, /* 0x405D6B7A, 0xDA1884A9 */
            8.36463893371618283368e+00, /* 0x4020BAB1, 0xF44E5192 */
        ];
        
        let xabs = self.abs();
        let (p, q) = if xabs >= Self::from_uint(8u8)
        {
            (&PR8, &PS8)
        }
        else if xabs >= Self::from(4.5454)
        {
            (&PR5, &PS5)
        }
        else if xabs >= Self::from(2.857)
        {
            (&PR3, &PS3)
        }
        else
        {
            (&PR2, &PS2)
        };
        let one = Self::one();
        let z = self.squared().recip();

        let mut y_s = false;
        let mut numer = util::polynomial(p, z, SIGN_BIT, false);
        let mut denom = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && numer.is_nan()
        {
            y_s = !y_s;
            numer = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && denom.is_nan()
        {
            denom = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false);
            if denom > one
            {
                y_s = !y_s;
                denom -= one
            }
            else
            {
                denom = one - denom
            }
        }
        
        let y = numer/denom;
        one.add_with_sign_extra_sign(false, y, y_s)
    }

    fn bessel1_q(self) -> (Self, bool)
    {
        /* For x >= 8, the asymptotic expansions of qone is
        *      3/8 s - 105/1024 s^3 - ..., where s = 1/x.
        * We approximate pone by
        *      qone(x) = s*(0.375 + (R/S))
        * where  R = qr1*s^2 + qr2*s^4 + ... + qr5*s^10
        *        S = 1 + qs1*s^2 + ... + qs6*s^12
        * and
        *      | qone(x)/s -0.375-R/S | <= 2  ** ( -61.13)
        */

        const QR8: [f64; 6] = [
            /* for x in [inf, 8]=1/[0,0.125] */
            0.00000000000000000000e+00,  /* 0x00000000, 0x00000000 */
            -1.02539062499992714161e-01, /* 0xBFBA3FFF, 0xFFFFFDF3 */
            -1.62717534544589987888e+01, /* 0xC0304591, 0xA26779F7 */
            -7.59601722513950107896e+02, /* 0xC087BCD0, 0x53E4B576 */
            -1.18498066702429587167e+04, /* 0xC0C724E7, 0x40F87415 */
            -4.84385124285750353010e+04, /* 0xC0E7A6D0, 0x65D09C6A */
        ];
        const QS8: [f64; 6] = [
            1.61395369700722909556e+02,  /* 0x40642CA6, 0xDE5BCDE5 */
            7.82538599923348465381e+03,  /* 0x40BE9162, 0xD0D88419 */
            1.33875336287249578163e+05,  /* 0x4100579A, 0xB0B75E98 */
            7.19657723683240939863e+05,  /* 0x4125F653, 0x72869C19 */
            6.66601232617776375264e+05,  /* 0x412457D2, 0x7719AD5C */
            -2.94490264303834643215e+05, /* 0xC111F969, 0x0EA5AA18 */
        ];

        const QR5: [f64; 6] = [
            /* for x in [8,4.5454]=1/[0.125,0.22001] */
            -2.08979931141764104297e-11, /* 0xBDB6FA43, 0x1AA1A098 */
            -1.02539050241375426231e-01, /* 0xBFBA3FFF, 0xCB597FEF */
            -8.05644828123936029840e+00, /* 0xC0201CE6, 0xCA03AD4B */
            -1.83669607474888380239e+02, /* 0xC066F56D, 0x6CA7B9B0 */
            -1.37319376065508163265e+03, /* 0xC09574C6, 0x6931734F */
            -2.61244440453215656817e+03, /* 0xC0A468E3, 0x88FDA79D */
        ];
        const QS5: [f64; 6] = [
            8.12765501384335777857e+01,  /* 0x405451B2, 0xFF5A11B2 */
            1.99179873460485964642e+03,  /* 0x409F1F31, 0xE77BF839 */
            1.74684851924908907677e+04,  /* 0x40D10F1F, 0x0D64CE29 */
            4.98514270910352279316e+04,  /* 0x40E8576D, 0xAABAD197 */
            2.79480751638918118260e+04,  /* 0x40DB4B04, 0xCF7C364B */
            -4.71918354795128470869e+03, /* 0xC0B26F2E, 0xFCFFA004 */
        ];

        const QR3: [f64; 6] = [
            -5.07831226461766561369e-09, /* 0xBE35CFA9, 0xD38FC84F */
            -1.02537829820837089745e-01, /* 0xBFBA3FEB, 0x51AEED54 */
            -4.61011581139473403113e+00, /* 0xC01270C2, 0x3302D9FF */
            -5.78472216562783643212e+01, /* 0xC04CEC71, 0xC25D16DA */
            -2.28244540737631695038e+02, /* 0xC06C87D3, 0x4718D55F */
            -2.19210128478909325622e+02, /* 0xC06B66B9, 0x5F5C1BF6 */
        ];
        const QS3: [f64; 6] = [
            4.76651550323729509273e+01,  /* 0x4047D523, 0xCCD367E4 */
            6.73865112676699709482e+02,  /* 0x40850EEB, 0xC031EE3E */
            3.38015286679526343505e+03,  /* 0x40AA684E, 0x448E7C9A */
            5.54772909720722782367e+03,  /* 0x40B5ABBA, 0xA61D54A6 */
            1.90311919338810798763e+03,  /* 0x409DBC7A, 0x0DD4DF4B */
            -1.35201191444307340817e+02, /* 0xC060E670, 0x290A311F */
        ];

        const QR2: [f64; 6] = [
            /* for x in [2.8570,2]=1/[0.3499,0.5] */
            -1.78381727510958865572e-07, /* 0xBE87F126, 0x44C626D2 */
            -1.02517042607985553460e-01, /* 0xBFBA3E8E, 0x9148B010 */
            -2.75220568278187460720e+00, /* 0xC0060484, 0x69BB4EDA */
            -1.96636162643703720221e+01, /* 0xC033A9E2, 0xC168907F */
            -4.23253133372830490089e+01, /* 0xC04529A3, 0xDE104AAA */
            -2.13719211703704061733e+01, /* 0xC0355F36, 0x39CF6E52 */
        ];
        const QS2: [f64; 6] = [
            2.95333629060523854548e+01,  /* 0x403D888A, 0x78AE64FF */
            2.52981549982190529136e+02,  /* 0x406F9F68, 0xDB821CBA */
            7.57502834868645436472e+02,  /* 0x4087AC05, 0xCE49A0F7 */
            7.39393205320467245656e+02,  /* 0x40871B25, 0x48D4C029 */
            1.55949003336666123687e+02,  /* 0x40637E5E, 0x3C3ED8D4 */
            -4.95949898822628210127e+00, /* 0xC013D686, 0xE71BE86B */
        ];
        
        let xabs = self.abs();
        let (p, q) = if xabs >= Self::from_uint(8u8)
        {
            (&QR8, &QS8)
        }
        else if xabs >= Self::from(4.5454)
        {
            (&QR5, &QS5)
        }
        else if xabs >= Self::from(2.857)
        {
            (&QR3, &QS3)
        }
        else
        {
            (&QR2, &QS2)
        };

        let one = Self::one();
        let z = xabs.squared().recip();
        let mut y_s = false;
        let mut numer = util::polynomial(p, z, SIGN_BIT, false);
        let mut denom = one + z*util::polynomial(q, z, SIGN_BIT, false);
        if !SIGN_BIT && numer.is_nan()
        {
            y_s = !y_s;
            numer = util::polynomial(&p.map(Neg::neg), z, SIGN_BIT, false);
        }
        if !SIGN_BIT && denom.is_nan()
        {
            y_s = !y_s;
            denom = z*util::polynomial(&q.map(Neg::neg), z, SIGN_BIT, false) - one;
        }

        let y = numer/denom;
        let ofs = Self::from(0.375);
        let (sum, s) = ofs.add_with_sign_extra_sign(false, y, y_s);
        (sum/self, s)
    }
    
    fn bessel1_common(self, y1: bool, sign: bool) -> (Self, bool)
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01;

        /*
        * j1(x) = sqrt(2/(pi*x))*(p1(x)*cos(x-3pi/4)-q1(x)*sin(x-3pi/4))
        * y1(x) = sqrt(2/(pi*x))*(p1(x)*sin(x-3pi/4)+q1(x)*cos(x-3pi/4))
        *
        * sin(x-3pi/4) = -(sin(x) + cos(x))/sqrt(2)
        * cos(x-3pi/4) = (sin(x) - cos(x))/sqrt(2)
        * sin(x) +- cos(x) = -cos(2x)/(sin(x) -+ cos(x))
        */
        let (s, mut s_s) = self.sin_extra_sign();
        s_s ^= y1;
        let (c, c_s) = self.cos_extra_sign();
        
        let (mut cc, mut cc_s) = s.add_with_sign_extra_sign(s_s, c, !c_s);

        let two = Self::from(2u8);
        let xabs = self.abs();
        if xabs < Self::max_value()/two
        {
            /* avoid overflow in 2*x */
            let (mut ss, mut ss_s) = s.add_with_sign_extra_sign(!s_s, c, !c_s);
            
            let (z, z_s) = (self + self).cos_extra_sign();
            
            if (s*c).is_sign_negative() ^ s_s ^ c_s
            {
                cc_s = z_s^ss_s;
                cc = z / ss;
            }
            else
            {
                ss_s = z_s^cc_s;
                ss = z / cc;
            }

            if xabs < Self::from((EXP_BASE as f64).powf(17.0/127.0*Self::exp_bias().to_f64().unwrap()))
            {
                ss_s ^= y1;
                
                let (mut p, mut p_s) = self.bessel1_p();
                p *= cc;
                p_s ^= cc_s;

                let (mut q, mut q_s) = self.bessel1_q();
                q *= ss;
                q_s ^= ss_s;

                (cc, cc_s) = p.add_with_sign_extra_sign(p_s, q, !q_s);
            }
        }
        cc_s ^= sign;

        (Self::from(INVSQRTPI)*cc/xabs.sqrt(), cc_s)
    }

    fn j1_extra_sign(self) -> (Self, bool)
    {
        match self.classify()
        {
            FpCategory::Infinite => (Self::zero(), false),
            FpCategory::Nan => (self.squared().recip(), false),
            FpCategory::Zero => (self, false),
            FpCategory::Normal | FpCategory::Subnormal => {
                let sign = self.is_sign_negative();
                let xabs = self.abs();
                let two = Self::from_uint(2u8);

                if xabs >= two
                {
                    /* |x| >= 2 */
                    return xabs.bessel1_common(false, sign)
                }
                let half = two.recip();
                let mut z;
                let mut b = false;
                if xabs >= Self::from((EXP_BASE as f64).powf(-13.0/127.0*Self::exp_bias().to_f64().unwrap()))
                {
                    /* R0/S0 on [0,2] */
                    const R0: [f64; 4] = [
                        -6.25000000000000000000e-02, /* 0xBFB00000, 0x00000000 */
                        1.40705666955189706048e-03, /* 0x3F570D9F, 0x98472C61 */
                        -1.59955631084035597520e-05, /* 0xBEF0C5C6, 0xBA169668 */
                        4.96727999609584448412e-08 /* 0x3E6AAAFA, 0x46CA0BD9 */
                    ];
                    const S0: [f64; 5] = [
                        1.91537599538363460805e-02, /* 0x3F939D0B, 0x12637E53 */
                        1.85946785588630915560e-04, /* 0x3F285F56, 0xB9CDF664 */
                        1.17718464042623683263e-06, /* 0x3EB3BFF8, 0x333F8498 */
                        5.04636257076217042715e-09, /* 0x3E35AC88, 0xC97DFF2C */
                        1.23542274426137913908e-11 /* 0x3DAB2ACF, 0xCFB97ED8 */
                    ];

                    /* |x| >= 2**-127 */
                    z = xabs.squared();

                    let mut r = z*util::polynomial(&R0, z, SIGN_BIT, false);
                    let mut s = Self::one() + z*util::polynomial(&S0, z, SIGN_BIT, false);
                    if !SIGN_BIT && r.is_nan()
                    {
                        b = !b;
                        r = z*util::polynomial(&R0.map(Neg::neg), z, SIGN_BIT, false);
                    }
                    if !SIGN_BIT && s.is_nan()
                    {
                        let one = Self::one();
                        s = z*util::polynomial(&S0.map(Neg::neg), z, SIGN_BIT, false);
                        if s > one
                        {
                            b = !b;
                            s -= one
                        }
                        else
                        {
                            s = one - s
                        }
                    }

                    z = r/s;
                }
                else
                {
                    /* avoid underflow, raise inexact if x!=0 */
                    z = self;
                }

                let (sum, y_s) = half.add_with_sign_extra_sign(false, z, b);

                (sum*self, y_s)
            }
        }
    }

    /// Bessel function of the first kind with α = 1.
    pub fn j1(self) -> Self
    {
        let (y, s) = self.j1_extra_sign();
        y.xor_sign(s)
    }

    fn y1_extra_sign(self) -> (Self, bool)
    {
        const TPI: f64 = 6.36619772367581382433e-01; /* 0x3FE45F30, 0x6DC9C883 */
        
        const U0: [f64; 5] = [
            -1.96057090646238940668e-01, /* 0xBFC91866, 0x143CBC8A */
            5.04438716639811282616e-02,  /* 0x3FA9D3C7, 0x76292CD1 */
            -1.91256895875763547298e-03, /* 0xBF5F55E5, 0x4844F50F */
            2.35252600561610495928e-05,  /* 0x3EF8AB03, 0x8FA6B88E */
            -9.19099158039878874504e-08, /* 0xBE78AC00, 0x569105B8 */
        ];
        const V0: [f64; 5] = [
            1.99167318236649903973e-02, /* 0x3F94650D, 0x3F4DA9F0 */
            2.02552581025135171496e-04, /* 0x3F2A8C89, 0x6C257764 */
            1.35608801097516229404e-06, /* 0x3EB6C05A, 0x894E8CA6 */
            6.22741452364621501295e-09, /* 0x3E3ABF1D, 0x5BA69A86 */
            1.66559246207992079114e-11, /* 0x3DB25039, 0xDACA772A */
        ];

        /* y1(nan)=nan, y1(<0)=nan, y1(0)=-inf, y1(inf)=0 */
        match self.classify()
        {
            FpCategory::Nan => (self, false),
            FpCategory::Zero => (Self::infinity(), true),
            _ if self.is_sign_negative() => (Self::nan(), false),
            FpCategory::Infinite => (Self::zero(), false),
            FpCategory::Normal | FpCategory::Subnormal => {
                if self >= Self::from_uint(2u8)
                {
                    /* x >= 2 */
                    return self.bessel1_common(true, false)
                }
                if self < Self::from((EXP_BASE as f64).powf(-3.125*EXP_SIZE as f64))
                {
                    return (Self::from(TPI)/self, true)
                }
        
                let one = Self::one();
                let z = self.squared();
        
                let mut uv_s = false;
                let mut numer = util::polynomial(&U0, z, SIGN_BIT, false);
                let mut denom = Self::one() + z*util::polynomial(&V0, z, SIGN_BIT, false);
                if !SIGN_BIT && numer.is_nan()
                {
                    uv_s = !uv_s;
                    numer = util::polynomial(&U0.map(Neg::neg), z, SIGN_BIT, false)
                }
                if !SIGN_BIT && denom.is_nan()
                {
                    denom = z*util::polynomial(&V0.map(Neg::neg), z, SIGN_BIT, false);
                    if denom > one
                    {
                        uv_s = !uv_s;
                        denom -= one
                    }
                    else
                    {
                        denom = one - denom
                    }
                }
                
                let uv = numer/denom;
                let xuv = self*uv;
        
                let (xj1, xj1_s) = self.j1_extra_sign();
                let xln_s = self < one;
                let xinv = self.recip();
                let xln_arg = if xln_s
                {
                    xinv
                }
                else
                {
                    self
                };

                let xln = xln_arg.ln();
                let xj1mxln = xj1*xln;
                let xj1mxln_s = xj1_s ^ xln_s;

                let (diff, diff_s) = xj1mxln.add_with_sign_extra_sign(xj1mxln_s, xinv, true);
                xuv.add_with_sign_extra_sign(uv_s, Self::from(TPI)*diff, diff_s)
            }
        }
    }
    
    /// Bessel function of the second kind with α = 1.
    pub fn y1(self) -> Self
    {
        let (y, s) = self.y1_extra_sign();
        y.xor_sign(s)
    }

    /// Bessel function of the first kind with α = `n`.
    pub fn jn(mut self, n: i32) -> Self
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01; /* 0x3FE20DD7, 0x50429B6D */

        let class = self.classify();
        match class
        {
            FpCategory::Nan => self,
            FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal | FpCategory::Zero => {
                let mut sign = self.is_sign_negative();
        
                /* J(-n,x) = (-1)^n * J(n, x), J(n, -x) = (-1)^n * J(n, x)
                * Thus, J(-n,x) = J(n,-x)
                */
                /* nm1 = |n|-1 is used instead of |n| to handle n==INT_MIN */
                if n == 0
                {
                    return self.j0()
                }
                let nm1 = if n < 0
                {
                    self = -self;
                    sign = !sign;
                    -(n + 1)
                }
                else
                {
                    n - 1
                };
                if nm1 == 0
                {
                    return self.j1();
                }
        
                sign &= (n & 1) != 0; /* even n: 0, odd n: signbit(x) */
                self = self.abs();
        
                let one = Self::one();
                let two = one + one;
                let h = two/self;
        
                let (b, b_s) = if matches!(class, FpCategory::Zero | FpCategory::Infinite)
                {
                    /* if x is 0 or inf */
                    (Self::zero(), false)
                }
                else if Self::from_int(nm1) < self
                {
                    /* Safe to use J(n+1,x)=2n/x *J(n,x)-J(n-1,x) */
                    if self >= Self::from(2.1359870359209100823950217061696e96)
                    {
                        /* x > 2**302 */
                        /* (x >> n**2)
                        *      Jn(x) = cos(x-(2n+1)*pi/4)*sqrt(2/x*pi)
                        *      Yn(x) = sin(x-(2n+1)*pi/4)*sqrt(2/x*pi)
                        *      Let s=sin(x), c=cos(x),
                        *          xn=x-(2n+1)*pi/4, sqt2 = sqrt(2),then
                        *
                        *             n    sin(xn)*sqt2    cos(xn)*sqt2
                        *          ----------------------------------
                        *             0     s-c             c+s
                        *             1    -s-c            -c+s
                        *             2    -s+c            -c-s
                        *             3     s+c             c-s
                        */
                        let (s, s_s) = self.sin_extra_sign();
                        let (c, c_s) = self.cos_extra_sign();
                        let (n0, n1) = (nm1 & 0b01 == 0, nm1 & 0b10 == 0);
                        let s1 = c_s ^ n1;
                        let s2 = s_s ^ n0 ^ n1;
                        let (sum, b_s) = c.add_with_sign_extra_sign(s1, s, s2);
                        (Self::from(INVSQRTPI)*sum/self.sqrt(), b_s)
                    }
                    else
                    {
                        let (mut a, mut a_s) = self.j0_extra_sign();
                        let (mut b, mut b_s) = self.j1_extra_sign();
                        for i in 1..=nm1
                        {
                            (b, b_s) = (b*Self::from_int(i)*h).add_with_sign_extra_sign(
                                b_s,
                                core::mem::replace(&mut a, b),
                                !core::mem::replace(&mut a_s, b_s)
                            ); /* avoid underflow */
                        }
                        (b, b_s)
                    }
                }
                else if self < Self::from(0.00000000186264514923095703125)
                {
                    /* x < 2**-29 */
                    /* x is tiny, return the first Taylor expansion of J(n,x)
                    * J(n,x) = 1/n!*(x/2)^n  - ...
                    */
                    if nm1 > 32
                    {
                        /* underflow */
                        (Self::zero(), false)
                    }
                    else
                    {
                        let x_half = h.recip();
                        let (mut b, b_s) = (x_half, false);
                        let mut a = one;
                        for i in 2..=nm1 + 1
                        {
                            a *= Self::from_int(i); /* a = n! */
                            b *= x_half; /* b = (x/2)^n */
                        }
                        b /= a;
                        (b, b_s)
                    }
                }
                else
                {
                    /* use backward recurrence */
                    /*                      x      x^2      x^2
                    *  J(n,x)/J(n-1,x) =  ----   ------   ------   .....
                    *                      2n  - 2(n+1) - 2(n+2)
                    *
                    *                      1      1        1
                    *  (for large x)   =  ----  ------   ------   .....
                    *                      2n   2(n+1)   2(n+2)
                    *                      -- - ------ - ------ -
                    *                       x     x         x
                    *
                    * Let w = 2n/x and h=2/x, then the above quotient
                    * is equal to the continued fraction:
                    *                  1
                    *      = -----------------------
                    *                     1
                    *         w - -----------------
                    *                        1
                    *              w+h - ---------
                    *                     w+2h - ...
                    *
                    * To determine how many terms needed, let
                    * Q(0) = w, Q(1) = w(w+h) - 1,
                    * Q(k) = (w+k*h)*Q(k-1) - Q(k-2),
                    * When Q(k) > 1e4      good for single
                    * When Q(k) > 1e9      good for double
                    * When Q(k) > 1e17     good for quadruple
                    */
                    /* determine k */
        
                    let nf = Self::from_int(nm1) + one;
                    let w = h*nf;
                    let mut z = w + h;
                    let wz = w*z;
        
                    let (mut q0, mut q0_s) = (w, false);
                    let (mut q1, mut q1_s) = wz.sub_extra_sign(one);
        
                    let mut k = 1;
                    let q1_min = Self::from((EXP_BASE as f64).powf(13.287712379549449391481277717958/8.0*EXP_SIZE as f64));
                    while q1 < q1_min
                    {
                        k += 1;
                        z += h;
                        (q1, q1_s) = (z*q1).add_with_sign_extra_sign(
                            q1_s,
                            core::mem::replace(&mut q0, q1),
                            !core::mem::replace(&mut q0_s, q1_s)
                        );
                    }
                    let mut t = Self::zero();
                    let mut i = k;
                    while i >= 0
                    {
                        t = ((Self::from_int(i) + nf)*h - t).recip();
                        assert!(!t.is_nan());
                        i -= 1;
                    }
                    let (mut a, mut a_s) = (t, false);
                    let (mut b, mut b_s) = (one, false);

                    /*  estimate log((2/x)^n*n!) = n*log(2/x)+n*ln(n)
                    *  Hence, if n*(log(2n/x)) > ...
                    *  single 8.8722839355e+01
                    *  double 7.09782712893383973096e+02
                    *  long double 1.1356523406294143949491931077970765006170e+04
                    *  then recurrent value may overflow and the result is
                    *  likely underflow to zero
                    */
                    let wabs = w.abs();
                    let tmp = nf*wabs.ln();
                    if tmp < Self::from(0.69860503429133858267716535433071*Self::exp_bias().to_f64().unwrap()) || wabs < one
                    {
                        for i in (1..=nm1).rev()
                        {
                            (b, b_s) = (b*Self::from_int(i)*h).add_with_sign_extra_sign(
                                b_s,
                                core::mem::replace(&mut a, b),
                                !core::mem::replace(&mut a_s, b_s)
                            );
                        }
                    }
                    else
                    {
                        let x1p500 = Self::from((EXP_BASE as f64).powf(0.47244094488188976377952755905512*Self::exp_bias().to_f64().unwrap())); // 0x1p500 == 2^500
        
                        for i in (1..=nm1).rev()
                        {
                            (b, b_s) = (b*Self::from_int(i)*h).add_with_sign_extra_sign(
                                b_s,
                                core::mem::replace(&mut a, b),
                                !core::mem::replace(&mut a_s, b_s)
                            );
        
                            /* scale b to avoid spurious overflow */
                            if b.abs() > x1p500
                            {
                                a /= b;
                                t /= b;
                                b = one;
                            }
                        }
                    }
                    let (z, z_s) = self.j0_extra_sign();
                    let (w, w_s) = self.j1_extra_sign();
                    if z.abs() >= w.abs()
                    {
                        b_s ^= z_s;
                        b = t * z / b;
                    }
                    else
                    {
                        b_s = w_s^a_s;
                        b = t * w / a;
                    }

                    (b, b_s)
                };

                b.xor_sign(sign^b_s)
            }
        }
    }

    /// Bessel function of the second kind with α = `n`.
    pub fn yn(self, n: i32) -> Self
    {
        const INVSQRTPI: f64 = 5.64189583547756279280e-01; /* 0x3FE20DD7, 0x50429B6D */

        match self.classify()
        {
            FpCategory::Nan => self,
            FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal if self.is_sign_negative() => Self::nan(),
            FpCategory::Infinite => Self::zero(),
            FpCategory::Zero | FpCategory::Normal | FpCategory::Subnormal => {
                if n == 0
                {
                    return self.y0();
                }
                let (nm1, sign) = if n < 0
                {
                    (-(n + 1), (n & 1) != 0)
                }
                else
                {
                    (n - 1, false)
                };
                if nm1 == 0
                {
                    let (y1, y1_s) = self.y1_extra_sign();
                    return y1.xor_sign(sign^y1_s)
                }
        
                let (b, b_s) = if self > Self::from(8.1481439053379443450737827536375e90)
                {
                    /* x > 2**302 */
                    /* (x >> n**2)
                     *      Jn(x) = cos(x-(2n+1)*pi/4)*sqrt(2/x*pi)
                     *      Yn(x) = sin(x-(2n+1)*pi/4)*sqrt(2/x*pi)
                     *      Let s=sin(x), c=cos(x),
                     *          xn=x-(2n+1)*pi/4, sqt2 = sqrt(2),then
                     *
                     *             n    sin(xn)*sqt2    cos(xn)*sqt2
                     *          ----------------------------------
                     *             0     s-c             c+s
                     *             1    -s-c            -c+s
                     *             2    -s+c            -c-s
                     *             3     s+c             c-s
                     */
                    let (s, s_s) = self.sin_extra_sign();
                    let (c, c_s) = self.cos_extra_sign();
                    let (n1, n2) = (nm1 & 0b01 == 0, nm1 & 0b10 == 0);
                    let s1 = s_s ^ n2;
                    let s2 = c_s ^ n1 ^ n2;
                    let (sc, b_s) = s.add_with_sign_extra_sign(s1, c, s2);
                    (Self::from(INVSQRTPI)*sc/self.sqrt(), b_s)
                }
                else
                {
                    let h = Self::from(2u8)/self;
                    let (mut a, mut a_s) = self.y0_extra_sign();
                    let (mut b, mut b_s) = self.y1_extra_sign();
                    /* quit if b is -inf */
                    for i in 1..=nm1
                    {
                        if !b.is_finite()
                        {
                            break
                        }
                        (b, b_s) = (b*Self::from(i)*h).add_with_sign_extra_sign(
                            b_s,
                            core::mem::replace(&mut a, b),
                            !core::mem::replace(&mut a_s, b_s)
                        );
                    }
                    (b, b_s)
                };
            
                b.xor_sign(sign^b_s)
            }
        }
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> + ConstZero
{
    /// The additive identity element of `Self`, `0`.
    pub const ZERO: Self = Self::from_bits(U::ZERO);
}

#[cfg(test)]
mod test
{
    #![allow(unused)]

    use std::process::Termination;

    use num::zero;
    use test::Bencher;

    use crate::{ieee754::{FpDouble, FpHalf}, tests::{self, bench_op1, test_op1, test_op1_for}, Fp};

    #[test]
    fn test_ln_gamma()
    {
        test_op1!("ln_gamma", |x| x.ln_gamma().0, None, Some(-4.5..20.0));
    }
    #[bench]
    fn bench_ln_gamma(bencher: &mut Bencher) -> impl Termination
    {
        test_ln_gamma();
        bench_op1!(bencher, Fp::ln_gamma)
    }

    #[test]
    fn test_gamma()
    {
        test_op1!("gamma", |x| x.gamma(), None, Some(-4.5..7.0))
    }
    #[bench]
    fn bench_gamma(bencher: &mut Bencher) -> impl Termination
    {
        test_gamma();
        bench_op1!(bencher, Fp::gamma)
    }
    
    #[test]
    fn test_j0()
    {
        test_op1!("j0", libm::j0f, Fp::j0, Some(0.1), Some(-20.0..20.0))
    }
    #[bench]
    fn bench_j0(bencher: &mut Bencher) -> impl Termination
    {
        test_j0();
        bench_op1!(bencher, Fp::j0)
    }
    
    #[test]
    fn test_y0()
    {
        test_op1!("y0", libm::y0f, Fp::y0, Some(0.1), Some(0.01..20.0))
    }
    #[bench]
    fn bench_y0(bencher: &mut Bencher) -> impl Termination
    {
        test_y0();
        bench_op1!(bencher, Fp::y0)
    }
    
    #[test]
    fn test_j1()
    {
        test_op1!("j1", libm::j1f, Fp::j1, Some(0.1), Some(-20.0..20.0))
    }
    #[bench]
    fn bench_j1(bencher: &mut Bencher) -> impl Termination
    {
        test_j1();
        bench_op1!(bencher, Fp::j1)
    }
    
    #[test]
    fn test_y1()
    {
        test_op1!("y1", libm::y1f, Fp::y1, Some(0.1), Some(0.1..20.0))
    }
    #[bench]
    fn bench_y1(bencher: &mut Bencher) -> impl Termination
    {
        test_y1();
        bench_op1!(bencher, Fp::y1)
    }
    
    #[test]
    fn test_j2()
    {
        test_op1!("j2", |x| libm::jnf(2, x), |x| x.jn(2), Some(0.1), Some(-20.0..20.0))
    }
    #[bench]
    fn bench_j2(bencher: &mut Bencher) -> impl Termination
    {
        test_j2();
        bench_op1!(bencher, |x| x.jn(2))
    }
    
    #[test]
    fn test_y2()
    {
        test_op1!("y2", |x| libm::ynf(2, x), |x| x.yn(2), Some(0.1), Some(1.0..20.0))
    }
    #[bench]
    fn bench_y2(bencher: &mut Bencher) -> impl Termination
    {
        test_y2();
        bench_op1!(bencher, |x| x.yn(2))
    }
    
    #[test]
    fn test_j3()
    {
        test_op1!("j3", |x| libm::jnf(3, x), |x| x.jn(3), Some(0.1), Some(-20.0..20.0))
    }
    #[bench]
    fn bench_j3(bencher: &mut Bencher) -> impl Termination
    {
        test_j3();
        bench_op1!(bencher, |x| x.jn(3))
    }
    
    #[test]
    fn test_y3()
    {
        test_op1!("y3", |x| libm::ynf(3, x), |x| x.yn(3), Some(0.1), Some(1.0..20.0))
    }
    #[bench]
    fn bench_y3(bencher: &mut Bencher) -> impl Termination
    {
        test_y3();
        bench_op1!(bencher, |x| x.yn(3))
    }
    
    #[test]
    fn test_erf()
    {
        test_op1!("erf", libm::erff, Fp::erf, Some(0.1), Some(-5.0..5.0))
    }
    #[bench]
    fn bench_erf(bencher: &mut Bencher) -> impl Termination
    {
        test_erf();
        bench_op1!(bencher, Fp::erf)
    }

    #[test]
    fn test_erfc_once()
    {
        type F = FpDouble;

        let x = F::from(0.9);

        let y = x.erfc();

        println!("{y}")
    }
    
    #[test]
    fn test_erfc()
    {
        test_op1!("erfc", libm::erfcf, Fp::erfc, Some(0.1), Some(-5.0..5.0))
    }
    #[bench]
    fn bench_erfc(bencher: &mut Bencher) -> impl Termination
    {
        test_erfc();
        bench_op1!(bencher, Fp::erfc)
    }

    #[test]
    fn test_next_up_down()
    {
        test_op1!("next_up", |x| x.next_up(), None, Some(-10.0..10.0));
        test_op1!("next_down", |x| x.next_down(), None, Some(-10.0..10.0));

        type F = Fp<u8, true, 3, 1, 3, {usize::MAX}>;

        let mut x = F::neg_infinity();

        loop
        {
            let y = x.next_up();
            if !y.is_zero()
            {
                if x == y
                {
                    let y = x.next_up();
                }
                assert_ne!(x, y);
            }
            assert_eq!(-(-x).next_down(), y);
            if !(x == y.next_down())
            {
                let yy = x.next_up();
                let xx = y.next_down();
                println!("{x} ^ {yy} v {xx}")
            }
            assert_eq!(x, y.next_down());
            x = y;
            if !x.is_finite()
            {
                break
            }
        }

        loop
        {
            let y = x.next_down();
            if !y.is_zero()
            {
                assert_ne!(x, y);
            }
            assert_eq!(-(-x).next_up(), y);
            assert_eq!(x, y.next_up());
            x = y;
            if !x.is_finite()
            {
                break
            }
        }
    }
    #[bench]
    fn bench_next_up_down(bencher: &mut Bencher) -> impl Termination
    {
        test_next_up_down();
        bench_op1!(bencher, Fp::next_up);
        bench_op1!(bencher, Fp::next_down)
    }
}