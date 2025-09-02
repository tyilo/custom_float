use core::{cmp::Ordering, num::FpCategory};

use num_traits::Pow;

use crate::{util, Fp, FpRepr, Int, UInt};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Raises a number to an integer power.
    ///
    /// Using this function is generally faster than using `powf`.
    /// It might have a different sequence of rounding operations than `powf`,
    /// so the results are not guaranteed to agree.
    /// 
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.23);
    /// 
    /// let abs_difference = (x.powi(2) - x*x).abs();
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// let abs_difference = (x.powi(-2) - (x*x).recip()).abs();
    /// assert!(abs_difference < FpDouble::from(1e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn powi<I: Int>(self, n: I) -> Self
    {
        util::powi(self, n)
    }
    
    /// Raise a number to an unsigned integer power.
    ///
    /// Using this function is generally faster than using `powf`.
    /// It might have a different sequence of rounding operations than `powf`,
    /// so the results are not guaranteed to agree.
    /// 
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.23);
    /// let n = 2u32;
    /// let abs_difference = (x.powu(n) - x*x).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn powu<I: UInt>(self, n: I) -> Self
    {
        util::powu(self, n)
    }

    pub(crate) fn powf_generic<const SPECIALIZE: bool>(self, n: Self) -> Self
    {
        let explode = |s| match s
        {
            true => Self::zero(),
            false => Self::infinity()
        };

        let pow_sign = |x: Self, s| match s
        {
            true => x.recip(),
            false => x
        };

        let n_s = n.is_sign_negative();
        let nabs = n.abs();
        let one = Self::one();
        let uone = U::one();
        let utwo = uone + uone;

        match (self.classify(), n.classify())
        {
            (_, FpCategory::Zero) => one,
            (FpCategory::Normal, _) if self == one => self,
            (FpCategory::Nan, _) | (_, FpCategory::Nan) => self.add_nan(n),
            (_, FpCategory::Normal) if nabs == one => pow_sign(self, n_s),
            (FpCategory::Zero, FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal) => explode(!n_s),
            /*(FpCategory::Zero, FpCategory::Normal | FpCategory::Subnormal) => {
                let x_s = n_s;

                // if n is an odd integer
                let n_d = (nabs.trunc() % utwo + uone)/utwo;
                let noi = n_d.trunc() == n_d;
                let sign = noi && x_s;

                explode(!n_s).with_sign(sign)
            },*/
            (FpCategory::Infinite, FpCategory::Infinite | FpCategory::Normal | FpCategory::Subnormal) => explode(n_s),
            (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Infinite) => {
                let xabs = self.abs();

                match xabs.total_cmp(one)
                {
                    Ordering::Greater => explode(n_s),
                    Ordering::Equal => one,
                    Ordering::Less => explode(!n_s)
                }
            }
            /*(FpCategory::Infinite, FpCategory::Normal | FpCategory::Subnormal) => {
                let x_s = self.is_sign_negative();
                
                // if n is an odd integer
                let n_d = (nabs.trunc() % utwo + uone)/utwo;
                let noi = n_d.trunc() == n_d;
                let sign = noi && x_s;

                explode(n_s).with_sign(sign)
            },*/
            (FpCategory::Normal | FpCategory::Subnormal, FpCategory::Normal | FpCategory::Subnormal) => {
                let x_s = self.is_sign_negative();
                let xabs = self.abs();
                
                // if n is an odd integer
                let n_d = (nabs.trunc() % utwo + uone)/utwo;
                let noi = n_d.trunc() == n_d;

                let edge_x = {
                    Self::from_bits(Self::shift_exp(Self::max_exponent_bits() - U::one()) | Self::shift_int(U::one()))
                };
                let edge_n = U::from(util::exp2_ilog(FRAC_SIZE + INT_SIZE, EXP_BASE))
                    .map(|exp_frac| {
                        Self::from_bits(Self::shift_exp(Self::exp_bias() + exp_frac) | Self::shift_int(U::one()))
                    });
        
                if (xabs >= edge_x) || edge_n.is_some_and(|edge_n| nabs >= edge_n)
                {
                    if noi && let Some(i) = {
                        let n_clamp = Self::from_int(U::max_value());
                        nabs.minimum(n_clamp).to_int::<u128>()
                    }
                    {
                        pow_sign(self, n_s).powu(i)
                    }
                    else
                    {
                        pow_sign(self, n_s && !x_s).sqrt()
                    }
                }
                else if x_s && n.trunc() != n
                {
                    Self::snan()
                }
                else if SPECIALIZE && (nabs/utwo).is_one()
                {
                    pow_sign(self, n_s && !x_s).sqrt()
                }
                else if SPECIALIZE && (nabs/(utwo + uone)).is_one()
                {
                    pow_sign(self, n_s && !x_s).cbrt()
                }
                else if !SIGN_BIT && xabs < one
                {
                    let xabs_log = xabs.recip().logb();
        
                    let n_xabs_log = xabs_log*n;
                    
                    n_xabs_log.expb().recip()
                }
                else
                {
                    let xabs_log = xabs.logb();
        
                    let n_xabs_log = xabs_log*n;
                    
                    n_xabs_log.expb()
                }
            }
        }
    }

    /// Raises a number to a floating point power.
    ///
    /// This implementation is based on the [Apple Libm-315 implementation of powf](https://opensource.apple.com/source/Libm/Libm-315/Source/ARM/powf.c)
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.23);
    /// 
    /// let abs_difference = (x.powf(FpDouble::from(2.0/3.0)) - (x*x).cbrt()).abs();
    /// assert!(abs_difference < FpDouble::from(1e-4));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn powf(self, n: Self) -> Self
    {
        self.powf_generic::<true>(n)
    }

    /// Returns `self` to the power `rhs`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(1.23);
    /// 
    /// let abs_difference = (x.pow(FpDouble::from(2.0/3.0)) - (x*x).cbrt()).abs();
    /// assert!(abs_difference < FpDouble::from(1e-4));
    /// 
    /// let abs_difference = (x.pow(2u32) - x*x).abs();
    /// assert!(abs_difference < FpDouble::from(1e-5));
    /// 
    /// let abs_difference = (x.pow(-2i32) - (x*x).recip()).abs();
    /// assert!(abs_difference < FpDouble::from(1e-3));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    pub fn pow<Rhs>(self, n: Rhs) -> <Self as Pow<Rhs>>::Output
    where
        Self: Pow<Rhs>
    {
        Pow::pow(self, n)
    }
}

macro_rules! impl_powu {
    ($($i:ty),*) => {
        $(
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Pow<$i> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                type Output = Self;
            
                #[inline]
                fn pow(self, rhs: $i) -> Self::Output
                {
                    self.powu(rhs)
                }
            }
        )*
    };
}

impl_powu!(u8, u16, u32, usize, u64, u128/*, U256*/);

macro_rules! impl_powi {
    ($($i:ty),*) => {
        $(
            impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Pow<$i> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            where
                U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
            {
                type Output = Self;
            
                #[inline]
                fn pow(self, rhs: $i) -> Self::Output
                {
                    self.powi(rhs)
                }
            }
        )*
    };
}

impl_powi!(i8, i16, i32, isize, i64, i128/*, I256*/);

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Pow<Self> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn pow(self, rhs: Self) -> Self::Output
    {
        self.powf(rhs)
    }
}