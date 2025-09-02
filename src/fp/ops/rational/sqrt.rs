use core::num::FpCategory;

use crate::{fp::NEWTON_RT, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the square root of a number.
    ///
    /// Returns NaN if `self` is a negative number other than `-0.0`.
    ///
    /// This implementation is based on the fast sqrt described in: https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Approximations_that_depend_on_the_floating_point_representation
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let positive = FpDouble::from(4.0);
    /// let negative = FpDouble::from(-4.0);
    ///
    /// let abs_difference = (positive.sqrt() - FpDouble::from(2.0)).abs();
    ///
    /// assert!(abs_difference < FpDouble::from(1e-10));
    /// assert!(negative.sqrt().is_nan());
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn sqrt(self) -> Self
    {
        trait SqrtSpec: Sized
        {
            fn _y(self) -> Self; 
        }
        impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SqrtSpec for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
        where
            U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
        {
            default fn _y(self) -> Self
            {
                self.powf_generic::<false>(Self::from(0.5))
            }
        }
        #[allow(clippy::identity_op)]
        impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SqrtSpec for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
        where
            U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> + FpRepr<SIGN_BIT, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>
        {
            fn _y(self) -> Self
            {
                // TODO: Work for bases that are squares of 2
                if EXP_BASE != 2
                {
                    self.powf_generic::<false>(Self::from(0.5))
                }
                else if !Self::IS_INT_IMPLICIT
                {
                    let x = Fp::<U, SIGN_BIT, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::from_fp(self);
                    Self::from_fp(Fp::<U, SIGN_BIT, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::from_bits(
                        ((Fp::<U, SIGN_BIT, EXP_SIZE, 0, FRAC_SIZE, EXP_BASE>::exp_bias() + U::one()) << (FRAC_SIZE - 1))
                        + (x.to_bits() >> 1usize)
                        - (U::one() << (FRAC_SIZE - 1))
                    ))
                }
                else
                {
                    Self::from_bits(
                        ((Self::exp_bias() + U::one()) << (FRAC_SIZE - 1))
                        + (self.to_bits() >> 1usize)
                        - (U::one() << (FRAC_SIZE - 1))
                    )
                }
            }
        }

        match self.classify()
        {
            FpCategory::Zero => self.abs(),
            FpCategory::Nan => self,
            FpCategory::Infinite => {
                if self.is_sign_negative()
                {
                    return Self::snan()
                }
                self
            },
            FpCategory::Normal | FpCategory::Subnormal => {
                if self.is_sign_negative()
                {
                    return Self::snan()
                }
                let y = self._y();
        
                const NEWTON: usize = NEWTON_RT;
                let half = Self::from(0.5);
                let mut y = y;
                for _ in 0..NEWTON
                {
                    let y_ = half*(y + self/y);
                    if !y_.is_finite()
                    {
                        break
                    }
                    y = y_;
                }
                y
            }
        }
    }
}