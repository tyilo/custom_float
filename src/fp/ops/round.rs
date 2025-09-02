use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    /// Returns the largest integer less than or equal to `self`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.99);
    /// let g = FpDouble::from(3.0);
    ///
    /// assert_eq!(f.floor(), FpDouble::from(3.0));
    /// assert_eq!(g.floor(), FpDouble::from(3.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn floor(self) -> Self
    {
        //TODO: Don't do this!
        if !self.is_finite()
        {
            return self
        }
        let one = U::one();
        let mut m = self % one;
        if self.is_sign_negative() && !m.is_zero()
        {
            m += one
        }
        self - m
    }

    /// Returns the smallest integer greater than or equal to `self`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.01);
    /// let g = FpDouble::from(4.0);
    ///
    /// assert_eq!(f.ceil(), FpDouble::from(4.0));
    /// assert_eq!(g.ceil(), FpDouble::from(4.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn ceil(mut self) -> Self
    {
        //TODO: Don't do this!
        if !self.is_finite()
        {
            return self
        }
        let one = U::one();
        let mut m = self % one;
        if self.is_sign_positive() && !m.is_zero()
        {
            if !SIGN_BIT
            {
                self += one
            }
            else
            {
                m -= one
            }
        }
        self - m
    }

    /// Returns the nearest integer to `self`. If a value is half-way between two
    /// integers, round away from `0.0`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.3);
    /// let g = FpDouble::from(-3.3);
    ///
    /// assert_eq!(f.round(), FpDouble::from(3.0));
    /// assert_eq!(g.round(), FpDouble::from(-3.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn round(mut self) -> Self
    {
        //TODO: Don't do this!
        if !self.is_finite()
        {
            return self
        }
        let one = U::one();
        let half = Self::from(0.5);
        let mut m = self % one;
        if self.is_sign_positive() && m >= half
        {
            if !SIGN_BIT
            {
                self += one
            }
            else
            {
                m -= one
            }
        }
        if self.is_sign_negative() && m <= -half
        {
            m += one
        }
        self - m
    }
    
    /// Returns the nearest integer to a number. Rounds half-way cases to the number
    /// with an even least significant digit.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpSingle;
    ///
    /// let f = FpSingle::from(3.3);
    /// let g = FpSingle::from(-3.3);
    /// let h = FpSingle::from(3.5);
    /// let i = FpSingle::from(4.5);
    ///
    /// assert_eq!(f.round_ties_even(), FpSingle::from(3.0));
    /// assert_eq!(g.round_ties_even(), FpSingle::from(-3.0));
    /// assert_eq!(h.round_ties_even(), FpSingle::from(4.0));
    /// assert_eq!(i.round_ties_even(), FpSingle::from(4.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn round_ties_even(mut self) -> Self
    {
        //TODO: Don't do this!
        if !self.is_finite()
        {
            return self
        }
        let onef = Self::one();
        let one = U::one();
        let two = one + one;
        let t = self % two;
        let mut m = self % one;
        let half = Self::from(0.5);
        if self.is_sign_positive() && m >= half && (m != half || t > onef)
        {
            if !SIGN_BIT
            {
                self += one
            }
            else
            {
                m -= one
            }
        }
        if self.is_sign_negative() && m <= -half && (m != -half || t < -onef)
        {
            m += one
        }
        self - m
    }

    /// Returns the integer part of `self`.
    /// This means that non-integer numbers are always truncated towards zero.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let f = FpDouble::from(3.3);
    /// let g = FpDouble::from(-3.7);
    ///
    /// assert_eq!(f.trunc(), FpDouble::from(3.0));
    /// assert_eq!(g.trunc(), FpDouble::from(-3.0));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn trunc(self) -> Self
    {
        //TODO: Don't do this!
        if !self.is_finite()
        {
            return self
        }
        let one = U::one();
        let m = self % one;
        self - m
    }

    /// Returns the fractional part of a number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use custom_float::ieee754::FpDouble;
    ///
    /// let x = FpDouble::from(3.5);
    /// let y = FpDouble::from(-3.5);
    /// let abs_difference_x = (x.fract() - FpDouble::from(0.5)).abs();
    /// let abs_difference_y = (y.fract() - FpDouble::from(-0.5)).abs();
    ///
    /// assert!(abs_difference_x < FpDouble::from(1e-10));
    /// assert!(abs_difference_y < FpDouble::from(1e-10));
    /// ```
    #[must_use = "method returns a new number and does not mutate the original value"]
    #[inline]
    pub fn fract(self) -> Self
    {
        //TODO: Don't do this!
        self - self.trunc()
    }
}