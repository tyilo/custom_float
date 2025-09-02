use num_traits::Float;

use crate::{fp::Fps, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    pub(crate) fn polynomial<T>(mut self, p: &[T]) -> Self
    where
        T: Into<Self> + Float
    {
        let mut y = Self::zero();
        let mut zn = Self::one();
        for p in p.iter()
        {
            y += p.into()*zn;
            zn *= self;
        }
        y
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    pub(crate) fn polynomial<T>(self, p: &[T]) -> Self
    where
        T: Into<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> + Float
    {
        Fps::from(self).polynomial(p).into()
    }
}