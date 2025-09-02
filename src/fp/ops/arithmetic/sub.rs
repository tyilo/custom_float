use core::ops::{Sub, SubAssign};

use crate::{AnyInt, Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub(crate) fn sub_extra_sign(self, rhs: Self) -> (Self, bool)
    {
        self.add_with_sign_extra_sign(false, rhs, true)
    }

    pub(crate) fn sub_int_extra_sign<I: AnyInt>(self, rhs: I) -> (Self, bool)
    {
        self.add_int_with_sign_extra_sign(false, rhs, true)
    }

    pub fn sub_int<I: AnyInt>(self, rhs: I) -> Self
    {
        self.add_int_with_sign(false, rhs, true)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Sub for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output
    {
        self.add_with_sign(false, rhs, true)
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Sub<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: U) -> Self::Output
    {
        self.sub_int(rhs)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SubAssign for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self)
    {
        *self = *self - rhs
    }
}
impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> SubAssign<U> for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn sub_assign(&mut self, rhs: U)
    {
        *self = *self - rhs
    }
}

#[cfg(test)]
mod test
{
    use core::ops::Sub;
    use test::Bencher;

    use crate::tests::{bench_op2, test_op2};

    #[test]
    fn test_sub()
    {
        test_op2!("sub", Sub::sub, None)
    }
    #[bench]
    fn bench_sub(bencher: &mut Bencher)
    {
        test_sub();
        bench_op2!(bencher, Sub::sub)
    }
}