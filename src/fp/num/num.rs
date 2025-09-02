use core::str::FromStr;

use num_traits::Num;

use crate::{Fp, FpRepr};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Num for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type FromStrRadixErr = <Self as FromStr>::Err;

    fn from_str_radix(src: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr>
    {
        Self::from_str_radix(src, radix)
    }
}