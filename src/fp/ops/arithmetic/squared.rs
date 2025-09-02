use core::num::FpCategory;

use crate::{util::Saturation, Fp, FpRepr, fp::as_lossless};

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub(super) fn mantissa_squared(mut mantissa: U, exp: &mut U) -> Result<U, Self>
    {
        let mut shifts = FRAC_SIZE;
        shifts -= Self::shr_mantissa_without_loss::<_, usize>(&mut mantissa, Some(shifts/2), 1, None)*2;
        mantissa = Self::integral_mul(mantissa, mantissa, exp)?;
        match Self::mantissa_shr(&mut mantissa, shifts, exp)
        {
            Ok(()) => (),
            Err(sat) => match sat
            {
                Saturation::Overflow => return Err(Self::infinity()),
                Saturation::Underflow => return Err(Self::zero())
            }
        }
        Ok(mantissa)
    }

    pub fn squared(mut self) -> Self
    {
        as_lossless!(
            [self],
            |[x]| [x*x],
            {
                self = self.abs();
                match self.classify()
                {
                    FpCategory::Nan | FpCategory::Infinite | FpCategory::Zero => self,
                    FpCategory::Normal | FpCategory::Subnormal => {
                        if self.is_one()
                        {
                            return self
                        }

                        let mut e: U = self.exp_bits();
                        let mut f: U = self.mantissa_bits();
                
                        e = match Self::exponent_add(e, e, &mut f, None)
                        {
                            Ok(e) => e,
                            Err(done) => return done
                        };
                        f = match Self::mantissa_squared(f, &mut e)
                        {
                            Ok(e) => e,
                            Err(done) => return done
                        };
                
                        Self::normalize_mantissa(&mut e, &mut f, None);
                        Self::from_exp_mantissa(e, f)
                    }
                }
            }
        )
    }
}