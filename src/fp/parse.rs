use core::str::FromStr;

use num_traits::ParseFloatError;

use crate::{Fp, FpRepr};

fn str_to_ascii_lower_eq_str(a: &str, b: &str) -> bool
{
    a.len() == b.len()
        && a.bytes().zip(b.bytes()).all(|(a, b)| {
            let a_to_ascii_lower = a | ((a.is_ascii_uppercase() as u8) << 5);
            a_to_ascii_lower == b
        })
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> FromStr for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    type Err = ParseFloatError;

    fn from_str(src: &str) -> Result<Self, ParseFloatError>
    {
        Self::from_str_radix(src, 10)
    }
}

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    pub fn from_str_radix(src: &str, radix: u32) -> Result<Self, <Self as FromStr>::Err>
    {
        use num_traits::FloatErrorKind::*;
        use num_traits::ParseFloatError as PFE;

        #[cfg(any(
            feature = "use_std_float",
            all(debug_assertions, test)
        ))]
        macro_rules! from_str_radix_as_float {
            ($($t:ty),*) => {
                $(
                    if crate::util::is_float_conversion_lossless::<$t, Self>()
                    {
                        Some(<$t as num_traits::Num>::from_str_radix(src, radix).map(Self::from))
                    }
                    else
                )*
                {
                    None
                }
            };
        }

        #[cfg(any(
            feature = "use_std_float",
            all(debug_assertions, test)
        ))]
        macro_rules! parse_as_float {
            ($($t:ty),*) => {
                $(
                    if crate::util::is_float_conversion_lossless::<$t, Self>()
                    {
                        Some(
                            src.parse::<$t>()
                                .map(Self::from)
                                .map_err(|_| PFE {
                                    kind: if src.is_empty() { Empty } else { Invalid },
                                })
                        )
                    }
                    else
                )*
                {
                    None
                }
            };
        }

        #[cfg(any(
            feature = "use_std_float",
            all(debug_assertions, test)
        ))]
        let _as_lossless = if radix == 10
        {
            parse_as_float!(
                f16,
                f32,
                f64
            )
        }
        else
        {
            from_str_radix_as_float!(
                f32,
                f64
            )
        };
        #[cfg(all(not(test), feature = "use_std_float"))]
        if let Some(result) = _as_lossless
        {
            return result
        }
        let y = (|| {
            // Special values
            if str_to_ascii_lower_eq_str(src, "inf")
                || str_to_ascii_lower_eq_str(src, "infinity")
            {
                return Ok(Self::infinity());
            }
            else if str_to_ascii_lower_eq_str(src, "-inf")
                || str_to_ascii_lower_eq_str(src, "-infinity")
            {
                return Ok(Self::neg_infinity());
            }
            else if str_to_ascii_lower_eq_str(src, "nan")
            {
                return Ok(Self::nan());
            }
            else if str_to_ascii_lower_eq_str(src, "-nan")
            {
                return Ok(-Self::nan());
            }

            fn slice_shift_char(src: &str) -> Option<(char, &str)>
            {
                let mut chars = src.chars();
                Some((chars.next()?, chars.as_str()))
            }

            let (sign, src) =  match slice_shift_char(src)
            {
                None             => return Err(PFE { kind: Empty }),
                Some(('-', ""))  => return Err(PFE { kind: Empty }),
                Some(('-', src)) => (true, src),
                Some((_, _))     => (false,  src),
            };

            // The significand to accumulate
            let mut sig = Self::zero();
            // Necessary to detect overflow
            let mut prev_sig = sig;
            let mut cs = src.chars().enumerate();
            // Exponent prefix and exponent index offset
            let mut exp_info = None::<(char, usize)>;

            // Parse the integer part of the significand
            for (i, c) in cs.by_ref()
            {
                match c.to_digit(radix)
                {
                    Some(digit) => {
                        // shift significand one digit left
                        // add/subtract current digit depending on sign
                        sig = sig.mul_int(radix).add_int(digit);

                        // Detect overflow by comparing to last value, except
                        // if we've not seen any non-zero digits.
                        if !prev_sig.is_zero()
                        {
                            if sig <= prev_sig
                            {
                                return Ok(Self::infinity().with_sign(sign))
                            }

                            // Detect overflow by reversing the shift-and-add process
                            if prev_sig != sig.sub_int(digit).div_int(radix)
                            {
                                return Ok(Self::infinity().with_sign(sign))
                            }
                        }
                        prev_sig = sig;
                    },
                    None => match c {
                        'e' | 'E' | 'p' | 'P' => {
                            exp_info = Some((c, i + 1));
                            break;  // start of exponent
                        },
                        '.' => {
                            break;  // start of fractional part
                        },
                        _ => {
                            return Err(PFE { kind: Invalid });
                        },
                    },
                }
            }

            // If we are not yet at the exponent parse the fractional
            // part of the significand
            if exp_info.is_none()
            {
                let mut power = Self::one();
                for (i, c) in cs.by_ref()
                {
                    match c.to_digit(radix)
                    {
                        Some(digit) => {
                            // Decrease power one order of magnitude
                            power /= power.div_int(radix);
                            // add/subtract current digit depending on sign
                            sig += power.mul_int(digit);

                            // Detect overflow by comparing to last value
                            if sig < prev_sig
                            {
                                return Ok(Self::infinity().with_sign(sign))
                            }
                            prev_sig = sig;
                        },
                        None => match c {
                            'e' | 'E' | 'p' | 'P' => {
                                exp_info = Some((c, i + 1));
                                break; // start of exponent
                            },
                            _ => {
                                return Err(PFE { kind: Invalid });
                            },
                        },
                    }
                }
            }

            sig = sig.with_sign(sign);

            // Parse and calculate the exponent
            if let Some((c, offset)) = exp_info
            {
                let base = match c
                {
                    'E' | 'e' if radix == 10 => Self::from_int(10u8),
                    'P' | 'p' if radix == 16 => Self::from_int(2u8),
                    _ => return Err(PFE { kind: Invalid }),
                };

                // Parse the exponent as decimal integer
                let src = &src[offset..];
                let (sign_exp, exp) = match slice_shift_char(src)
                {
                    Some(('-', src)) => (true, src.parse::<usize>()),
                    Some(('+', src)) => (false,  src.parse::<usize>()),
                    Some((_, _))     => (false,  src.parse::<usize>()),
                    None             => return Err(PFE { kind: Invalid }),
                };

                sig *= match (sign_exp, exp) {
                    (false, Ok(exp)) => base.powi(exp as i32),
                    (true,  Ok(exp)) => base.powi(exp as i32).recip(),
                    (_, Err(_))      => return Err(PFE { kind: Invalid }),
                }
            };

            Ok(sig)
        })();
        #[cfg(any(
            feature = "use_std_float",
            all(debug_assertions, test)
        ))]
        if let Some(as_lossless) = _as_lossless
        {
            match (&y, as_lossless)
            {
                (&Ok(_yy), Err(_)) | (&Err(_), Ok(_yy)) => {
                    #[cfg(all(debug_assertions, test))]
                    panic!("Parse error mismatch");
                    #[cfg(all(not(test), feature = "use_std_float"))]
                    return Ok(_yy)
                },
                (&Ok(yy), Ok(yl)) => {
                    #[cfg(all(debug_assertions, test))]
                    if (!yy.is_finite() || yl.is_finite()) && !yy.approx_eq(yl.midpoint(yy)) // EXTRA FUZZY (i don't care that much)
                    {
                        debug_assert_eq!(yy, yl, "Error is too big!")
                    }
                    #[cfg(all(not(test), feature = "use_std_float"))]
                    return Ok(yl)
                },
                _ => ()
            }
        }
        y
    }
}

#[cfg(test)]
mod test
{
    use num_traits::Num;
    use test::Bencher;

    use crate::{ieee754::FpSingle, tests::{bench_op1_premap_for, for_all_floats, test_op1_premap}};

    #[test]
    fn test_from_str_radix_once()
    {
        type F = FpSingle;

        let x = F::from(1);
        let s = x.to_string();
        print!("{s} @ {x}");
        let y = F::from_str_radix(&s, 10).unwrap();
        println!(" = {y}")
    }

    #[test]
    fn test_from_str_radix()
    {
        test_op1_premap!(
            "from_str_radix",
            |x| x.to_string(),
            |s: String| Num::from_str_radix(&s, 10).unwrap(),
            Some(0.01),
            Some(-10.0..10.0)
        );
    }

    #[bench]
    fn bench_from_str_radix(bencher: &mut Bencher)
    {
        test_from_str_radix();
        bench_op1_premap_for!(
            F,
            bencher,
            |x| x.to_string(),
            |s| F::from_str_radix(s, 10).unwrap()
        );
    }

    #[test]
    fn from_str_radix_multi_byte_fail()
    {
        for_all_floats!(
            F, _DIR, {
                // Ensure parsing doesn't panic, even on invalid sign characters
                assert!(F::from_str_radix("™0.2", 10).is_err());
        
                // Even when parsing the exponent sign
                assert!(F::from_str_radix("0.2E™1", 10).is_err());
            }
        )
    }

    #[test]
    fn from_str_radix_ignore_case()
    {
        for_all_floats!(
            F, _DIR, {
                assert_eq!(
                    F::from_str_radix("InF", 16).unwrap(),
                    F::infinity()
                );
                assert_eq!(
                    F::from_str_radix("InfinitY", 16).unwrap(),
                    F::infinity()
                );
                assert_eq!(
                    F::from_str_radix("-InF", 8).unwrap(),
                    F::neg_infinity()
                );
                assert_eq!(
                    F::from_str_radix("-InfinitY", 8).unwrap(),
                    F::neg_infinity()
                );
                assert!(F::from_str_radix("nAn", 4).unwrap().is_nan());
                assert!(F::from_str_radix("-nAn", 4).unwrap().is_nan());
            }
        )
    }
}