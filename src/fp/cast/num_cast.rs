use num_traits::{NumCast, ToPrimitive};

use crate::{fp::Fps, Fp, FpRepr, UInt, AnyInt};

macro_rules! spec_override {
    (
        $(<{$($g:tt)+}>)? $trait:ident($self:ident) for $for:ty
        $(where {$($where:tt)*})?
        {
            fn $cast:ident
            $cast_body:block

            fn $cast_extra_sign:ident
            $cast_extra_sign_body:block
        }
    ) => {
        impl<
            U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize,
            $($($g)*)?
        > $trait<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> for $for
        where
            U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>,
            Self: ToPrimitive,
            $($($where)*)?
        {
            fn $cast($self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
            $cast_body
            
            fn $cast_extra_sign($self) -> Option<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
            $cast_extra_sign_body
        }
    };
    (
        $(<{$($g:tt)+}>)? $trait:ident($self:ident)
        $(where {$($where:tt)*})?
        {
            fn $cast:ident
            $cast_body:block

            fn $cast_extra_sign:ident
            $cast_extra_sign_body:block
        }
    ) => {
        spec_override!(
            <{$($($g)*)?}> $trait($self) for T
            $(where {$($where)*})?
            {
                fn $cast
                $cast_body

                fn $cast_extra_sign
                $cast_extra_sign_body
            }
        );
    };
}

macro_rules! specs {
    (
        $vis:vis $trait:ident($self:ident)
        {
            fn $cast:ident
            $cast_body:block

            fn $cast_extra_sign:ident
            $cast_extra_sign_body:block
        }
    ) => {
        $vis trait $trait<
            U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize
        >: ToPrimitive
        where
            U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
        {
            fn $cast($self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>;
            fn $cast_extra_sign($self) -> Option<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>;
        }
        impl<
            T: ToPrimitive, U: UInt, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize
        > $trait<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE> for T
        where
            U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
        {
            default fn $cast($self) -> Option<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
            $cast_body
            
            default fn $cast_extra_sign($self) -> Option<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
            $cast_extra_sign_body
        }
    };
    (
        $vis:vis $trait:ident($self:ident)
        {
            fn $cast:ident
            $cast_body:block

            fn $cast_extra_sign:ident
            $cast_extra_sign_body:block
        }
        $(<{$($o_g:tt)+}>)? $(for $o_for:ty)?
        $(where {$($o_where:tt)+})?
        {
            $o_cast_body:block

            $o_cast_extra_sign_body:block
        }
    ) => {
        specs!(
            $vis $trait($self)
            {
                fn $cast
                $cast_body
    
                fn $cast_extra_sign
                $cast_extra_sign_body
            }
        );
        spec_override!(
            $(<{$($o_g)+}>)? $trait($self) $(for $o_for)?
            $(where {$($o_where)+})?
            {
                fn $cast
                $o_cast_body
    
                fn $cast_extra_sign
                $o_cast_extra_sign_body
            }
        );
    };
    (
        $vis:vis $trait:ident($self:ident)
        {
            fn $cast:ident
            $cast_body:block

            fn $cast_extra_sign:ident
            $cast_extra_sign_body:block
        }
        $(<{$($n_g:tt)+}>)? $(for $n_for:ty)?
        $(where {$($n_where:tt)+})?
        {
            $n_cast_body:block

            $n_cast_extra_sign_body:block
        }
        $($(<{$($o_g:tt)+}>)? $(for $o_for:ty)?
        $(where {$($o_where:tt)+})?
        {
            $o_cast_body:block

            $o_cast_extra_sign_body:block
        })+
    ) => {
        mod private
        {
            use super::*;
            specs!(
                pub(super) $trait($self)
                {
                    fn $cast
                    $cast_body
        
                    fn $cast_extra_sign
                    $cast_extra_sign_body
                }
                $($(<{$($o_g)+}>)? $(for $o_for)?
                $(where {$($o_where)+})?
                {
                    $o_cast_body
        
                    $o_cast_extra_sign_body
                })+
            );
        }
        specs!(
            $vis $trait($self)
            {
                fn $cast
                {
                    <T as private::$trait::<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>::$cast($self)
                }
    
                fn $cast_extra_sign
                {
                    <T as private::$trait::<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>::$cast_extra_sign($self)
                }
            }
            $(<{$($n_g)+}>)? $(for $n_for)?
            $(where {$($n_where)+})?
            {
                $n_cast_body
    
                $n_cast_extra_sign_body
            }
        );
    };
}

mod into_primitive
{
    mod private
    {
        use core::convert::Infallible;

        use num_traits::AsPrimitive;

        pub(super) trait IntoPrimitiveSpec<To>
        {
            type IsImpl;
    
            fn into_primitive(self) -> To;
        }
        impl<From, To> IntoPrimitiveSpec<To> for From
        {
            default type IsImpl = Infallible;

            default fn into_primitive(self) -> To
            {
                unreachable!()
            }
        }
        impl<From, To> IntoPrimitiveSpec<To> for From
        where
            From: AsPrimitive<To>,
            To: Copy + 'static
        {
            type IsImpl = ();

            fn into_primitive(self) -> To
            {
                self.as_()
            }
        }
    }

    pub(super) trait IntoPrimitiveSpec<To>
    {
        type IsImpl;

        fn into_primitive(self) -> To;
    }
    impl<From, To> IntoPrimitiveSpec<To> for From
    {
        default type IsImpl = <From as private::IntoPrimitiveSpec<To>>::IsImpl;

        default fn into_primitive(self) -> To
        {
            <From as private::IntoPrimitiveSpec<To>>::into_primitive(self)
        }
    }
    impl<From, To> IntoPrimitiveSpec<To> for From
    where
        From: Into<To>
    {
        type IsImpl = ();

        fn into_primitive(self) -> To
        {
            self.into()
        }
    }
}

trait IntoPrimitive<To>
{
    fn into_primitive(self) -> To;
}

impl<From, To> IntoPrimitive<To> for From
where
    From: into_primitive::IntoPrimitiveSpec<To, IsImpl = ()>
{
    fn into_primitive(self) -> To
    {
        <From as into_primitive::IntoPrimitiveSpec<To>>::into_primitive(self)
    }
}

specs!(
    pub NumCastSpec(self)
    {
        fn cast
        {
            self.to_f64().map(|f| f.into())
        }
        fn cast_extra_sign
        {
            self.to_f64().map(|f| f.into())
        }
    }
    <{T}> where {
        T: IntoPrimitive<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>> + IntoPrimitive<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    }
    {
        {
            Some(self.into_primitive())
        }
        {
            Some(self.into_primitive())
        }
    }
    <{T}> where {
        T: IntoPrimitive<Fps<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    }
    {
        {
            Some(self.into_primitive().collect())
        }
        {
            Some(self.into_primitive())
        }
    }
    <{T}> where {
        T: IntoPrimitive<Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>>
    }
    {
        {
            Some(self.into_primitive())
        }
        {
            Some(self.into_primitive().extra_sign(false))
        }
    }
    <{T}> where {
        T: AnyInt
    }
    {
        {
            Some(Fp::from_int(self))
        }
        {
            Some(Fps::from_int(self))
        }
    }
    <{T, const S: bool, const E: usize, const I: usize, const F: usize, const B: usize}> for Fp<T, S, E, I, F, B>
    where {
        T: FpRepr<S, E, I, F, B>
    }
    {
        {
            Some(Fp::from_fp(self))
        }
        {
            Some(Fps::from_fp(self))
        }
    }
    <{T}> where {
        T: IntoPrimitive<f64>
    }
    {
        {
            Some(Fp::from(self.into_primitive()))
        }
        {
            Some(Fps::from(self.into_primitive()))
        }
    }
    <{T}> where {
        T: IntoPrimitive<f128>
    }
    {
        {
            Some(Fp::from(self.into_primitive()))
        }
        {
            Some(Fps::from(self.into_primitive()))
        }
    }
);

impl<U, const SIGN_BIT: bool, const EXP_SIZE: usize, const INT_SIZE: usize, const FRAC_SIZE: usize, const EXP_BASE: usize> NumCast for Fp<U, SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
where
    U: FpRepr<SIGN_BIT, EXP_SIZE, INT_SIZE, FRAC_SIZE, EXP_BASE>
{
    #[inline]
    fn from<T: ToPrimitive>(n: T) -> Option<Self>
    {
        n.cast()
    }
}