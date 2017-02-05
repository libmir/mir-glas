/++
Copyright: Copyright © 2016-, Ilya Yaroshenko.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas.internal.l1_;

import mir.ndslice.slice: Slice, SliceKind;
import ldc.intrinsics: llvm_expect;
import glas.ndslice;
import glas.fortran;
import glas.internal.utility;

package(glas) enum L1(Type) =
q{
    pragma(LDC_no_moduleinfo);

    import glas.ndslice;
    import glas.fortran;
    import mir.ndslice.slice: Slice, SliceKind;
    import ldc.intrinsics: llvm_expect;
    import ldc.attributes: fastmath;
    import glas.internal.utility;
    import glas.precompiled.utility;
    import glas.internal.l1;

    private alias T = } ~ Type.stringof ~ q{;

    export extern(C) @system nothrow @nogc @fastmath pragma(inline, false):

    void _glas_} ~ prefix!Type ~ q{scal(
        T a,
        size_t length,
        size_t stride,
        T* ptr,
        )
    {
        glas.internal.l1.scal(a, length, stride, ptr);
    }

    void glas_} ~ prefix!Type ~ q{scal(
        T a,
        Slice!(SliceKind.universal, [1], T*) xsl,
        )
    {
        glas.ndslice.scal(a, xsl.length, xsl._stride, xsl._iterator);
    }

    int } ~ prefix!Type ~ q{scal_(
        ref const FortranInt n,
        ref const T a,
            T* x,
        ref const FortranInt incx,
        )
    {
        glas.ndslice.scal(a, n, incx, x);
        return 0;
    }

    static if (isComplex!T)
    {
        alias R = realType!T;
        alias I = imaginaryType!T;
    }

    static if (isComplex!T)
    void _glas_} ~ prefix!Type ~ prefix!(realType!Type) ~ q{scal(
        R a,
        size_t length,
        size_t stride,
        T* ptr,
        )
    {
        glas.internal.l1.scal(a, length, stride, ptr);
    }

    static if (isComplex!T)
    void glas_} ~ prefix!Type ~ prefix!(realType!Type) ~ q{scal(
        R a,
        Slice!(SliceKind.universal, [1], T*) xsl,
        )
    {
        glas.ndslice.scal(a, xsl.length, xsl._stride, xsl._iterator);
    }

    static if (isComplex!T)
    int } ~ prefix!Type ~ prefix!(realType!Type) ~ q{scal_(
        ref const FortranInt n,
        ref const R a,
            T* x,
        ref const FortranInt incx,
        )
    {
        glas.ndslice.scal(a, n, incx, x);
        return 0;
    }

    static if (isComplex!T)
    void _glas_} ~ prefix!Type ~ prefix!(realType!Type) ~ "I" ~ q{scal(
        I a,
        size_t length,
        size_t stride,
        T* ptr,
        )
    {
        glas.internal.l1.scal(a, length, stride, ptr);
    }

    static if (isComplex!T)
    void glas_} ~ prefix!Type ~ prefix!(realType!Type) ~ q{scal(
        I a,
        Slice!(SliceKind.universal, [1], T*) xsl,
        )
    {
        glas.ndslice.scal(a, xsl.length, xsl._stride, xsl._iterator);
    }

    static if (isComplex!T)
    int } ~ prefix!Type ~ prefix!(realType!Type) ~ "I" ~ q{scal_(
        ref const FortranInt n,
        ref const I a,
            T* x,
        ref const FortranInt incx,
        )
    {
        glas.ndslice.scal(a, n, incx, x);
        return 0;
    }
};