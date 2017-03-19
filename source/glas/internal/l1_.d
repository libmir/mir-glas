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

    static if (isComplex!T)
    {
        alias R = realType!T;
        alias I = imaginaryType!T;
    }

    //////////// copy
    void _glas_} ~ prefix!Type ~ q{copy(
        size_t n,
        ptrdiff_t incx,
        const(T)* x,
        ptrdiff_t incy,
        T* y,
        )
    {
        glas.internal.l1.copy(n, incx, x, incy, y);
    }

    void glas_} ~ prefix!Type ~ q{copy(
        Slice!(SliceKind.universal, [1], const(T)*) xsl,
        Slice!(SliceKind.universal, [1], T*) ysl,
        )
    {
        glas.ndslice.copy(xsl.length, xsl._stride, xsl._iterator, ysl._stride, ysl._iterator);
    }

    int } ~ prefix!Type ~ q{copy_(
        ref const FortranInt n,
        const(T)* x,
        ref const FortranInt incx,
        T* y,
        ref const FortranInt incy,
        )
    {
        if(incx < 0)
            x -= (n - 1) * incx;
        if(incy < 0)
            y -= (n - 1) * incy;

        glas.ndslice.copy(n, incx, x, incy, y);
        return 0;
    }

    //////////// swap
    void _glas_} ~ prefix!Type ~ q{swap(
        size_t n,
        ptrdiff_t incx,
        T* x,
        ptrdiff_t incy,
        T* y,
        )
    {
        glas.internal.l1.swap(n, incx, x, incy, y);
    }

    void glas_} ~ prefix!Type ~ q{swap(
        Slice!(SliceKind.universal, [1], T*) xsl,
        Slice!(SliceKind.universal, [1], T*) ysl,
        )
    {
        glas.ndslice.swap(xsl.length, xsl._stride, xsl._iterator, ysl._stride, ysl._iterator);
    }

    int } ~ prefix!Type ~ q{swap_(
        ref const FortranInt n,
        T* x,
        ref const FortranInt incx,
        T* y,
        ref const FortranInt incy,
        )
    {
        if(incx < 0)
            x -= (n - 1) * incx;
        if(incy < 0)
            y -= (n - 1) * incy;

        glas.ndslice.swap(n, incx, x, incy, y);
        return 0;
    }

    //////////// axpy
    void _glas_} ~ prefix!Type ~ q{axpy(
        T a,
        size_t n,
        ptrdiff_t incx,
        const(T)* x,
        ptrdiff_t incy,
        T* y,
        )
    {
        glas.internal.l1.axpy(a, n, incx, x, incy, y);
    }

    void glas_} ~ prefix!Type ~ q{axpy(
        T a,
        Slice!(SliceKind.universal, [1], const(T)*) xsl,
        Slice!(SliceKind.universal, [1], T*) ysl,
        )
    {
        glas.ndslice.axpy(a, xsl.length, xsl._stride, xsl._iterator, ysl._stride, ysl._iterator);
    }

    int } ~ prefix!Type ~ q{axpy_(
        ref const FortranInt n,
        ref const T a,
        const(T)* x,
        ref const FortranInt incx,
        T* y,
        ref const FortranInt incy,
        )
    {
        if(incx < 0)
            x -= (n - 1) * incx;
        if(incy < 0)
            y -= (n - 1) * incy;

        glas.ndslice.axpy(a, n, incx, x, incy, y);
        return 0;
    }

    //////////// scal
    void _glas_} ~ prefix!Type ~ q{scal(
        T a,
        size_t length,
        ptrdiff_t stride,
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
    void _glas_} ~ prefix!Type ~ prefix!(realType!Type) ~ q{scal(
        R a,
        size_t length,
        ptrdiff_t stride,
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
        ptrdiff_t stride,
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