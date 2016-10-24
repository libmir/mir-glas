module glas.internal.utility;
pragma(LDC_no_moduleinfo);

import std.traits: Unqual;
import std.meta: AliasSeq;
import glas.common;

enum bool isComplex(C)
     = is(Unqual!C == creal)
    || is(Unqual!C == cdouble)
    || is(Unqual!C == cfloat);

enum bool isReal(C)
     = is(Unqual!C == real)
    || is(Unqual!C == double)
    || is(Unqual!C == float);

enum bool isImaginary(C)
     = is(Unqual!C == ireal)
    || is(Unqual!C == idouble)
    || is(Unqual!C == ifloat);


template realType(C)
{
    static if (isComplex!C)
        alias realType = typeof(Unqual!C.init.re);
    else
    static if (isImaginary!C)
        alias realType = typeof(Uqual!C.init.im);
    else
        alias realType = Unqual!C;
}

template imaginaryType(C)
{
    static if (isComplex!C)
        alias imaginaryType = typeof(Unqual!C.init.re * 1fi);
    else
    static if (isImaginary!C)
        alias imaginaryType = Unqual!C;
    else
        alias imaginaryType = typeof(Unqual!C.init * 1fi);
}

alias Iota(size_t j) = Iota!(0, j);

template Iota(size_t i, size_t j)
{
    static assert(i <= j, "Iota: i should be less than or equal to j");
    static if (i == j)
        alias Iota = AliasSeq!();
    else
        alias Iota = AliasSeq!(i, Iota!(i + 1, j));
}


static if (__VERSION__ < 2072)
{
    import std.experimental.ndslice.slice: Slice;
    import ldc.attributes : fastmath;

    @fastmath:
    pragma(inline, true)
    @property T* ptr(T)(Slice!(2, T*) slice)
    {
        return &(slice.front.front());
    }

    static if (__VERSION__ < 2072)
    pragma(inline, true)
    @property T* ptr(T)(Slice!(1, T*) slice)
    {
        return &(slice.front());
    }

    pragma(inline, true)
    auto _toSlice(size_t N, T)(size_t[N] lengths, sizediff_t[N] strides, T ptr)
    {
        static union U
        {
            Slice!(N, T) slice;
            struct
            {
                size_t[N] lengths;
                sizediff_t[N] strides;
                T ptr;
            }
        }
        U ret = void;
        ret.lengths = lengths;
        ret.strides = strides;
        ret.ptr = ptr;
        return ret.slice;
    }
}

template prefix(T)
{
    static if (is(T == float))
        enum prefix = "s";
    else
    static if (is(T == double))
        enum prefix = "d";
    else
    static if (is(T == cfloat))
        enum prefix = "c";
    else
    static if (is(T == cdouble))
        enum prefix = "z";
    else static assert(0);
}
