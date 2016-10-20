module glas.precompiled.l3_;

version(LDC)
{
    version(unittest) {} else
    {
        pragma(LDC_no_moduleinfo);
    }
}

import glas.precompiled.utility;
import mir.ndslice.slice: Slice;
import glas.common;
import ldc.attributes: fastmath;
import ldc.intrinsics: llvm_expect;
import glas.precompiled.utility : integer;

package enum L3(string p, string T) =
q{
    version(LDC)
    {
        version(unittest) {} else
        {
            pragma(LDC_no_moduleinfo);
        }
    }

    import glas.common;
    import mir.internal.utility: isComplex;
    import mir.ndslice.slice: Slice;
    import mir.ndslice.iteration: transposed;
    import ldc.attributes: fastmath;
    import ldc.intrinsics: llvm_expect;
    import glas.precompiled.utility;

    private alias T = } ~ T ~ q{;

    export extern(C) @system nothrow @nogc @fastmath pragma(inline, false):

    void glas_} ~ p ~ q{gemm
        (
            T alpha,
                Slice!(2, const(T)*) asl,
                Slice!(2, const(T)*) bsl,
            T beta,
                Slice!(2,        T*) csl,
            Conjugated conja = Conjugated.no,
            Conjugated conjb = Conjugated.no,
        )
    {
        import glas.internal.gemm: gemm_impl, SL3;
        SL3!(T, T, T) arg = void;
        arg.asl = asl;
        arg.bsl = bsl;
        arg.csl = csl;
        arg.alpha_beta[0] = alpha;
        arg.alpha_beta[1] = beta;
        gemm_impl(arg, conja, conjb);
    }

    void glas_} ~ p ~ q{symm
        (
            Side side,
            Uplo uplo,
            T alpha,
                Slice!(2, const(T)*) asl,
                Slice!(2, const(T)*) bsl,
            T beta,
                Slice!(2,        T*) csl,
            Conjugated conja = Conjugated.no,
            Conjugated conjb = Conjugated.no,
        )
    {
        import glas.internal.gemm: SL3;
        import glas.internal.symm: symm_impl;
        SL3!(T, T, T) arg = void;
        arg.asl = asl;
        arg.bsl = bsl;
        arg.csl = csl;
        arg.alpha_beta[0] = alpha;
        arg.alpha_beta[1] = beta;
        symm_impl(arg, side, uplo, conja, conjb);
    }

    int } ~ p ~ q{gemm_(
        ref const char transa,
        ref const char transb,
        ref const integer m,
        ref const integer  n,
        ref const integer k,
        ref const T alpha,
            const(T)* a,
            ref const integer lda,
            const(T)* b,
            ref const integer ldb,
        ref const T beta,
            T* c,
            ref const integer ldc,
        )
    {
        auto  tra = toUpper(transa);
        auto  trb = toUpper(transb);

        auto nota = tra == 'N';
        auto notb = trb == 'N';
        auto conja = cast(Conjugated) (tra == 'C');
        auto conjb = cast(Conjugated) (trb == 'C');

        integer info = void;
        
        if (llvm_expect(!nota && !conja && tra != 'T', false))
            info = 1;
        else
        if (llvm_expect(!notb && !conjb && trb != 'T', false))
            info = 2;
        else
        if (llvm_expect(m < 0, false))
            info = 3;
        else
        if (llvm_expect(n < 0, false))
            info = 4;
        else
        if (llvm_expect(k < 0, false))
            info = 5;
        else
        if (llvm_expect(lda < max(1, nota ? m : k), false))
            info = 8;
        else
        if (llvm_expect(ldb < max(1, notb ? k : n), false))
            info = 10;
        else
        if (llvm_expect(ldc < max(1, m), false))
            info = 13;
        else
        {
            static if (!isComplex!T)
                conja = conjb = Conjugated.no;
            import glas.internal.gemm: gemm_impl, SL3;
            SL3!(T, T, T) arg = void;
            arg.asl = Slice!(2, const(T)*)([m, k], nota ? [1, lda] : [lda, 1], a);
            arg.bsl = Slice!(2, const(T)*)([k, n], notb ? [1, ldb] : [ldb, 1], b);
            arg.csl = Slice!(2,       T *)([m, n],        [1, ldc],            c);
            arg.alpha_beta[0] = alpha;
            arg.alpha_beta[1] = beta;
            gemm_impl(arg, conja, conjb);
            return 0;
        }
        enum name = prefix!T ~ "GEMM";
        xerbla_(name.ptr, &info);
        return 0;
    }

    int } ~ p ~ q{symm_(
        ref const char side,
        ref const char uplo,
        ref const integer m,
        ref const integer  n,
        ref const T alpha,
            const(T)* a,
            ref const integer lda,
            const(T)* b,
            ref const integer ldb,
        ref const T beta,
            T* c,
            ref const integer ldc,
        )
    {
        return symm_impl_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, Conjugated.no);
    }

    static if (isComplex!T)
    int } ~ p ~ q{hemm_(
        ref const char side,
        ref const char uplo,
        ref const integer m,
        ref const integer  n,
        ref const T alpha,
            const(T)* a,
            ref const integer lda,
            const(T)* b,
            ref const integer ldb,
        ref const T beta,
            T* c,
            ref const integer ldc,
        )
    {
        return symm_impl_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, Conjugated.yes);
    }
};


pragma(inline, true)
package auto toUpper()(dchar c)
{
    return dchar(c & 0b1101_1111);
}

pragma(inline, true)
package auto max()(integer a, integer b)
{
    return a > b ? a : b;
}

@system nothrow @nogc pragma(inline, true)
package int symm_impl_(T)(
    ref const char side,
    ref const char uplo,
    ref const integer m,
    ref const integer  n,
    ref const T alpha,
        const(T)* a,
        ref const integer lda,
        const(T)* b,
        ref const integer ldb,
    ref const T beta,
        T* c,
        ref const integer ldc,
    Conjugated conj,
    )
{
    auto  _side = toUpper(side);
    auto  _uplo = toUpper(uplo);
    auto s = _side == 'L' ? Side.left : Side.right;
    auto u = _uplo == 'L' ? Uplo.lower : Uplo.upper;

    auto k = s == Side.left ? m : n;

    integer info = 0;

    if (llvm_expect(s != Side.left && _side != 'R', false))
        info = 1;
    else
    if (llvm_expect(u != Uplo.lower && _uplo != 'U', false))
        info = 2;
    else
    if (llvm_expect(m < 0, false))
        info = 3;
    else
    if (llvm_expect(n < 0, false))
        info = 4;
    else
    if (llvm_expect(lda < max(1, k), false))
        info = 7;
    else
    if (llvm_expect(ldb < max(1, m), false))
        info = 9;
    else
    if (llvm_expect(ldc < max(1, m), false))
        info = 12;
    else
    {
        import glas.internal.gemm: SL3;
        import glas.internal.symm: symm_impl;
        SL3!(T, T, T) arg = void;
        arg.asl = Slice!(2, const(T)*)([k, k], [1, lda], a);
        arg.bsl = Slice!(2, const(T)*)([m, n], [1, ldb], b);
        arg.csl = Slice!(2,       T *)([m, n], [1, ldc], c);
        arg.alpha_beta[0] = alpha;
        arg.alpha_beta[1] = beta;
        symm_impl(arg, s, u, conj, Conjugated.no);
        return 0;
    }
    enum nameSYMM = prefix!T ~ "SYMM";
    enum nameHEMM = prefix!T ~ "HEMM";
    //printf("info = %i\n", info);
    xerbla_(conj ? nameHEMM.ptr : nameSYMM.ptr, &info);
    return 0;
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
