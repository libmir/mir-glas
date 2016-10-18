module glas.precompiled.l3_;

version(LDC)
{
    version(unittest) {} else
    {
        pragma(LDC_no_moduleinfo);
    }
}

package enum string L3(string p, string T) =
q{
    version(LDC)
    {
        version(unittest) {} else
        {
            pragma(LDC_no_moduleinfo);
        }
    }

    import glas.common;
    import mir.ndslice.slice: Slice;
    import ldc.attributes: fastmath;

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
        import glas.internal.gemm: gemm_impl;
        gemm_impl(alpha, asl, bsl, beta, csl, conja, conjb);
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
        import glas.internal.symm: symm_impl;
        symm_impl(side, uplo, alpha, asl, bsl, beta, csl, conja, conjb);
    }
};
