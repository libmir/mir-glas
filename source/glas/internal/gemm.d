module glas.internal.gemm;

version(LDC)
{
    version(unittest) {} else
    {
        pragma(LDC_no_moduleinfo);
    }
}

import std.traits;
import std.meta;

public import glas.common;
import mir.ndslice.slice : Slice;
import mir.internal.utility;
import glas.internal.blocking;
import glas.internal.copy;
import glas.internal.config;

import ldc.attributes;
import ldc.intrinsics;
@fastmath:

version = PREFETCH;

template SL3(A, B, C)
{
    pragma(LDC_no_typeinfo)
    struct SL3
    {
        Slice!(2, const(A)*) asl = void;
        Slice!(2, const(B)*) bsl = void;
        Slice!(2, C*) csl = void;
        C[2] alpha_beta = void;
    }
}

pragma(inline, false)
void gemm_fast_path(A, B, C)(ref SL3!(A, B, C) abc)
{with(abc){
    mixin prefix3;
    mixin RegisterConfig!(PA, PB, PC, T);
    if (alpha_beta[1] == 0)
    {
        do {
            setZero(cast(T[])(csl.front!1.toDense)); // memset
            csl.popFront!1;
        }
        while (csl.length!1);
        return;
    }
    if (alpha_beta[1] == 1)
        return;
    do {
        scale(csl.front!1.toDense, alpha_beta[1]);
        csl.popFront!1;
    } 
    while (csl.length!1);
}}

pragma(inline, false)
@optStrategy("optsize")
nothrow @nogc
void gemm_impl(A, B, C)
(
    ref SL3!(A, B, C) abc,
    Conjugated conja,
    Conjugated conjb,
)
{with(abc){
    assert(asl.length!1 == bsl.length!0, "constraint: asl.length!1 == bsl.length!0");
    assert(csl.length!0 == asl.length!0, "constraint: csl.length!0 == asl.length!0");
    assert(csl.length!1 == bsl.length!1, "constraint: csl.length!1 == bsl.length!1");
    assert(csl.stride!0 == +1
        || csl.stride!0 == -1
        || csl.stride!1 == +1
        || csl.stride!1 == -1, "constraint: csl.stride!0 or csl.stride!1 must be equal to +/-1");

    mixin prefix3;
    mixin RegisterConfig!(PA, PB, PC, T);
    import mir.ndslice.iteration: reversed, transposed;
    //#########################################################
    if (llvm_expect(csl.anyEmpty, false))
        return;
    if (llvm_expect(csl.stride!0 < 0, false))
    {
        csl = csl.reversed!0;
        asl = asl.reversed!0;
    }
    if (llvm_expect(csl.stride!1 < 0, false))
    {
        csl = csl.reversed!1;
        bsl = bsl.reversed!1;
    }
    // change row based to column based
    if (csl.stride!0 != 1)
    {
        static if (is(A == B))
        {
            auto tsl = asl;
            asl = bsl.transposed;
            bsl = tsl.transposed;
            csl = csl.transposed;
            auto conjt = conja;
            conja = conjb;
            conjb = conjt;
        }
        else
        {
            SL3!(B, A, C) tr = void;
            tr.asl = bsl;
            tr.bsl = asl;
            tr.csl = csl;
            tr.alpha_beta[0] = alpha_beta[0];
            tr.alpha_beta[1] = alpha_beta[1];
            gemm_impl!(B, A, C)(tr, conjb, conja);
            return;
        }
    }
    assert(csl.stride!0 == 1);
    if (llvm_expect(asl.empty!1 || alpha_beta[0] == 0, false))
    {
        gemm_fast_path(abc);
        return;
    }
    //#########################################################
    PackKernel!(B, T)[mr_chain.length] pack_a_kernels = void;
    PackKernel!(B, T)[nr_chain.length] pack_b_kernels = void;
    Kernel!(PC, T)   [nr_chain.length]   beta_kernels = void;
    Kernel!(PC, T)   [nr_chain.length]    one_kernels = void;
    Kernel!(PC, T)*                           kernels = void;

    static if (PA == 2)
    {
        if (conja)
        foreach (mri, mr; mr_chain)
            pack_a_kernels[mri] = &pack_a_nano!(mr, PA, 1, A, T);
        else
        foreach (mri, mr; mr_chain)
            pack_a_kernels[mri] = &pack_a_nano!(mr, PA, 0, A, T);
    }
    else
    {
        foreach (mri, mr; mr_chain)
            pack_a_kernels[mri] = &pack_a_nano!(mr, PA, 0, A, T);
    }
    static if (PB == 2)
    {
        if (conjb)
        foreach (nri, nr; nr_chain)
            pack_b_kernels[nri] = &pack_b_nano!(nr, PB, 1, B, T);
        else
        foreach (nri, nr; nr_chain)
            pack_b_kernels[nri] = &pack_b_nano!(nr, PB, 0, B, T);
    }
    else
    {
        foreach (nri, nr; nr_chain)
            pack_b_kernels[nri] = &pack_b_nano!(nr, PB, 0, B, T);
    }
    foreach (nri, nr; nr_chain)
    {
        one_kernels[nri] = &gemv_reg!(BetaType.one, PA, PB, PC, nr, T);
    }
    kernels = one_kernels.ptr;
    if (alpha_beta[1] == 0)
    {
        foreach (nri, nr; nr_chain)
            beta_kernels[nri] = &gemv_reg!(BetaType.zero, PA, PB, PC, nr, T);
        kernels = beta_kernels.ptr;
    }
    else
    if (alpha_beta[1] != 1)
    {
        foreach (nri, nr; nr_chain)
            beta_kernels[nri] = &gemv_reg!(BetaType.beta, PA, PB, PC, nr, T);
        kernels = beta_kernels.ptr;
    }
    //#########################################################
    with(blocking!(PA, PB, PC, T)(asl.length!0, bsl.length!1, asl.length!1))
    {
        sizediff_t incb;
        if (mc < asl.length!0)
            incb = kc;
        do
        {
            if (asl.length!1 < kc)
                kc = asl.length!1;
            ////////////////////////
            auto aslp = asl[0 .. $, 0 .. kc];
            auto bsl_ptr = bsl.ptr;
            auto cslm = csl;
            auto mc = mc;
            //======================
            do
            {
                if (aslp.length!0 < mc)
                    mc = aslp.length!0;
                ////////////////////////
                pack_a!(A, T)(aslp[0 .. mc], a, pack_a_kernels.ptr, main_mr);
                //======================
                gebp!(PA, PB, PC, T, B)(
                    mc,
                    bsl.length!1,
                    kc,
                    a,
                    b,
                    incb,
                    bsl_ptr,
                    bsl.stride!0,
                    bsl.stride!1,
                    cast(T*) cslm.ptr,
                    cslm.stride!1,
                    *cast(T[PC][2]*)&alpha_beta,
                    pack_b_kernels.ptr,
                    kernels,
                    main_nr,
                    );
                ////////////////////////
                bsl_ptr = null;
                cslm.popFrontExactly!0(mc);
                aslp.popFrontExactly!0(mc);
            }
            while (aslp.length!0);
            ////////////////////////
            kernels = one_kernels.ptr;
            bsl.popFrontExactly!0(kc);
            asl.popFrontExactly!1(kc);
        }
        while (asl.length!1);
    }
}}

pragma(inline, false)
void setZero(T)(T[] a)
{
    assert(a.length);
    do
    {
        a[0] = 0;
        a = a[1 .. $];
    }
    while(a.length);
}

pragma(inline, false)
void scale(T)(T[] a, T c)
{
    assert(a.length);
    do
    {
        a[0] *= c;
        a = a[1 .. $];
    }
    while(a.length);
}

pragma(inline, true)
void gebp(size_t PA, size_t PB, size_t PC, T, B)(
    size_t mc,
    size_t nc,
    size_t kc,
    const(T)* a,
    T* b,
    sizediff_t incb,
    const(B)* ptrb,
    sizediff_t ldb,
    sizediff_t ldbe,
    T* c,
    sizediff_t ldc,
    ref const T[PC][2] alpha_beta,
    PackKernel!(B, T)* pack_b_kernels,
    Kernel!(PC, T)* kernels,
    size_t nr,
    )
{
    mixin RegisterConfig!(PA, PB, PC, T);
    do
    {
        if (nc >= nr) do
        {
            if (ptrb)
            {
                pack_b_kernels[0](kc, ldb, ldbe, ptrb, b);
                ptrb += nr * ldbe;
            }
            kernels[0](mc, kc, a, b, c, ldc, alpha_beta);
            b +=  nr * PB * incb;
            nc -= nr;
            c += nr * PC * ldc;
        }
        while (nc >= nr);
        //import core.bitop: bsr;
        pack_b_kernels++;
        kernels++;
        import ldc.intrinsics: llvm_ctlz;
        auto newNr = size_t(1) << (size_t.sizeof * 8 - 1 - llvm_ctlz(nr, true));
        nr = nr == newNr ? newNr / 2 : newNr;
    }
    while (nr);
}

alias Kernel(size_t P, T) =
    pure nothrow @nogc
    void function(
        size_t mc,
        size_t kc,
        const(T)* a,
        const(T)* b,
        T* c,
        sizediff_t ldc,
        ref const T[P][2] alpha_beta,
    );

pragma(inline, false)
pure nothrow @nogc
void gemv_reg (
    BetaType beta_type,
    size_t PA,
    size_t PB,
    size_t PC,
    size_t N,
    F,
)
(
    size_t mc,
    size_t kc,
    const(F)* a,
    const(F)* b,
    F* c,
    sizediff_t ldc,
    ref const F[PC][2] alpha_beta,
)
{
    mixin RegisterConfig!(PA, PB, PC, F);
    foreach (mri, mr; mr_chain)
    if (mc >= mr) do
    {
        enum M = Mi!(mri);
        alias V = Vi!(mri);
        auto as = a;
        a = cast(typeof(a)) dot_reg!beta_type(cast(const(V[M][PA])*)a, cast(const(F[PB][N])*)b, kc, cast(F[PC]*)c, ldc, alpha_beta);
        mc -= mr;
        c += mr * PC;
    }
    while (!mri && mc >= mr);
}

enum BetaType
{
    zero,
    one,
    beta,
}

pragma(inline, true)
pure nothrow @nogc
const(V[M][PA])* dot_reg (
    BetaType beta_type,
    size_t PA,
    size_t PB,
    size_t PC,
    size_t M,
    size_t N,
    V,
    F,
)
(
    const(V[M][PA])* a,
    const(F[PB][N])* b,
    size_t length,
    F[PC]* c,
    sizediff_t ldc,
    ref const F[PC][2] alpha_beta,
)
{
    prefetch_w!(V[M][PC].sizeof, N, 1)(c, ldc * c[0].sizeof);
    V[M][PC][N] reg = void;
    a = dot_reg_basic(a, b, length, reg);
    scale_nano(alpha_beta[0], reg);
    static if (beta_type == BetaType.zero)
        save_nano(reg, c, ldc);
    else
    static if (beta_type == BetaType.one)
        save_add_nano(reg, c, ldc);
    else
        save_madd_nano(reg, alpha_beta[1], c, ldc);
    return a;
}

pragma(inline, true)
//nothrow @nogc
const(V[M][PA])*
dot_reg_basic (
    size_t PA,
    size_t PB,
    size_t PC,
    size_t M,
    size_t N,
    V,
    F,
)
(
    const(V[M][PA])* a,
    const(F[PB][N])* b,
    size_t length,
    ref V[M][PC][N] c,
)
    if (is(V == F) || isSIMDVector!V)
{
    V[M][PC][N] reg = void;

    foreach (n; Iota!N)
    foreach (p; Iota!PC)
    foreach (m; Iota!M)
        reg[n][p][m] = 0;
    do
    {
        V[M][PA] ai = void;
        V[PB][N] bi = void;

        prefetch_r!(V[M][PA].sizeof, 1, 8, prefetchShift)(cast(void*)a, 0);

        foreach (p; Iota!PA)
        foreach (m; Iota!M)
            ai[p][m] = a[0][p][m];

        enum AB = PA + PB == 4;
        enum CA = PC + PA == 4;
        enum CB = PC + PB == 4;

        foreach (u; Iota!(N/2 + N%2))
        //foreach (u; Iota!(N))
        {
            alias um = Iota!(2*u, 2*u + 2 > N ? 2*u + 1 : 2*u + 2);
            //alias um = AliasSeq!(u);
            foreach (n; um)
            foreach (p; Iota!PB)
                bi[n][p] = b[0][n][p];
            foreach (n; um)
            foreach (m; Iota!M)
            {
                reg[n][0][m] += ai[0][m] * bi[n][0];
 static if (CB) reg[n][1][m] += ai[0][m] * bi[n][1];
 static if (AB) reg[n][0][m] -= ai[1][m] * bi[n][1];
 static if (CA) reg[n][1][m] += ai[1][m] * bi[n][0];
            }
        }
        a++;
        b++;
    }
    while (--length);
    load_nano(c, reg);
    return a;
}

pragma(LDC_intrinsic, "llvm.prefetch")
    pure nothrow @nogc
    void llvm_prefetch(void* ptr, uint rw, uint locality, uint cachetype);

pragma(inline, true)
pure
void prefetch_w(size_t M, size_t N, size_t rem = 1)(void* ptr, sizediff_t ld)
{
    version(PREFETCH)
    {
        foreach (n; Iota!N)
        {
            foreach (m; Iota!(M / 64 + bool(M % 64 >= rem)))
                llvm_prefetch(ptr + m * 64, 1, 3, 1);
            ptr += ld;
        }
    }
}

pragma(inline, true)
pure
void prefetch_r(size_t M, size_t N, size_t rem, size_t shift)(void* ptr, sizediff_t ld)
{
    version(PREFETCH)
    {
        foreach (n; Iota!N)
        {
            foreach (m; Iota!(M / 64 + bool(M % 64 >= rem)))
            {
                llvm_prefetch(ptr + m * 64 + shift + ld * n, 0, 3, 1);
            }
        }
    }
}

pragma(inline, true)
void scale_nano(size_t M, size_t P, size_t N, V, F)(ref const F[P] alpha, ref V[M][P][N] c)
{
    V[P] s = void;
    V[M][P][N] reg = void;
    load_nano(s, alpha);
    load_nano(reg, c);
    foreach (n; Iota!N)
    foreach (m; Iota!M)
    {
        static if (P == 1)
        {
            reg[n][0][m] *= s[0];
        }
        else
        {
            auto re = s[0] * reg[n][0][m];
            auto im = s[0] * reg[n][1][m];
            re -= s[1] * reg[n][1][m];
            im += s[1] * reg[n][0][m];
            reg[n][0][m] = re;
            reg[n][1][m] = im;
        }
    }
    load_nano(c, reg);
}

pragma(inline, true)
ref auto castByRef(C)(return ref C val)
{
    static if (isComplex!C)
        alias R = typeof(val.re)[2];
    else
        alias R = C[1];
    return *cast(R*) &val;
}
