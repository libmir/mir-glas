/++
Copyright: Ilya Yaroshenko 2016-.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas.internal.symm;
pragma(LDC_no_moduleinfo);

import std.traits;
import std.meta;
import std.experimental.ndslice.slice : Slice;

import ldc.attributes;
import ldc.intrinsics;

import glas.common;
import glas.internal.utility;
import glas.internal.blocking;
import glas.internal.copy;
import glas.internal.config;
import glas.internal.gemm;


@fastmath:

pragma(inline, false)
@optStrategy("optsize")
nothrow @nogc
void symm_impl(A, B, C)
(
    ref SL3!(A, B, C) abc,
    Side side,
    Uplo uplo,
    Conjugated conja,
    Conjugated conjb,
)
{with(abc){
    assert(asl.length!0 == asl.length!1, "constraint: asl.length!0 == asl.length!1");
    assert(asl.length!1 == bsl.length!0, "constraint: asl.length!1 == bsl.length!0");
    assert(csl.length!0 == asl.length!0, "constraint: csl.length!0 == asl.length!0");
    assert(csl.length!1 == bsl.length!1, "constraint: csl.length!1 == bsl.length!1");
    assert(csl.stride!0 == +1
        || csl.stride!0 == -1
        || csl.stride!1 == +1
        || csl.stride!1 == -1, "constraint: csl.stride!0 or csl.stride!1 must be equal to +/-1");

    mixin prefix3;
    import std.experimental.ndslice.iteration: reversed, transposed;
    mixin RegisterConfig!(PA, PB, PC, T);
    //#########################################################
    if (llvm_expect(csl.empty!0 || csl.empty!1, false))
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
    if (csl.stride!0 != 1)
    {
        asl = asl.transposed;
        bsl = bsl.transposed;
        csl = csl.transposed;
        uplo = swap(uplo);
        side = swap(side);
    }
    static if (PA == 2)
    {
        int hem = conja;
    }
    if (uplo ^ side)
    {
        asl = asl.transposed;
        static if (PA == 2)
        {
            hem = -hem;
        }
    }
    assert(csl.stride!0 == 1);
    if (llvm_expect(asl.empty!1 || alpha_beta[0] == 0, false))
    {
        gemm_fast_path(abc);
        return;
    }
    //#########################################################
    Kernel!(PC, T)[nr_chain.length] beta_kernels = void;
    Kernel!(PC, T)[nr_chain.length] one_kernels  = void;
    Kernel!(PC, T)* kernels                      = void;

    static if (PA == PB)
    {
        foreach (nri, nr; nr_chain)
            one_kernels [nri] = &gemv_reg!(BetaType.one, PA, PB, PC, nr, T);

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
    }
    //#########################################################
    if (!side)
    {
        //#########################################################
        PackKernel!(B, T)[PA][mr_chain.length] pack_a_kernels = void;
        PackKernel!(B, T)    [nr_chain.length] pack_b_kernels = void;
        PackKernelTri!(A, T) pack_a_tri_kernel = void;

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
        static if (PA != PB)
        {
            foreach (nri, nr; nr_chain)
                one_kernels [nri] = &gemv_reg!(BetaType.one, PA, PB, PC, nr, T);

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
        }
        static if (PA == 2)
        {
            if (hem == 0)
            {
                foreach (mri, mr; mr_chain)
                {
                    pack_a_kernels[mri][0] = &pack_a_nano!(mr, PA, 0, A, T);
                    pack_a_kernels[mri][1] = &pack_a_nano!(mr, PA, 0, A, T);
                }
                pack_a_tri_kernel = &pack_a_tri!(PA, A, T, 0);
            }
            else
            if (hem < 0)
            {
                foreach (mri, mr; mr_chain)
                {
                    pack_a_kernels[mri][0] = &pack_a_nano!(mr, PA, 1, A, T);
                    pack_a_kernels[mri][1] = &pack_a_nano!(mr, PA, 0, A, T);
                }
                pack_a_tri_kernel = &pack_a_tri!(PA, A, T, -1);
            }
            else
            {
                foreach (mri, mr; mr_chain)
                {
                    pack_a_kernels[mri][0] = &pack_a_nano!(mr, PA, 0, A, T);
                    pack_a_kernels[mri][1] = &pack_a_nano!(mr, PA, 1, A, T);
                }
                pack_a_tri_kernel = &pack_a_tri!(PA, A, T, +1);
            }
        }
        else
        {
            foreach (mri, mr; mr_chain)
                pack_a_kernels[mri][0] = &pack_a_nano!(mr, PA, 0, A, T);
            pack_a_tri_kernel = &pack_a_tri!(PA, A, T, 0);
        }
        //#########################################################
        with(blocking!(PA, PB, PC, T)(asl.length!0, bsl.length!1, asl.length!0))
        {
            size_t j;
            sizediff_t incb;
            if (mc  < asl.length!0)
                incb = kc;
            do
            {
                if (asl.length!0 - j < kc)
                    kc = asl.length!0 - j;
                ////////////////////////
                size_t i;
                auto bsl_ptr = bsl.ptr;
                auto cslm = csl;
                auto mc = mc;
                //======================
                do
                {
                    if (asl.length!0 - i < mc)
                        mc = asl.length!0 - i;
                    pack_a_sym!(PA, A, T)(asl.ptr, asl.stride!0, asl.stride!1, i, j, mc, kc, a, pack_a_kernels.ptr, pack_a_tri_kernel, main_mr);
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
                    i += mc;
                }
                while (i < asl.length!0);
                ////////////////////////
                kernels = one_kernels.ptr;
                bsl.popFrontExactly!0(kc);
                j += kc;
            }
            while (j < asl.length!0);
        }
    }
    else
    {
        //#########################################################
        PackKernel!(B, T)    [mr_chain.length] pack_a_kernels = void;
        PackKernel!(A, T)[PB][nr_chain.length] pack_b_kernels = void;
        PackKernelTri!(A, T) pack_b_tri_kernel = void;
        static if (PB == 2)
        {
            if (conjb)
            foreach (mri, mr; mr_chain)
                pack_a_kernels[mri] = &pack_a_nano!(mr, PB, 1, B, T);
            else
            foreach (mri, mr; mr_chain)
                pack_a_kernels[mri] = &pack_a_nano!(mr, PB, 0, B, T);
        }
        else
        {
            foreach (mri, mr; mr_chain)
                pack_a_kernels[mri] = &pack_a_nano!(mr, PB, 0, B, T);
        }
        static if (PA == 2)
        {
            if (hem == 0)
            {
                foreach (nri, nr; nr_chain)
                {
                    pack_b_kernels[nri][0] = &pack_b_nano!(nr, PA, 0, A, T);
                    pack_b_kernels[nri][1] = &pack_b_nano!(nr, PA, 0, A, T);
                }
                pack_b_tri_kernel = &pack_b_tri!(PA, A, T, 0);
            }
            else
            if (hem < 0)
            {
                foreach (nri, nr; nr_chain)
                {
                    pack_b_kernels[nri][0] = &pack_b_nano!(nr, PA, 1, A, T);
                    pack_b_kernels[nri][1] = &pack_b_nano!(nr, PA, 0, A, T);
                }
                pack_b_tri_kernel = &pack_b_tri!(PA, A, T, -1);
            }
            else
            {
                foreach (nri, nr; nr_chain)
                {
                    pack_b_kernels[nri][0] = &pack_b_nano!(nr, PA, 0, A, T);
                    pack_b_kernels[nri][1] = &pack_b_nano!(nr, PA, 1, A, T);
                }
                pack_b_tri_kernel = &pack_b_tri!(PA, A, T, +1);
            }
        }
        else
        {
            foreach (nri, nr; nr_chain)
                pack_b_kernels[nri][0] = &pack_b_nano!(nr, PA, 0, A, T);
            pack_b_tri_kernel = &pack_b_tri!(PA, A, T, 0);
        }
        static if (PA != PB)
        {
            foreach (nri, nr; nr_chain)
                one_kernels [nri] = &gemv_reg!(BetaType.one, PB, PA, PC, nr, T);

            kernels = one_kernels.ptr;
            if (beta == 0)
            {
                foreach (nri, nr; nr_chain)
                    beta_kernels[nri] = &gemv_reg!(BetaType.zero, PB, PA, PC, nr, T);
                kernels = beta_kernels.ptr;
            }
            else
            if (beta != 1)
            {
                foreach (nri, nr; nr_chain)
                    beta_kernels[nri] = &gemv_reg!(BetaType.beta, PB, PA, PC, nr, T);
                kernels = beta_kernels.ptr;
            }
        }
        //#########################################################
        with(blocking!(PB, PA, PC, T)(bsl.length!0, asl.length!0, bsl.length!1))
        {
            sizediff_t incb;
            if (mc < bsl.length!0)
                incb = kc;
            size_t j;
            do
            {
                if (bsl.length!1 < kc)
                    kc = bsl.length!1;
                ////////////////////////
                auto bslp = bsl[0 .. $, 0 .. kc];
                auto asl_ptr = asl.ptr;
                auto cslm = csl;
                auto mc = mc;
                //======================
                do
                {
                    if (bslp.length!0 < mc)
                        mc = bslp.length!0;
                    ////////////////////////
                    pack_a!(B, T)(bslp[0 .. mc], a, pack_a_kernels.ptr, main_mr);
                    //======================
                    sybp!(PB, PA, PC, T, A)(
                        mc,
                        asl.length!0,
                        kc,
                        a,
                        b,
                        incb,
                        asl_ptr,
                        asl.stride!0,
                        asl.stride!1,
                        j,
                        cast(T*) cslm.ptr,
                        cslm.stride!1,
                        *cast(T[PC][2]*)&alpha_beta,
                        pack_b_tri_kernel,
                        pack_b_kernels.ptr,
                        kernels,
                        main_nr,
                        );
                    ////////////////////////
                    asl_ptr = null;
                    cslm.popFrontExactly!0(mc);
                    bslp.popFrontExactly!0(mc);
                }
                while (bslp.length!0);
                ////////////////////////
                j += kc;
                kernels = one_kernels.ptr;
                bsl.popFrontExactly!1(kc);
            }
            while (bsl.length!1);
        }
    }
}}

pragma(inline, true)
void sybp(size_t PA, size_t PB, size_t PC, T, B)(
    size_t mc,
    size_t nc,
    size_t kc,
    scope const(T)* a,
    scope T* b,
    sizediff_t incb,
    const(B)* ptrb,
    sizediff_t str0b,
    sizediff_t str1b,
    size_t js,
    scope T* c,
    sizediff_t ldc,
    ref const T[PC][2] alpha_beta,
    PackKernelTri!(B, T) pack_b_tri_kernel,
    PackKernel!(B, T)[PB]* pack_b_kernels,
    Kernel!(PC, T)* kernels,
    size_t nr,
    )
{
    mixin RegisterConfig!(PA, PB, PC, T);
    size_t i;
    do
    {
        if (nc >= nr) do
        {
            if (ptrb)
            {
                size_t j = js;
                size_t length = kc;
                auto to = b;
                {
                    sizediff_t len = i - j;
                    static if (PB == 1)
                        len++;
                    if (len > 0)
                    {
                        if (len > length)
                            len = length;
                        to = pack_b_kernels[0][0](len, str0b, str1b, ptrb + j * str0b + i * str1b, to);
                        j += len;
                        length -= len;
                    }
                }
                {
                    sizediff_t start = j - i;
                    sizediff_t len = nr - start;
                    static if (PB == 1)
                        len--;
                    if (len > length)
                        len = length;
                    if (len > 0)
                    {
                        to = pack_b_tri_kernel(ptrb + i * str0b + j * str1b, str0b, str1b, to, start, len, nr);
                        length -= len;
                        j += len;
                    }
                }
                if(length)
                {
                    pack_b_kernels[0][$-1](length, str1b, str0b, ptrb + i * str0b + j * str1b, to);
                }
                i += nr;
            }
            kernels[0](mc, kc, a, b, c, ldc, alpha_beta);
            b +=  nr * PB * incb;
            nc -= nr;
            c += nr * PC * ldc;
        }
        while (nc >= nr);
        pack_b_kernels++;
        kernels++;
        import ldc.intrinsics: llvm_ctlz;
        auto newNr = size_t(1) << (size_t.sizeof * 8 - 1 - llvm_ctlz(nr, true));
        nr = nr == newNr ? newNr / 2 : newNr;
    }
    while(nr);
}
