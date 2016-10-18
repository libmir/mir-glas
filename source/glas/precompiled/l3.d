module glas.precompiled.l3;

version(LDC)
{
    version(unittest) {} else
    {
        pragma(LDC_no_moduleinfo);
    }
}

import glas.precompiled.l3s;
import glas.precompiled.l3d;
import glas.precompiled.l3c;
import glas.precompiled.l3z;

import glas.common;
import mir.ndslice.slice: Slice;
import ldc.attributes: fastmath;

//extern(C) @system nothrow @nogc @fastmath pragma(inline, true):

alias glas_gemm = glas_sgemm;
alias glas_gemm = glas_dgemm;
alias glas_gemm = glas_cgemm;
alias glas_gemm = glas_zgemm;

alias glas_symm = glas_ssymm;
alias glas_symm = glas_dsymm;
alias glas_symm = glas_csymm;
alias glas_symm = glas_zsymm;

void foo(Slice!(2, float*) a, Slice!(2, float*) b, Slice!(2, float*) c)
{
	import glas.l3;
	gemm(1.0, a, b, 0.0, c);
}
