module glas.precompiled.l3;

version(LDC)
{
    version(unittest) {} else
    {
        pragma(LDC_no_moduleinfo);
    }
}

public import glas.precompiled.l3s;
public import glas.precompiled.l3d;
public import glas.precompiled.l3c;
public import glas.precompiled.l3z;

alias glas_gemm = glas_sgemm;
alias glas_gemm = glas_dgemm;
alias glas_gemm = glas_cgemm;
alias glas_gemm = glas_zgemm;

alias glas_symm = glas_ssymm;
alias glas_symm = glas_dsymm;
alias glas_symm = glas_csymm; // includes hemm
alias glas_symm = glas_zsymm; // includes hemm

alias gemm_ = sgemm_;
alias gemm_ = dgemm_;
alias gemm_ = cgemm_;
alias gemm_ = zgemm_;

alias symm_ = ssymm_;
alias symm_ = dsymm_;
alias symm_ = csymm_;
alias symm_ = zsymm_;

alias hemm_ = chemm_;
alias hemm_ = zhemm_;
