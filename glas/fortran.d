/++
Copyright: Ilya Yaroshenko 2016-.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas.fortran;

public import glas.common;

export extern(C) nothrow @nogc @system:

int sscal_(ref const FortranInt n, ref const float a, float* x, ref const FortranInt incx);
int dscal_(ref const FortranInt n, ref const double a, double* x, ref const FortranInt incx);
int csscal_(ref const FortranInt n, ref const float a, cfloat* x, ref const FortranInt incx);
int cscal_(ref const FortranInt n, ref const cfloat a, cfloat* x, ref const FortranInt incx);
int zdscal_(ref const FortranInt n, ref const double a, cdouble* x, ref const FortranInt incx);
int zscal_(ref const FortranInt n, ref const cdouble a, cdouble* x, ref const FortranInt incx);
alias scal = sscal_;
alias scal = dscal_;
alias scal = csscal_;
alias scal = cscal_;
alias scal = zdscal_;
alias scal = zscal_;

int sgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const float alpha, const(float)* a, ref const FortranInt lda, const(float)* b, ref const FortranInt ldb, ref const float beta, float* c, ref const FortranInt ldc);
int dgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const double alpha, const(double)* a, ref const FortranInt lda, const(double)* b, ref const FortranInt ldb, ref const double beta, double* c, ref const FortranInt ldc);
int cgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc);
int zgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc);
alias gemm_ = sgemm_;
alias gemm_ = dgemm_;
alias gemm_ = cgemm_;
alias gemm_ = zgemm_;

int ssymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const float alpha, const(float)* a, ref const FortranInt lda, const(float)* b, ref const FortranInt ldb, ref const float beta, float* c, ref const FortranInt ldc);
int dsymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const double alpha, const(double)* a, ref const FortranInt lda, const(double)* b, ref const FortranInt ldb, ref const double beta, double* c, ref const FortranInt ldc);
int csymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc);
int zsymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc);
alias symm_ = ssymm_;
alias symm_ = dsymm_;
alias symm_ = csymm_;
alias symm_ = zsymm_;

int chemm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc);
int zhemm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc);
alias hemm_ = chemm_;
alias hemm_ = zhemm_;
