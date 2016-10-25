/++
$(H2 BLAS API for GLAS)

Please read $(LINK2 http://www.netlib.org/blas/ , Netlib BLAS) for more details.

Note: Standard (fortran) BLAS API is column major.

Copyright: Ilya Yaroshenko 2016-.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas.fortran;

version(LDC)
	pragma(LDC_no_moduleinfo);

/// Alias for Fortran77 integer.
alias FortranInt = int;

export extern(C) nothrow @nogc @system:

/++
 `gemm_`  performs one of the matrix-matrix operations

    `C := alpha*op( A )*op( B ) + beta*C`,

 where  `op( X )` is one of

    `op( X ) = X   or   op( X ) = X**T or op( X ) = X**H`,

 alpha and beta are scalars, and `A`, `B` and `C` are matrices, with `op( A )`
 an `m ⨉ k` matrix,  `op( B )`  a  `k ⨉ n` matrix and  `C` an `m ⨉ n` matrix.

Unified_alias: `gemm_`
+/
int sgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const float alpha, const(float)* a, ref const FortranInt lda, const(float)* b, ref const FortranInt ldb, ref const float beta, float* c, ref const FortranInt ldc);
/// ditto
int dgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const double alpha, const(double)* a, ref const FortranInt lda, const(double)* b, ref const FortranInt ldb, ref const double beta, double* c, ref const FortranInt ldc);
/// ditto
int cgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc);
/// ditto
int zgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc);

alias gemm_ = sgemm_;
alias gemm_ = dgemm_;
alias gemm_ = cgemm_;
alias gemm_ = zgemm_;

/++
 `symm_`  performs one of the matrix-matrix operations

    `C := alpha*A*B + beta*C`,

 where  `op( X )` is one of

    `C := alpha*B*A + beta*C`,

 alpha and beta are scalars, `A` is a symmetric matrix and `B` and
 `C` are `m ⨉ n` matrices..

Unified_alias: `symm_`
+/
int ssymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const float alpha, const(float)* a, ref const FortranInt lda, const(float)* b, ref const FortranInt ldb, ref const float beta, float* c, ref const FortranInt ldc);
/// ditto
int dsymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const double alpha, const(double)* a, ref const FortranInt lda, const(double)* b, ref const FortranInt ldb, ref const double beta, double* c, ref const FortranInt ldc);
/// ditto
int csymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc);
/// ditto
int zsymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc);

alias symm_ = ssymm_;
alias symm_ = dsymm_;
alias symm_ = csymm_;
alias symm_ = zsymm_;

/++
 `hemm_`  performs one of the matrix-matrix operations

    `C := alpha*A*B + beta*C`,

 where  `op( X )` is one of

    `C := alpha*B*A + beta*C`,

 alpha and beta are scalars, `A` is a hermitian matrix and `B` and
 `C` are `m ⨉ n` matrices..

Unified_alias: `hemm_`
+/
int chemm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc);
/// ditto
int zhemm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc);

alias hemm_ = chemm_;
alias hemm_ = zhemm_;

/++
`scal_` scales a vector by a constant.

Unified_alias: `scal_`
+/
int sscal_(ref const FortranInt n, ref const float a, float* x, ref const FortranInt incx);
/// ditto 
int dscal_(ref const FortranInt n, ref const double a, double* x, ref const FortranInt incx);
/// ditto 
int csscal_(ref const FortranInt n, ref const float a, cfloat* x, ref const FortranInt incx);
/// ditto 
int cscal_(ref const FortranInt n, ref const cfloat a, cfloat* x, ref const FortranInt incx);
/// ditto 
int csIscal_(ref const FortranInt n, ref const ifloat a, cfloat* x, ref const FortranInt incx);
/// ditto 
int zdscal_(ref const FortranInt n, ref const double a, cdouble* x, ref const FortranInt incx);
/// ditto 
int zscal_(ref const FortranInt n, ref const cdouble a, cdouble* x, ref const FortranInt incx);
/// ditto 
int zdIscal_(ref const FortranInt n, ref const idouble a, cdouble* x, ref const FortranInt incx);

alias scal_ = sscal_;
alias scal_ = dscal_;
alias scal_ = csscal_;
alias scal_ = cscal_;
alias scal_ = csIscal_;
alias scal_ = zdscal_;
alias scal_ = zscal_;
alias scal_ = zdIscal_;

/++
 XERBLA  is an error handler for the LAPACK routines.
 It is called by an LAPACK routine if an input parameter has an
 invalid value.  A message is printed and execution stops.

 Installers may consider modifying the STOP statement in order to
 call system-specific exception-handling facilities.
+/
int xerbla_(in char* srname, ref FortranInt info);
