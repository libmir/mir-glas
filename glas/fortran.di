module glas.fortran;

public import glas.common;

/++
CGEMM  performs one of the matrix-matrix operations
`C := alpha*op( A )*op( B ) + beta*C`,
where  `op( X )` is one of
`op( X ) = X`   or   `op( X ) = X**T`   or   `op( X ) = X**H`,
alpha and beta are scalars, and `A`, `B` and `C` are matrices, with `op( A )`
an m by k matrix,  `op( B )`  a  k by n matrix and  C an m by n matrix.
+/
int sgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const float alpha, const(float)* a, ref const FortranInt lda, const(float)* b, ref const FortranInt ldb, ref const float beta, float* c, ref const FortranInt ldc)
/// ditto
int dgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const double alpha, const(double)* a, ref const FortranInt lda, const(double)* b, ref const FortranInt ldb, ref const double beta, double* c, ref const FortranInt ldc)
/// ditto
int cgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc)
/// ditto
int zgemm_(ref const char transa, ref const char transb, ref const FortranInt m, ref const FortranInt  n, ref const FortranInt k, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc)
/// ditto
alias gemm_ = sgemm_;
/// ditto
alias gemm_ = dgemm_;
/// ditto
alias gemm_ = cgemm_;
/// ditto
alias gemm_ = zgemm_;

/++
SSYMM  performs one of the matrix-matrix operations
    `C := alpha*A*B + beta*C`,
or
    `C := alpha*B*A + beta*C`,
where alpha and beta are scalars,  A is a symmetric matrix and  B and
C are  m by n matrices.
+/
int ssymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const float alpha, const(float)* a, ref const FortranInt lda, const(float)* b, ref const FortranInt ldb, ref const float beta, float* c, ref const FortranInt ldc)
/// ditto
int dsymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const double alpha, const(double)* a, ref const FortranInt lda, const(double)* b, ref const FortranInt ldb, ref const double beta, double* c, ref const FortranInt ldc)
/// ditto
int csymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc)
/// ditto
int zsymm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc)
/// ditto
alias symm_ = ssymm_;
/// ditto
alias symm_ = dsymm_;
/// ditto
alias symm_ = csymm_;
/// ditto
alias symm_ = zsymm_;

/++
SSYMM  performs one of the matrix-matrix operations
    `C := alpha*A*B + beta*C`,
or
    `C := alpha*B*A + beta*C`,
where alpha and beta are scalars,  A is a hermitian matrix and  B and
C are  m by n matrices.
+/
int chemm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cfloat alpha, const(cfloat)* a, ref const FortranInt lda, const(cfloat)* b, ref const FortranInt ldb, ref const cfloat beta, cfloat* c, ref const FortranInt ldc)
/// ditto
int zhemm_(ref const char side, ref const char uplo, ref const FortranInt m, ref const FortranInt  n, ref const cdouble alpha, const(cdouble)* a, ref const FortranInt lda, const(cdouble)* b, ref const FortranInt ldb, ref const cdouble beta, cdouble* c, ref const FortranInt ldc)
/// ditto
alias hemm_ = chemm_;
/// ditto
alias hemm_ = zhemm_;
