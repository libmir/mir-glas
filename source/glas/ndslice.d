/++
$(H2 GLAS API)

Copyright: Copyright © 2016-, Ilya Yaroshenko.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko

$(H4 Transposition)
GLAS does not require transposition parameters.
Use $(LINK2 http://dlang.org/phobos/std_experimental_ndslice_iteration.html#transposed, transposed)
to perform zero cost ndslice transposition.

Note: $(LINK2 , ndslice) uses is row major representation.
    $(BR)
+/
module glas.ndslice;

version(LDC)
    pragma(LDC_no_moduleinfo);

version(Have_mir)
{
    import mir.ndslice.slice: Slice;
}
else
{
    import std.experimental.ndslice.slice: Slice;
}

extern(C) nothrow @nogc @system:

/++
Specifies if the matrix `asl` stores conjugated elements.
+/
enum ulong ConjA = 0x1;
/++
Specifies if the matrix `bsl` stores conjugated elements.
+/
enum ulong ConjB = 0x2;
/++
Specifies if the lower  triangular
part of the symmetric matrix A is to be referenced.

The lower triangular part of the matrix `asl`
must contain the lower triangular part of the symmetric / hermitian
matrix A and the strictly upper triangular part of `asl` is not
referenced.

Note: Lower is default flag.
+/
enum ulong Lower = 0x0;
/++
Specifies if the upper  triangular
part of the symmetric matrix A is to be referenced.

The upper triangular
part of the matrix `asl`  must contain the upper triangular part
of the symmetric / hermitian matrix A and the strictly lower triangular
part of `asl` is not referenced.
+/
enum ulong Upper = 0x0100;
/++
Specifies if the symmetric/hermitian matrix A
    appears on the left in the  operation.

Note: Left is default flag.
+/
enum ulong Left = 0x0;
/++
Specifies if the symmetric/hermitian matrix A
    appears on the left in the  operation.
+/
enum ulong Right = 0x0200;

/++
Performs general matrix-matrix multiplication.

Pseudo_code: `C := alpha A × B + beta C`.

Params:
    alpha = scalar
    asl = `m ⨉ k` matrix
    bsl = `k ⨉ n` matrix
    beta = scalar. When  `beta`  is supplied as zero then the matrix `csl` need not be set on input.
    csl = `m ⨉ n` matrix with one stride equal to `±1`.
    settings = Operation settings. Allowed flags are $(LREF ConjA) and $(LREF ConjB).

Unified_alias: `gemm`

BLAS: SGEMM, DGEMM, CGEMM, ZGEMM
+/
void glas_sgemm(float alpha, Slice!(2, const(float)*) asl, Slice!(2, const(float)*) bsl, float beta, Slice!(2, float*) csl, ulong settings = 0);
/// ditto
void glas_dgemm(double alpha, Slice!(2, const(double)*) asl, Slice!(2, const(double)*) bsl, double beta, Slice!(2, double*) csl, ulong settings = 0);
/// ditto
void glas_cgemm(cfloat alpha, Slice!(2, const(cfloat)*) asl, Slice!(2, const(cfloat)*) bsl, cfloat beta, Slice!(2, cfloat*) csl, ulong settings = 0);
/// ditto
void glas_zgemm(cdouble alpha, Slice!(2, const(cdouble)*) asl, Slice!(2, const(cdouble)*) bsl, cdouble beta, Slice!(2, cdouble*) csl, ulong settings = 0);

alias gemm = glas_sgemm;
alias gemm = glas_dgemm;
alias gemm = glas_cgemm;
alias gemm = glas_zgemm;

/++
Performs symmetric or hermitian matrix-matrix multiplication.

Pseudo_code: `C := alpha A × B + beta C` or `C := alpha B × A + beta C`,
    where  `alpha` and `beta` are scalars, `A` is a symmetric or hermitian matrix and `B` and
    `C` are `m × n` matrices.

Params:
    alpha = scalar
    asl = `k ⨉ k` matrix, where `k` is `n`  when $(LREF Right) flag is set
           and is `m` otherwise.
    bsl = `m ⨉ n` matrix
    beta = scalar. When  `beta`  is supplied as zero then the matrix `csl` need not be set on input.
    csl = `m ⨉ n` matrix with one stride equals to `±1`.
    settings = Operation settings.
        Allowed flags are
            $(LREF Left), $(LREF Right),
            $(LREF Lower), $(LREF Upper),
            $(LREF ConjA), $(LREF ConjB).
            $(LREF ConjA) flag specifies if the matrix A is hermitian.

Unified_alias: `symm`

BLAS: SSYMM, DSYMM, CSYMM, ZSYMM, SHEMM, DHEMM, CHEMM, ZHEMM
+/
void glas_ssymm(float alpha, Slice!(2, const(float)*) asl, Slice!(2, const(float)*) bsl, float beta, Slice!(2, float*) csl, ulong settings = 0);
/// ditto
void glas_dsymm(double alpha, Slice!(2, const(double)*) asl, Slice!(2, const(double)*) bsl, double beta, Slice!(2, double*) csl, ulong settings = 0);
/// ditto
void glas_csymm(cfloat alpha, Slice!(2, const(cfloat)*) asl, Slice!(2, const(cfloat)*) bsl, cfloat beta, Slice!(2, cfloat*) csl, ulong settings = 0);
/// ditto
void glas_zsymm(cdouble alpha, Slice!(2, const(cdouble)*) asl, Slice!(2, const(cdouble)*) bsl, cdouble beta, Slice!(2, cdouble*) csl, ulong settings = 0);

alias symm = glas_ssymm;
alias symm = glas_dsymm;
alias symm = glas_csymm;
alias symm = glas_zsymm;

/++
`scal` scales a vector by a constant.

Pseudo_code: `x := a x`.

Unified_alias: `scal`

BLAS: SSCSCAL, DSCSCAL, CSCSCAL, ZSCSCAL, CSCAL, ZSCAL
+/
void glas_sscal(float a, Slice!(1, float*) xsl);
/// ditto
void glas_dscal(double a, Slice!(1, double*) xsl);
/// ditto
void glas_csscal(float a, Slice!(1, cfloat*) xsl);
/// ditto
void glas_cscal(cfloat a, Slice!(1, cfloat*) xsl);
/// ditto
void glas_csIscal(ifloat a, Slice!(1, cfloat*) xsl);
/// ditto
void glas_zdscal(double a, Slice!(1, cdouble*) xsl);
/// ditto
void glas_zscal(cdouble a, Slice!(1, cdouble*) xsl);
/// ditto
void glas_zdIscal(idouble a, Slice!(1, cdouble*) xsl);

/// ditto
void _glas_sscal(float a, size_t n, size_t incx, float* x);
/// ditto
void _glas_dscal(double a, size_t n, size_t incx, double* x);
/// ditto
void _glas_csscal(float a, size_t n, size_t incx, cfloat* x);
/// ditto
void _glas_cscal(cfloat a, size_t n, size_t incx, cfloat* x);
/// ditto
void _glas_csIscal(ifloat a, size_t n, size_t incx, cfloat* x);
/// ditto
void _glas_zdscal(double a, size_t n, size_t incx, cdouble* x);
/// ditto
void _glas_zscal(cdouble a, size_t n, size_t incx, cdouble* x);
/// ditto
void _glas_zdIscal(idouble a, size_t n, size_t incx, cdouble* x);

alias scal = glas_sscal;
alias scal = glas_dscal;
alias scal = glas_csscal;
alias scal = glas_cscal;
alias scal = glas_csIscal;
alias scal = glas_zdscal;
alias scal = glas_zscal;
alias scal = glas_zdIscal;

alias scal = _glas_sscal;
alias scal = _glas_dscal;
alias scal = _glas_csscal;
alias scal = _glas_cscal;
alias scal = _glas_csIscal;
alias scal = _glas_zdscal;
alias scal = _glas_zscal;
alias scal = _glas_zdIscal;
