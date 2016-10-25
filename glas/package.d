/++
Copyright: Ilya Yaroshenko 2016-.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas;
pragma(LDC_no_moduleinfo);

version(Have_mir)
	import mir.ndslice: Slice;
else
	import std.experimental.ndslice: Slice;

export extern(C) nothrow @nogc @system:

enum ulong ConjA = 0x1;
enum ulong ConjB = 0x2;
enum ulong Lower = 0x0;
enum ulong Left = 0x0;
enum ulong Upper = 0x0100;
enum ulong Right = 0x0200;

void glas_sscal(float a, Slice!(1, float*) xsl);
void glas_dscal(double a, Slice!(1, double*) xsl);
void glas_csscal(float a, Slice!(1, cfloat*) xsl);
void glas_cscal(cfloat a, Slice!(1, cfloat*) xsl);
void glas_csIscal(ifloat a, Slice!(1, cfloat*) xsl);
void glas_zdscal(double a, Slice!(1, cdouble*) xsl);
void glas_zscal(cdouble a, Slice!(1, cdouble*) xsl);
void glas_zdIscal(idouble a, Slice!(1, cdouble*) xsl);

void _glas_sscal(float a, size_t n, size_t incx, float* x);
void _glas_dscal(double a, size_t n, size_t incx, double* x);
void _glas_csscal(float a, size_t n, size_t incx, cfloat* x);
void _glas_cscal(cfloat a, size_t n, size_t incx, cfloat* x);
void _glas_csIscal(ifloat a, size_t n, size_t incx, cfloat* x);
void _glas_zdscal(double a, size_t n, size_t incx, cdouble* x);
void _glas_zscal(cdouble a, size_t n, size_t incx, cdouble* x);
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

void glas_sgemm(float alpha, Slice!(2, const(float)*) asl, Slice!(2, const(float)*) bsl, float beta, Slice!(2, float*) csl, ulong settings = 0);
void glas_dgemm(double alpha, Slice!(2, const(double)*) asl, Slice!(2, const(double)*) bsl, double beta, Slice!(2, double*) csl, ulong settings = 0);
void glas_cgemm(cfloat alpha, Slice!(2, const(cfloat)*) asl, Slice!(2, const(cfloat)*) bsl, cfloat beta, Slice!(2, cfloat*) csl, ulong settings = 0);
void glas_zgemm(cdouble alpha, Slice!(2, const(cdouble)*) asl, Slice!(2, const(cdouble)*) bsl, cdouble beta, Slice!(2, cdouble*) csl, ulong settings = 0);

alias gemm = glas_sgemm;
alias gemm = glas_dgemm;
alias gemm = glas_cgemm;
alias gemm = glas_zgemm;

void glas_ssymm(float alpha, Slice!(2, const(float)*) asl, Slice!(2, const(float)*) bsl, float beta, Slice!(2, float*) csl, ulong settings = 0);
void glas_dsymm(double alpha, Slice!(2, const(double)*) asl, Slice!(2, const(double)*) bsl, double beta, Slice!(2, double*) csl, ulong settings = 0);
void glas_csymm(cfloat alpha, Slice!(2, const(cfloat)*) asl, Slice!(2, const(cfloat)*) bsl, cfloat beta, Slice!(2, cfloat*) csl, ulong settings = 0);
void glas_zsymm(cdouble alpha, Slice!(2, const(cdouble)*) asl, Slice!(2, const(cdouble)*) bsl, cdouble beta, Slice!(2, cdouble*) csl, ulong settings = 0);

alias symm = glas_ssymm;
alias symm = glas_dsymm;
alias symm = glas_csymm;
alias symm = glas_zsymm;
