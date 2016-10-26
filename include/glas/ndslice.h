#ifndef Have_GLAS
#define Have_GLAS
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

enum
{
	glas_ConjA = 0x1,
	glas_ConjB = 0x2,
	glas_Lower = 0x0,
	glas_Left = 0x0,
	glas_Upper = 0x0100,
	glas_Right = 0x0200,
};

struct glas_MatrixStructure
{
	size_t lengths[2];
	ptrdiff_t strides[2];
};

struct glas_MutVector
{
	size_t lengths[1];
	ptrdiff_t strides[1];
	void* ptr;
};

struct glas_ConstVector
{
	size_t lengths[1];
	ptrdiff_t strides[1];
	const void *ptr;
};

struct glas_MutMatrix
{
	size_t lengths[2];
	ptrdiff_t strides[2];
	void* ptr;
};

struct glas_ConstMatrix
{
	size_t lengths[2];
	ptrdiff_t strides[2];
	const void *ptr;
};

struct ErrorString
{
	size_t length;
	const char *ptr; // zero terminated
};

struct ErrorString glas_error(int error_code);
int glas_validate_gemm(struct glas_MatrixStructure as, struct glas_MatrixStructure bs, struct glas_MatrixStructure cs, unsigned long settings);
int glas_validate_symm(struct glas_MatrixStructure as, struct glas_MatrixStructure bs, struct glas_MatrixStructure cs, unsigned long settings);

void glas_sscal(float a, struct glas_MutVector xsl);
void glas_dscal(double a, struct glas_MutVector xsl);
void glas_csscal(float a, struct glas_MutVector xsl);
void glas_cscal(float _Complex a, struct glas_MutVector xsl);
void glas_csIscal(float a, struct glas_MutVector xsl);
void glas_zdscal(double a, struct glas_MutVector xsl);
void glas_zscal(double _Complex a, struct glas_MutVector xsl);
void glas_zdIscal(double a, struct glas_MutVector xsl);

void _glas_sscal(float a, size_t n, size_t incx, float* x);
void _glas_dscal(double a, size_t n, size_t incx, double* x);
void _glas_csscal(float a, size_t n, size_t incx, float _Complex* x);
void _glas_cscal(float _Complex a, size_t n, size_t incx, float _Complex* x);
void _glas_csIscal(float a, size_t n, size_t incx, float _Complex* x);
void _glas_zdscal(double a, size_t n, size_t incx, double _Complex* x);
void _glas_zscal(double _Complex a, size_t n, size_t incx, double _Complex* x);
void _glas_zdIscal(double a, size_t n, size_t incx, double _Complex* x);

void glas_sgemm(float alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, float beta, struct glas_MutMatrix csl, unsigned long settings);
void glas_dgemm(double alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, double beta, struct glas_MutMatrix csl, unsigned long settings);
void glas_cgemm(float _Complex alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, float _Complex beta, struct glas_MutMatrix csl, unsigned long settings);
void glas_zgemm(double _Complex alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, double _Complex beta, struct glas_MutMatrix csl, unsigned long settings);

void glas_ssymm(float alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, float beta, struct glas_MutMatrix csl, unsigned long settings);
void glas_dsymm(double alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, double beta, struct glas_MutMatrix csl, unsigned long settings);
void glas_csymm(float _Complex alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, float _Complex beta, struct glas_MutMatrix csl, unsigned long settings);
void glas_zsymm(double _Complex alpha, struct glas_ConstMatrix asl, struct glas_ConstMatrix bsl, double _Complex beta, struct glas_MutMatrix csl, unsigned long settings);

#ifdef __cplusplus
}
#endif
#endif
