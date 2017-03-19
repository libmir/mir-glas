/++
Copyright: Copyright © 2016-, Ilya Yaroshenko.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas.internal.l1;

import mir.ndslice.slice: Slice, SliceKind;
import ldc.attributes: fastmath;
import glas.internal.utility;
import ldc.intrinsics: llvm_expect;
static import glas.ndslice;

@fastmath:

pragma(inline, true):

void scal(A, T)(A a, size_t len, sizediff_t str, T* ptr)
{
	if (str == 1)
	{
		static if (isComplex!T && isReal!A)
		{
			len *= 2;
			glas.ndslice.scal(a, len, str, cast(A*)ptr);
		}
		else
		{
			while (len)
			{
				*ptr *= a;
				ptr++;
				len--;
			}
		}
	}
	else
	if (str != 0)
	{
		while (len)
		{
			*ptr *= a;
			ptr += str;
			len--;
		}
	}
}

void copy(T)(size_t len, sizediff_t incx, const(T)* x, sizediff_t incy, T* y)
{
	if (incx == 1 && incy == 1)
	{
		static if (isComplex!T)
		{
			alias R = realType!T;
			len *= 2;
			glas.ndslice.copy(len, incx, cast(const(R)*)x, incy, cast(R*)y);
		}
		else
		{
			for(; len; len--)
				*y++ = *x++;
		}
		return;
	}
	for(; len; len--)
	{
		*y = *x;
		x += incx;
		y += incy;
	}
}

void swap(T)(size_t len, sizediff_t incx, T* x, sizediff_t incy, T* y)
{
	if (incx == 1 && incy == 1)
	{
		static if (isComplex!T)
		{
			alias R = realType!T;
			len *= 2;
			glas.ndslice.swap(len, incx, cast(R*)x, incy, cast(R*)y);
		}
		else
		{
			for(; len; len--)
			{
				auto t = *x;
				*x = *y;
				*y = t;
				x++;
				y++;
			}
		}
		return;
	}
	for(; len; len--)
	{
		auto t = *x;
		*x = *y;
		*y = t;
		x += incx;
		y += incy;
	}
}

void axpy(T)(T a, size_t len, sizediff_t incx, const(T)* x, sizediff_t incy, T* y)
{
	if (incx == 1 && incy == 1)
	{
		for(; len; len--)
			*y++ += a * *x++;
		return;
	}
	for(; len; len--)
	{
		*y += a * *x;
		x += incx;
		y += incy;
	}
}

float nrm2(T : float)(size_t len, sizediff_t incx, const(T)* x)
{
	import ldc.intrinsics: llvm_sqrt, llvm_fabs;
	double ret = 0;
	if (len == 1)
		return llvm_fabs(*x);
	if (incx == 1)
	{
		for(; len; len--)
		{
			auto t = double(*x);
			ret += t * t;
			x++;
		}
	}
	else
	{
		for(; len; len--)
		{
			auto t = double(*x);
			ret += t * t;
			x += incx;
		}
	}
	return llvm_sqrt(ret);
}

float nrm2(T : cfloat)(size_t len, sizediff_t incx, const(T)* x)
{
	import ldc.intrinsics: llvm_sqrt;
	if (incx == 1)
	{
		len *= 2;
		return glas.ndslice.nrm2(len, incx, cast(const(float)*) x);
	}
	double ret_re = 0;
	double ret_im = 0;
	for(; len; len--)
	{
		cdouble t = *x;
		ret_re += t.re * t.re;
		ret_im += t.im * t.im;
		x += incx;
	}
	return llvm_sqrt(ret_re + ret_im);
}

double nrm2(T : double)(size_t len, sizediff_t incx, const(T)* x)
{
	import ldc.intrinsics: llvm_sqrt, llvm_fabs;
	if (len == 1)
		return llvm_fabs(*x);
	double scale = 0;
	double ssq = 1;
	for(; len; len--)
	{
		if (auto v = *x)
		{
			auto absxi = llvm_fabs(v);
			if (scale < absxi)
			{
				auto t = scale / absxi;
				ssq = 1 + ssq * (t * t);
				scale = absxi;
			}
			else
			{
				auto t = absxi / scale;
				ssq += t * t;
			}
		}
		x += incx;
	}
	return scale * llvm_sqrt(ssq);
}

double nrm2(T : cdouble)(size_t len, sizediff_t incx, const(T)* x)
{
	import ldc.intrinsics: llvm_sqrt, llvm_fabs;
	double scale = 0;
	double ssq = 1;
	for(; len; len--)
	{
		if (auto v = x.re)
		{
			auto absxi = llvm_fabs(v);
			if (scale < absxi)
			{
				auto t = scale / absxi;
				ssq = 1 + ssq * (t * t);
				scale = absxi;
			}
			else
			{
				auto t = absxi / scale;
				ssq += t * t;
			}
		}
		if (auto v = x.im)
		{
			auto absxi = llvm_fabs(v);
			if (scale < absxi)
			{
				auto t = scale / absxi;
				ssq = 1 + ssq * (t * t);
				scale = absxi;
			}
			else
			{
				auto t = absxi / scale;
				ssq += t * t;
			}
		}
		x += incx;
	}
	return scale * llvm_sqrt(ssq);
}


realType!T asum(T)(size_t len, sizediff_t incx, const(T)* x)
	if (isComplex!T)
{
	import ldc.intrinsics: llvm_fabs;
	if (incx == 1)
	{
		len *= 2;
		return glas.ndslice.asum(len, incx, cast(const(typeof(return))*)x);
	}
	typeof(return) ret_re = 0;
	typeof(return) ret_im = 0;
	for(; len; len--)
	{
		ret_re += llvm_fabs(x.re);
		ret_im += llvm_fabs(x.im);
		x += incx;
	}
	return ret_re + ret_im;
}

T asum(T)(size_t len, sizediff_t incx, const(T)* x)
	if (!isComplex!T)
{
	import ldc.intrinsics: llvm_fabs;
	typeof(return) ret = 0;
	if (incx == 1)
	{
		for(; len; len--)
		{
			ret += llvm_fabs(*x);
			x++;
		}
	}
	for(; len; len--)
	{
		ret += llvm_fabs(*x);
		x += incx;
	}
	return ret;
}

ptrdiff_t iamax(T)(size_t len, sizediff_t incx, const(T)* x)
{
	import ldc.intrinsics: llvm_fabs;
	typeof(return) ret = -1;
	auto smax = -realType!T.infinity;
	if (incx == 1)
	{
		foreach(ptrdiff_t i; 0 .. len)
		{
			static if (!isComplex!T)
				auto r = llvm_fabs(x[i]);
			else
				auto r = llvm_fabs(x[i].re) + llvm_fabs(x[i].im);
			if (smax < r)
			{
				smax = r;
				ret = i;
			}
		}
	}
	else
	{
		foreach(ptrdiff_t i; 0 .. len)
		{
			static if (!isComplex!T)
				auto r = llvm_fabs(*x);
			else
				auto r = llvm_fabs(x.re) + llvm_fabs(x.im);
			if (smax < r)
			{
				smax = r;
				ret = i;
			}
			x += incx;
		}
	}
	return ret;
}
