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

@fastmath:

pragma(inline, true)
void scal(A, T)(A a, size_t len, sizediff_t str, T* ptr)
{
	if (str == 1)
	{
		static if (isComplex!T && isReal!A)
		{
			static import glas.ndslice;
			glas.ndslice.scal(a, len * 2, 1, cast(A*)ptr);
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

pragma(inline, true)
void copy(T)(size_t len, sizediff_t incx, const(T)* x, sizediff_t incy, T* y)
{
	if (incx == 1 && incy == 1)
	{
		for(; len; len--)
			*y++ = *x++;
		return;
	}
	for(; len; len--)
	{
		*y = *x;
		x += incx;
		y += incy;
	}
}

pragma(inline, true)
void swap(T)(size_t len, sizediff_t incx, T* x, sizediff_t incy, T* y)
{
	if (incx == 1 && incy == 1)
	{
		for(; len; len--)
		{
			auto t = *x;
			*x = *y;
			*y = t;
			x++;
			y++;
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

pragma(inline, true)
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
