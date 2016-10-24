module glas.internal.l1;

import std.experimental.ndslice.slice: Slice;
import ldc.attributes: fastmath;
import glas.internal.utility;
import ldc.intrinsics: llvm_expect;

pragma(inline, true)
@fastmath void scal(A, T)(A a, size_t len, sizediff_t str, T* ptr)
{
	if (str == 1)
	{
		static if (isComplex!T && isReal!A)
		{
			static import glas;
			auto ysl = _toSlice!(1, A*)([len * 2], [1], cast(A*)ptr);
			glas.scal(a, ysl);
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
