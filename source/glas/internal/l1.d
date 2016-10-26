/++
Copyright: Copyright © 2016-, Ilya Yaroshenko.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
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
