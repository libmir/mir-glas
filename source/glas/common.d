/++
$(SCRIPT inhibitQuickIndex = 1;)

This is a submodule of $(MREF mir,glas).

Copyright: Ilya Yaroshenko 2016-.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas.common;
pragma(LDC_no_moduleinfo);

import ldc.attributes: fastmath;

@fastmath:

/++
Uplo specifies whether a matrix is an upper or lower triangular matrix.
+/
enum Uplo : bool
{
    /// upper triangular matrix.
    lower,
    /// lower triangular matrix
    upper,
}

/++
Convenient function to invert $(LREF Uplo) flag.
+/
Uplo swap()(Uplo type)
{
    return cast(Uplo) (type ^ 1);
}

///
unittest
{
    assert(swap(Uplo.upper) == Uplo.lower);
    assert(swap(Uplo.lower) == Uplo.upper);
}

/++
Diag specifies whether or not a matrix is unitriangular.
+/
enum Diag : bool
{
    /// a matrix assumed to be unit triangular
    unit,
    /// a matrix not assumed to be unit triangular
    nounit,
}

/++
On entry, `Side`  specifies whether  the  symmetric matrix  A
appears on the  left or right.
+/
enum Side : bool
{
    ///
    left,
    ///
    right,
}

/++
Convenient function to invert $(LREF Side) flag.
+/
Side swap()(Side type)
{
    return cast(Side) (type ^ 1);
}

///
enum Conjugated : bool
{
    ///
    no,
    ///
    yes,
}
