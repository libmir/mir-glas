/++
Copyright: Ilya Yaroshenko 2016-.
License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: Ilya Yaroshenko
+/
module glas.precompiled.utility;
pragma(LDC_no_moduleinfo);

import glas.common;
import ldc.attributes: weak;

export extern(C) @system nothrow @nogc pragma(inline, false):

/* -- LAPACK auxiliary routine (preliminary version) -- */
/* Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/* Courant Institute, Argonne National Lab, and Rice University */
/* February 29, 1992 */
/* .. Scalar Arguments .. */
/* .. */
/* Purpose */
/* ======= */
/* XERBLA is an error handler for the LAPACK routines. */
/* It is called by an LAPACK routine if an input parameter has an */
/* invalid value. A message is printed and execution stops. */
/* Installers may consider modifying the STOP statement in order to */
/* call system-specific exception-handling facilities. */
/* Arguments */
/* ========= */
/* SRNAME (input) CHARACTER*6 */
/* The name of the routine which called XERBLA. */
/* INFO (input) INTEGER */
/* The position of the invalid parameter in the parameter list */
/* of the calling routine. */
@weak int xerbla_(in char* srname, ref FortranInt info)
{
    import core.stdc.stdio;
    printf("** On entry to %6s, parameter number %2i had an illegal value\n", srname, info);
    return 0;
}
