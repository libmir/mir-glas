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

///
enum ulong ConjA = 0x1;
///
enum ulong ConjB = 0x2;
///
enum ulong Lower = 0x0;
///
enum ulong Left = 0x0;
///
enum ulong Upper = 0x0100;
///
enum ulong Right = 0x0200;
