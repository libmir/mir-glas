module glas.context;

version(LDC)
{
    version(unittest) {} else
    {
        pragma(LDC_no_moduleinfo);
    }
}

import glas.internal.memory;

package __gshared uint c1;
package __gshared uint c2;
package __gshared uint line;
package __gshared uint initialized;
package __gshared void[] _memory;


import ldc.attributes : fastmath;
@fastmath:

/++
Initialize GLAS Context. Optional.
+/
pragma(inline, false)
export extern(C) nothrow @nogc void glas_init()
{
    import std.range.primitives;
    import cpuid.unified;

    if(initialized)
        return;
    cpuid_init();
    auto dc = dCache;
    auto uc = uCache;

    while (!uc.empty && uc.back.size > (1024 * 64)) // > 64 MB is CPU memory
    {
        uc.popBack;
    }

    if (dc.length)
    {
        c1 = dc.front.size;
        line = dc.front.line;
        dc.popFront;
    }
    else
    if (uc.length)
    {
        c1 = uc.front.size;
        line = uc.front.line;
        uc.popFront;
    }
    else
    {
        c1 = 16;
    }

    if (uc.length)
    {
        c2 = uc.back.size;
    }
    else
    if (dc.length)
    {
        c2 = dc.back.size;
    }
    else
    {
        c1 = 256;
    }

    c1 <<= 10;
    c2 <<= 10;
    if (line == 0)
        line = 64;
    initialized = true;
}

/++
Releases memory and closes threads.
Optional.
+/
export extern(C) nothrow @nogc void glas_release()
{
    if (_memory !is null)
        deallocate(_memory);
}

// Returns: reused unaligned memory chunk
package nothrow @nogc void[] memory(size_t size)
{
    if (_memory.length < size)
    {
        auto f = _memory.length << 1;
        if (f > size)
            size = f;
        if (_memory !is null)
            deallocate(_memory);
        _memory = alignedAllocate(size, 4096);
    }
    return _memory[0 .. size];
}