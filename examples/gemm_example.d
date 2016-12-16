#!/usr/bin/env dub
/+ dub.json:
{
    "name": "gemm_example",
    "targetType":"executable",
    "systemDependencies": "Example requires libmir-cpuid and libmir-glas",
    "lflags": ["-L$MIR_GLAS_PACKAGE_DIR", "-L$MIR_CPUID_PACKAGE_DIR", "-L.."],
    "dependencies": {
        "cblas": "~>1.0.0",
        "mir-glas":{
            "path": "../"
        },
        "mir-cpuid": "~>0.4.2"
    },
    "configurations": [
        {
            "name": "std"
        },
        {
            "name": "mir",
            "dependencies": {
                "mir": "~>0.22.0"
            }
        }
    ]
}
+/
import glas.ndslice;

alias T = double;

version(Have_mir)
{
    pragma(msg, "Mir based configuration");
    import mir.ndslice;
}
else
{
    pragma(msg, "Phobos based configuration");
    import std.experimental.ndslice;
    static if (__VERSION__ < 2072)
    {
        version = OldAPI;
        alias ConstMatrix = Slice!(2, const(T)*);
        auto slice(T,  size_t N)(size_t[N] shape...)
        {
            size_t len = 1;
            foreach(l; shape)
                len *= l;
            return new T[len].sliced(shape);
        }
    }
    else version = HaveImplicitConstCast;
}

int main()
{
    auto a = slice!T(3, 5);
    a[] =
        [[-5,  1,  7, 7, -4],
         [-1, -5,  6, 3, -3],
         [-5, -2, -3, 6,  0]];

    auto b = slice!T(5, 4);
    b[] =
        [[-5.0, -3,  3,  1],
         [ 4.0,  3,  6,  4],
         [-4.0, -2, -2,  2],
         [-1.0,  9,  4,  8],
         [  9.0, 8,  3, -2]];

    auto c = slice!T(3, 4);

    auto alpha = 1.0;
    auto beta  = 0.0;

    if(auto error_code = validate_gemm(a.structure, b.structure, c.structure))
    {
        import core.stdc.stdio;
        puts(glas_error(error_code).ptr);
        return 1;
    }

    version (HaveImplicitConstCast)
        gemm(alpha, a, b, beta, c);
    else
        gemm(alpha, cast(ConstMatrix)a, cast(ConstMatrix)b, beta, c);

    return c ==
        [[-42.0,  35,  -7, 77],
         [-69.0, -21, -42, 21],
         [ 23.0,  69,   3, 29]] ? 0 : -1;
}
