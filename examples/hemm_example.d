#!/usr/bin/env dub
/+ dub.json:
{
    "name": "hemm_example",
    "targetType":"executable",
    "systemDependencies": "Example requires libmir-cpuid and libmir-glas",
    "libs": ["mir-cpuid", "mir-glas"],
    "lflags": ["-L../"],
    "dependencies": {
        "mir-glas":{
            "path": "../"
        }
    },
    "configurations": [
        {
            "name": "std"
        },
        {
            "name": "mir",
            "dependencies": {
                "mir": "~>0.20.2"
            }
        }
    ],
}
+/
import glas.ndslice;

alias T = cdouble;

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
}

int main()
{
    auto a = slice!T(3, 3);
    a[] =
        [[-2 + 0i, T.init, T.init],
         [+3 + 2i, -5 + 0i, T.init],
         [-4 + 7i, -2 + 3i, -3 + 0i]];

    auto b = slice!T(3, 4);
    b[] =
        [[-5 + 3i, -3 + 9i,  3 + 2i, 1 + 2i],
         [ 4 + 5i,  3 + 4i,  6 + 5i, 4 + 9i],
         [-4 + 2i, -2 + 2i, -2 + 7i, 2 + 6i]];

    auto c = slice!T(3, 4);
    auto d = slice!T(3, 4);
    auto alpha = 1 + 0i;
    auto beta  = 0 + 0i;

    version(OldAPI)
        symm(alpha, cast(ConstMatrix)a, cast(ConstMatrix)b, beta, c, ConjA | Left | Lower);
    else
        symm(alpha, a, b, beta, c, ConjA | Left | Lower);

    foreach(i; 0 .. a.length)
    foreach(j; i+1 .. a.length)
        a[i, j] = a[j, i].re - a[j, i].im * 1fi;

    version(OldAPI)
        gemm(alpha, cast(ConstMatrix)a, cast(ConstMatrix)b, beta, d);
    else
        gemm(alpha, a, b, beta, d);

    return c == d ? 0 : -1;
}
