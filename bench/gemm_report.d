#!/usr/bin/env dub
/+ dub.json:
{
	"name": "gemm_report",
	"dependencies": {"mir-glas": {"path": ".."}, },
	"libs": ["blas"],
    "lflags": ["-L$MIR_GLAS_PACKAGE_DIR", "-L.."],
    "dependencies": {
		"cblas": "~>1.0.0",
        "mir-glas":{
            "path": "../"
        }
        "mir-cpuid": "~>0.4.2",
    },
	"dflags-ldc": ["-mcpu=native"],
}
+/
	//"lflags": ["-L/opt/intel/mkl/lib"],
	//"libs": ["mkl_sequential", "mkl_core", "mkl_intel_lp64"],

// Set up your libblas to approporiate version, or just copy it to the benchmarks/glas folder.
// Note: GLAS is single thread for now.
// $ dub build --compiler=ldmd2 -b release --single gemm_report.d
// $ ./gemm_report
import std.math;
import std.traits;
import std.datetime;
import std.conv;
import std.algorithm.comparison;
import std.stdio;
import std.exception;
import std.getopt;
import std.experimental.ndslice;
import glas.ndslice;

alias C = float;
//alias C = double;
//alias C = cfloat;
//alias C = cdouble;

size_t[] reportValues = [
	10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
	200, 300, 500, 600, 700, 800, 900, 1000,
	1200, 1400, 1600, 1800, 2000];

void main(string[] args)
{
	size_t count = 10;
	auto helpInformation = 
	getopt(args,
		"count|c", "Iteration count. Default value is " ~ count.to!string, &count);
	if (helpInformation.helpWanted)
	{
		defaultGetoptPrinter("Parameters:", helpInformation.options);
		return;
	}

	writeln("m=n=k,GLAS(thread_count=1),BLAS(thread_count=?)");
	foreach(m; reportValues)
	{
		auto n = m;
		auto k = m;

		auto d = new C[m * n].sliced(m, n);
		auto c = new C[m * n].sliced(m, n);
		auto a = new C[m * k].sliced(m, k);
		auto b = new C[k * n].sliced(k, n);

		fillRNG(c);
		fillRNG(a);
		fillRNG(b);

		d[] = c[];

		static if(isComplex!C)
		{
			C alpha = 3 + 7i;
			C beta = 0 + 0i;
		}
		else
		{
			C alpha = 3;
			C beta = 0;
		}

		/// BLAS
		auto nsecsBLAS = double.max;
		foreach(_; 0..count) {
			StopWatch sw;
			sw.start;
			static if(!(is(C == real) || is(C == creal) || is(C : long)))
			{
				static import cblas;
				static if(isComplex!C)
				cblas.gemm(
					cblas.Order.RowMajor,
					cblas.Transpose.NoTrans,
					cblas.Transpose.NoTrans,
					cast(cblas.blasint) m,
					cast(cblas.blasint) n,
					cast(cblas.blasint) k,
					& alpha,
					a.ptr,
					cast(cblas.blasint) a.stride,
					b.ptr,
					cast(cblas.blasint) b.stride,
					& beta,
					d.ptr,
					cast(cblas.blasint) d.stride);
				else
				cblas.gemm(
					cblas.Order.RowMajor,
					cblas.Transpose.NoTrans,
					cblas.Transpose.NoTrans,
					cast(cblas.blasint) m,
					cast(cblas.blasint) n,
					cast(cblas.blasint) k,
					alpha,
					a.ptr,
					cast(cblas.blasint) a.stride,
					b.ptr,
					cast(cblas.blasint) b.stride,
					beta,
					d.ptr,
					cast(cblas.blasint) d.stride);

			}
			sw.stop;

			auto newns = sw.peek.to!Duration.total!"nsecs".to!double;
			//writefln("_BLAS (amount of threads is unknown): %5s GFLOPS", (m * n * k * 2) / newns);

			nsecsBLAS = min(newns, nsecsBLAS);

		}

		/// GLAS
		auto nsecsGLAS = double.max;
		foreach(_; 0..count)
		{
			StopWatch sw;
			sw.start;
			static if (__VERSION__ < 2072)
				gemm(alpha, cast(Slice!(2, const(C)*))a, cast(Slice!(2, const(C)*))b, beta, c);
			else
				gemm(alpha, a, b, beta, c);
			sw.stop;
			auto newns = sw.peek.to!Duration.total!"nsecs".to!double;
			//writefln("_GLAS (single thread)               : %5s GFLOPS", (m * n * k * 2) / newns);
			nsecsGLAS = min(newns, nsecsGLAS);
		}
		
		/// Result
		writefln("%s,%s,%s", m, (m * n * k * 2) / nsecsGLAS, (m * n * k * 2) / nsecsBLAS);
	}
}

enum bool isComplex(C)
     = is(Unqual!C == creal)
    || is(Unqual!C == cdouble)
    || is(Unqual!C == cfloat);

void fillRNG(T)(Slice!(2, T*) sl)
{
	import std.random;
	foreach(row; sl)
	foreach(ref e; row)
	{
		static if(isComplex!T)
		{
			e = uniform(-100, 100) + uniform(-100, 100) * 1i;
		}
		else
		{
			e = cast(T) uniform(-100, 100);
		}
	}
}

static if (__VERSION__ < 2072)
{
    import std.experimental.ndslice.slice: Slice;
    import ldc.attributes : fastmath;

    @fastmath:
    pragma(inline, true)
    @property T* ptr(T)(Slice!(2, T*) slice)
    {
        return &(slice.front.front());
    }

    static if (__VERSION__ < 2072)
    pragma(inline, true)
    @property T* ptr(T)(Slice!(1, T*) slice)
    {
        return &(slice.front());
    }
}