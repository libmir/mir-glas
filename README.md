[![Dub downloads](https://img.shields.io/dub/dt/mir-glas.svg)](http://code.dlang.org/packages/mir-glas)
[![License](https://img.shields.io/dub/l/mir-glas.svg)](http://code.dlang.org/packages/mir-glas)
[![Bountysource](https://www.bountysource.com/badge/team?team_id=145399&style=bounties_received)](https://www.bountysource.com/teams/libmir)
[![Gitter](https://img.shields.io/gitter/room/libmir/public.svg)](https://gitter.im/libmir/public)

[![Latest version](https://img.shields.io/dub/v/mir-glas.svg)](http://code.dlang.org/packages/mir-glas)

[![Circle CI](https://circleci.com/gh/libmir/mir-glas.svg?style=svg)](https://circleci.com/gh/libmir/mir-glas)
[![Build Status](https://travis-ci.org/libmir/mir-glas.svg?branch=master)](https://travis-ci.org/libmir/mir-glas)

[![Benchmarks](http://blog.mir.dlang.io/images/bench_csingle.svg)](http://blog.mir.dlang.io/glas/benchmark/openblas/2016/09/23/glas-gemm-benchmark.html)

# glas
LLVM-accelerated Generic Linear Algebra Subprograms (GLAS)

## Description
GLAS is a C library written in Dlang. No C++/D runtime is required, but libc that which is available everywhere.

It provides

 1. [BLAS](http://netlib.org/blas/) (Basic Linear Algebra Subprograms) API.
 2. GLAS (Generic Linear Algebra Subprograms) API.

CBLAS API can be provided by linking with [Natlib's CBLAS](http://netlib.org/blas/#_cblas) library.


## Usage

`mir-glas` can be used like common C libraries. It should be linked with `mir-cpuid`.
A compiler, for example GCC, may require `mir-cpuid` be passed after `mir-glas`: `-lmir-glas -lmir-cpuid`.


### Headers

C/C++ headers are located in [`include/`](https://github.com/libmir/mir-glas/tree/master/include/).
D headers are located in [`source/`](https://github.com/libmir/mir-glas/tree/master/source/).

There are two files:

 1. `glas/fortran.h` / `glas/fortran.d` - for Netilb's BLAS API
 2. `glas/ndslice.h` / `glas/ndslice.d` - for GLAS API

### GLAS API and Documentation

Documentation can be found at http://docs.glas.dlang.io/.

## Installation

#### Compiler installation

[LDC (LLVM D Compiler)](https://github.com/ldc-developers/ldc) >= `1.1.0` is required to build a project.
`1.1.0` version is not released yet.
You may want to build LDC from source or use [LDC 1.1.0 beta 2](https://github.com/ldc-developers/ldc/releases/tag/v1.1.0-beta2).
Beta 2 generates a lot of warnings that can be ignored. Beta 3 is not supported.

LDC binaries contains two compilers: ldc2 and ldmd2. The last one should be used.

Recent LDC packages come with [dub package manager](http://code.dlang.org/docs/commandline).
dub is used to build the project.

#### Mir CPUID
[Mir CPUID](https://github.com/libmir/mir-cpuid) is CPU Identification Routines.

Download `mir-cpuid`
```shell
dub fetch mir-cpuid --cache=local
```

Change the directory
```shell
cd mir-cpuid-<CURRENT-CPUID-VERSION>/mir-cpuid
```

Build `mir-cpuid`
```
dub build --build=release-nobounds --compiler=ldmd2 --build-mode=singleFile --parallel --force
```
You may need to add `--arch=x86_64`, if you use windows.

Copy `libmir-cpuid.a` to your project or add its directory to the library path.

#### Mir GLAS

Download `mir-glas`
```shell
dub fetch mir-glas --cache=local
```

Change the directory
```shell
cd mir-glas-<CURRENT-GLAS-VERSION>/mir-glas
```

Build `mir-glas`
```
dub build --config=static --build=target-native --compiler=ldmd2 --build-mode=singleFile --parallel --force
```
You may need to add `--arch=x86_64` if you use windows.

Copy `libmir-glas.a` to your project or add its directory to the library path.

## Status

 - [x] CI testing with Netlib's CBLAS test suite.
 - [ ] CI testing with Netlib's LAPACKE test suite.
 - [ ] Multi-threading
 - [ ] GPU back-end
 - [ ] Shared library support - requires only DUB configuration fixes.
 - [ ] Level 3 - matrix-matrix operations
   - [x] GEMM - matrix matrix multiply
   - [x] SYMM, HEMM - symmetric / hermitian matrix matrix multiply
   - [ ] SYRK, HERK, SYR2K, HER2K - symmetric / hermitian rank-k / rank-2k update to a matrix
   - [ ] TRMM - triangular matrix matrix multiply
   - [ ] TRSM - solving triangular matrix with multiple right hand sides
 - [ ] Level 2 - matrix-vector operations
   - [ ] GEMV - matrix vector multiply
   - [ ] GBMV - banded matrix vector multiply
   - [ ] HEMV - hermitian matrix vector multiply
   - [ ] HBMV - hermitian banded matrix vector multiply
   - [ ] HPMV - hermitian packed matrix vector multiply
   - [ ] TRMV - triangular matrix vector multiply
   - [ ] TBMV - triangular banded matrix vector multiply
   - [ ] TPMV - triangular packed matrix vector multiply
   - [ ] TRSV - solving triangular matrix problems
   - [ ] TBSV - solving triangular banded matrix problems
   - [ ] TPSV - solving triangular packed matrix problems
   - [ ] GERU - performs the rank 1 operation `A := alpha*x*y' + A`
   - [ ] GERC - performs the rank 1 operation `A := alpha*x*conjg( y' ) + A`
   - [ ] HER - hermitian rank 1 operation `A := alpha*x*conjg(x') + A`
   - [ ] HPR - hermitian packed rank 1 operation `A := alpha*x*conjg( x' ) + A`
   - [ ] HER2 - hermitian rank 2 operation
   - [ ] HPR2 - hermitian packed rank 2 operation
 - [ ] Level 1 - vector-vector and scalar operations. Note: [Mir](https://github.com/libmir/mir) already provides generic implementation.
   - [ ] ROTG - setup Givens rotation
   - [ ] ROTMG - setup modified Givens rotation
   - [ ] ROT - apply Givens rotation
   - [ ] ROTM - apply modified Givens rotation
   - [ ] SWAP - swap x and y
   - [x] SCAL - `x = a*x`. Note: requires addition optimization for complex numbers.
   - [ ] COPY - copy x into y
   - [ ] AXPY - `y = a*x + y`
   - [ ] DOT - dot product
   - [ ] NRM2 - Euclidean norm
   - [ ] ASUM - sum of absolute values
   - [ ] IAMAX - index of max abs value

## Questions & Answers

#### Why GLAS calls "Generic ..."?

 1. GLAS has generic internal implementation, which can be easily portable to any other architecture with minimal efforts (5 minutes).
 2. GLAS API provides more functionality comparing with BLAS.
 3. It is written is in Dlang using generic programming.

#### Why it is better then other BLAS Open Source Libraries like OpenBLAS and Eigen?

 1. It is [faster](http://blog.mir.dlang.io/glas/benchmark/openblas/2016/09/23/glas-gemm-benchmark.html).
 2. GLAS API is more user-friendly and does not require additional data coping.
 3. It does not require C++ runtime comparing with Eigen.
 4. It does not require platform specific optimizations like Eigen intrinsics micro kernels and OpenBLAS assembler macro kernels.
 5. It has simple implementation, which can easily ported and extended.

#### Why GLAS does not have Lazy Evaluation and Aliasing like Eigen?

GLAS is more low-level library comparing with Eigen. For example, GLAS can be Eigen BLAS back-end in the future.
Lazy Evaluation and Aliasing can be easily implemented in D.
Explicit composition operations can be found at [mir.ndslice.algorithm](http://docs.mir.dlang.io/latest/mir_ndslice_algorithm.html#mapSlice).