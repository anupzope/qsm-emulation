# qsm-emulation

This project is designed to verify Queue Streaming Model (QSM) by emulating it on traditional cache-based processor.  It demonstrates execution time predictability of a QSM algorithm using the SpMV multiplication as an example.

## Prerequisites

1. x86_64 processor with support for AVX instruction set
   * Note: Tested with Intel(R) Xeon(R) E5-2680 v2 processor, not sure of the other processors.
2. Intel(R) 2018 C++ compiler
3. CMake >= 3.9

## Building

Follow these steps:

1. Extract the source code into a directory, let us say, `$SRC`.
2. Create a build directory, let us say `$BUILD`
3. `cd $BUILD`
4. `cmake $SRC -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=Release`
5. `make`

## Running

The software has two programs:

1. `extract-matrix-structure`
2. `spmv`

`extract-matrix-structure` is used to extract structure of a sparse matrix from [Matrix Market](https://math.nist.gov/MatrixMarket/) file. Usage is as follows.

```
extract-matrix-structure CASE-NAME
```

Here, `CASE_NAME` is the matrix file name without the `.mtx` extension. If the input matrix is symmetric, it is converted into general format. Then it is converted into compressed-row-storage (CRS) format. The resulting structure is written in a file named `CASE-NAME.crs`. Format of `.crs` file is as follows.

```
nrows ncols nnz
ncolidxs1 colidx1 ...
ncolidxs2 colidx1 ...
...
```

After this, use the program `spmv` to read the `.crs` file and performs the SpMV multiplication `[X1-0, ..., X1-s] = A * [X0-0, ..., X0-s]`. Here, a number of input vectors are multiplied by the same sparse matrix `A` to produce output vectors. This type of SpMV multiplication is common in block conjugate gradient method or s-step iterative methods.

The types of matrix and vector elements can be adjusted by modifying the typedefs `VT` (for vector element type) and `MT` (for matrix element type) in file `spmv.cpp`.

```
using VT = <vector type>;
using MT = <matrix-type>;
```

Also, the number of vectors need manual adjustment. If `s == 1`, undefine the macro `SPMV_PACK_TEST` defined in `spmv.cpp`. If `s > 1`, define the macro and adjust the number of vectors by repeating `VT` `s` times in the template arguments to the typedef `PT` as follows.

```
using PT = pack_traits<VT, VT, VT, VT, VT>; // for s = 5
```

Before running `spmv`, make sure that the number of OpenMP threads and their affinity is set such that the NUMA effects do not interfere with the performance measurements.

Usage of `spmv`:

```
spmv [-?] [--batch=STRING] [--lookahead=NUM] [--lookback=NUM]
[--nbitsperstage=NUM] [--nrepeat=NUM] [--ntailbits=NUM] [--out]
[--help] [--usage] [--version=NUM] CASE-NAME
```

The parameters have following meaning.

 Argument | Meaning 
 -------- | ------- 
 `batch` | The file to which the performance metrics are appended. 
 `lookahead` | Look-ahead distance in terms of the values of type `VT`. 
 `lookback` | Look-back distance in terms of the values of type `VT`. 
 `nbitsperstage` | Number of bits processed by radix partitioning in each pass. 
 `ntailbits` | Number of bits processed by the radix partitioning from the L3 cache. 
 `out` | if present, the result of the multiplication is dumped into an output file. 
 `nrepeat` | Number of times the test is repeated (must be >= 2). 

To get meaningful results, set `lookback` to a value that is large enough to fill half of the L2 cache and the `lookahead` value half of the `lookback` value. That is,

```
lookback = ceil(0.5*L2size/(sizeof(VT)*s))
```

and

```
lookahead = lookback/2
```

Set the `ntailbits` to a value that is large enough to fill half of the L3 cache in the last pass of the radix partitioning used for copy stream generation. That is,

```
ntailbits = floor(log2(ceil(0.5*L3size/(sizeof(VT)*s))))
```

`nbitsperstage` is set to 5 for the Intel(R) Xeon(R) E5-2680 v2 processor.
