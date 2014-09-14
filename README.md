# SALT

## Introduction

The aim of this project is to develop a highly optimized library and toolkit for common algorithms and operations in biological sequence analysis. The new tool should have:

* open source code with an appropriate open source license
* 64-bit design supporting the most recent and popular SIMD architectures
* possibly a GPU port

## List of available algorithms

As a start the library should contain the following important algorithms:

* Smith-Waterman
* Needleman-Wunsch-Sellers
* Czech-Flouri algorithm for optimal prefix-suffix matching
* Optimal prefix-suffix matching (mismatches only).

The first three algorithms should be implemented with the affine gap penalties option. We should also implement both variants - scores and penalties.

## Implementation Details

Although this should be a library in the end, we should develop separate testbeds for every algorithm.

**Smith-Waterman**: A generic version of SW will be coded at first, mainly for correctness tests. Vectorized version will include a single-to-single sequence comparison (query/database), a single-to-multi sequence comparison (single database multiple queries), and finally a multiple single-to-single sequence comparison, i.e. multiple pairs in parallel. SSE and AVX version will be implemented. Moreover, we should implement two variations of the these methods, where the sizes of the matrix cells are bytes (short sequences), and words (lon sequences).

**Needleman-Wunsch-Sellers**: Same as SW.

**Czech-Flouri**: Same as SW.

**Optimal prefix-suffix matching (mismatches)**: Same as SW, but since the only dependency in the cells computation is the upper-left diagonal, we don't need the the three variations of the problem. For this one we should code only a highly-optimized single-to-single sequence comparison.
