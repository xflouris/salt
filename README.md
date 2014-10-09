# SALT (Sequence Analysis Library and Toolkit)

## Introduction

The aim of this project is to develop a highly optimized library and toolkit for common algorithms and operations in biological sequence analysis. The new tool should have:

* open source code with an appropriate open source license
* 64-bit design supporting the most recent and popular SIMD architectures
* possibly a GPU port

## List of available algorithms

As a start the library should contain the following important algorithms:

* Smith-Waterman
* Needleman-Wunsch-Sellers
* Optimal prefix-suffix (overlap) detection given a scoring matrix.
* Optimal prefix-suffix (overlap) detection given a scoring matrix and quality scores.

The first two algorithms should be implemented with the affine gap penalties option. We should also implement both variants - scores and penalties.

## Implementation Details

Although this should be a library in the end, we should develop separate testbeds for every algorithm.

**Smith-Waterman**: A generic version of SW will be coded at first, mainly for correctness tests. Vectorized version will include a single-to-single sequence comparison (query/database), a single-to-multi sequence comparison (single database multiple queries), and finally a multiple single-to-single sequence comparison, i.e. multiple pairs in parallel. SSE and AVX version will be implemented. Moreover, we should implement two variations of the these methods, where the sizes of the matrix cells are bytes (short sequences), and words (lon sequences).

**Needleman-Wunsch-Sellers**: Same as SW.

**Optimal prefix-suffix matching (mismatches)**: Same as SW, but since the only dependency in the cells computation is the upper-left diagonal, we don't need the the three variations of the problem. For this one we should code only a highly-optimized single-to-single sequence comparison.

## SALT Toolkit command line options

General options:

* `--help`
* `--version`

Merging reads:

* `--overlap <filename>`

## SALT license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

SALT binaries may include code from the [zlib](http://www.zlib.net) library copyright Jean-loup Gailly and Mark Adler.

SALT binaries may include code from the [bzip2](http://www.bzip.org) library copyright Julian R. Seward.

## Code

The code is written in C++ but most of it is actually C with some C++ syntax conventions.

    File     | Description
-------------|------
**maps.c** | Various character mapping arrays
**overlap_plain.c** | Detection of optimal overlap (prefix-suffix) between two sequences (Non-vectorized).
**overlap_plain_vec.c** | SIMD implementation of optimal overlap detection between two sequences.
**popcount.c** | SIMD implementation of the popcount instruction.
**query.cc** | Reads the fasta file containing the query sequences.
**salt.c** | Toolkit file, for testing the functions of SALT.
**util.c** | Various common utility functions.

## Bugs

SALT has not been tested comprehensively yet. All bug reports are highly appreciated.


## The SALT team

The following people have contributed to SALT:

* Tom&aacute;&scaron; Flouri
* Lucas Czech
* Kassian Kobert
* Jiajie Zhang
