/*
    Copyright (C) 2014 Tomas Flouri & Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <x86intrin.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <ctype.h>
#include <regex.h>
#include <fcntl.h>
#include <unistd.h>

#ifndef __APPLE__
#include <sys/sysinfo.h>
#endif

/* platform specific */

#ifdef _WIN32
#define SALT_EXPORT __declspec(dllexport)
#else
#define SALT_EXPORT
#endif


/* constants */

#define PROG_NAME "salt"
#define PROG_VERSION "v0.0.0"

#define SALT_ALIGNMENT_SSE 16
#define SALT_ALIGNMENT_AVX 32
#define SALT_ALIGNMENT_MAX 32 // used whenever it is yet unclear which alignment is needed

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#define MEMCHUNK 4096
#define LINEALLOC LINE_MAX

/* structures and data types */

typedef unsigned int UINT32;
typedef short WORD;
typedef unsigned char BYTE;

typedef struct
{
  FILE * fp;
  char line[LINEALLOC];

  long no;

  char * head;
  char * seq;

  long head_len;
  long seq_len;

  long head_alloc;
  long seq_alloc;

  long filesize;

  long lineno;

  long stripped_count;
  long stripped[256];

  long ofid;

  regex_t q_regexp;
} salt_fasta_t;


/* common data */

extern unsigned int chrstatus[256];
extern unsigned int chrmap_2bit[256];
extern unsigned int chrmap_4bit[256];
extern char chrmap_complement[256];
extern unsigned char chrmap_5bit_aa[256];

#ifdef __cplusplus
extern "C" {
#endif

/* functions in query.c */

SALT_EXPORT salt_fasta_t * salt_fasta_open(const char * filename);

SALT_EXPORT int salt_fasta_getnext(salt_fasta_t * fd, char ** head, long * head_len,
                                   char ** seq, long * seq_len, long * qno,
                                   long * qsize);

SALT_EXPORT void salt_fasta_close(salt_fasta_t * fd);

SALT_EXPORT long salt_fasta_getfilesize(salt_fasta_t * fd);

SALT_EXPORT long salt_fasta_getfilepos(salt_fasta_t * fd);

/* functions in util.c */

SALT_EXPORT long gcd(long a, long b);

SALT_EXPORT long roundup(long offset, long step);

SALT_EXPORT void fatal(const char * format, ...);

SALT_EXPORT void * xmalloc(size_t size, size_t alignment);

SALT_EXPORT void * xrealloc(void * ptr, size_t size);

SALT_EXPORT void xfree (void* ptr);

SALT_EXPORT char * xstrchrnul(char *s, int c);

SALT_EXPORT long getusec(void);

SALT_EXPORT void show_rusage();

SALT_EXPORT void * xstrdup_aligned(char * s, size_t alignment);

/* functions in overlap_nuc.c */

SALT_EXPORT void salt_overlap_nuc4(char * dseq, char * dend,
                                   char * qseq, char * qend,
                                   long * score_matrix,
                                   long * psmscore,
                                   long * overlaplen,
                                   long * matchcase);

/* functions in overlap_nuc4_avx2_8.c */

SALT_EXPORT void salt_overlap_avx2_8(BYTE * dseq,
                                     BYTE * dend,
                                     BYTE * qseq,
                                     BYTE * qend,
                                     char * score_matrix,
                                     long * psmscore,
                                     long * overlaplen,
                                     long * matchcase);

/* functions in overlap_nuc4_sse_16bit.c */

SALT_EXPORT void salt_overlap_nuc4_avx2_16(BYTE * dseq,
                                           BYTE * dend,
                                           BYTE * qseq,
                                           BYTE * qend,
                                           WORD * score_matrix,
                                           long * psmscore,
                                           long * overlaplen,
                                           long * matchcase);

/* functions in overlap_nuc4_sse_8.c */

SALT_EXPORT void salt_overlap_nuc4_sse_8(BYTE * dseq,
                                         BYTE * dend,
                                         BYTE * qseq,
                                         BYTE * qend,
                                         char * score_matrix,
                                         long * psmscore,
                                         long * overlaplen,
                                         long * matchcase);

/* functions in overlap_nuc4_sse_16.c */

SALT_EXPORT void salt_overlap_nuc4_sse_16(BYTE * dseq, BYTE * dend,
                                          BYTE * qseq, BYTE * qend,
                                          WORD * score_matrix,
                                          long * pmscore,
                                          long * overlaplen,
                                          long * matchcase);

/* functions in popcount.c */

SALT_EXPORT void pprint(__m128i x);

#ifdef __cplusplus
}
#endif
