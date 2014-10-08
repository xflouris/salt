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

/* constants */

#define PROG_NAME "salt"
#define PROG_VERSION "v0.0.0"

#define SALT_ALIGNMENT_SSE 16
#define SALT_ALIGNMENT_AVX 32

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

/* structures and data types */

typedef unsigned int UINT32;
typedef short WORD;
typedef unsigned char BYTE;

/* common data */

extern unsigned int chrstatus[256];
extern unsigned int chrmap_2bit[256];
extern unsigned int chrmap_4bit[256];
extern char chrmap_complement[256];
extern unsigned char chrmap_5bit_aa[256];

/* functions in query.c */

int salt_fasta_open(const char * filename);

int salt_fasta_getnext(int id, char ** head, long * head_len,
                  char ** seq, long * seq_len, long * qno,
		  long * qsize);

void salt_fasta_close(int id);

long salt_fasta_getfilesize(int id);

long salt_fasta_getfilepos(int id);

/* functions in util.c */

long gcd(long a, long b);
long roundup(long offset, long step);
void fatal(const char * format, ...);
void * xmalloc(size_t size, size_t alignment);
void * xrealloc(void * ptr, size_t size);
void xfree (void* ptr);
char * xstrchrnul(char *s, int c);
unsigned long hash_fnv_1a_64(unsigned char * s, unsigned long n);
long getusec(void);
void show_rusage();
void * xstrdup(char * s);

/* functions in score.c */

void score_chrmap_set(unsigned char * map);
void score_matrix_read_aa (const char * filename);
long score_int (int d, int q);
long score_chr (char d, char q);
void score_matrix_put();

/* functions in overlap_plain.c */

void salt_overlap_plain(char * dseq, char * dend,
                        char * qseq, char * qend,
                        long * score_matrix,
                        long * psmscore,
                        long * overlaplen,
                        long * matchcase);

/* functions in overlap_plain_vec.c */

void salt_overlap_plain16_sse(BYTE * dseq, BYTE * dend,
                              BYTE * qseq, BYTE * qend,
                              WORD * score_matrix,
                              long * pmscore,
                              long * overlaplen,
                              long * matchcase);

void salt_overlap_plain16_avx(BYTE * dseq, BYTE * dend,
                              BYTE * qseq, BYTE * qend,
                              WORD * score_matrix,
                              long * psmscore,
                              long * overlaplen,
                              long * matchcase);

void salt_overlap_plain16_avx2(BYTE * dseq,
                               BYTE * dend,
                               BYTE * qseq,
                               BYTE * qend,
                               WORD * score_matrix,
                               long * psmscore,
                               long * overlaplen,
                               long * matchcase);

/* functions in popcount.c */

void pprint(__m128i x);
