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

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

/* common data */

extern unsigned int chrstatus[256];
extern unsigned int chrmap_2bit[256];
extern unsigned int chrmap_4bit[256];
extern char chrmap_complement[256];
extern unsigned char chrmap_5bit_aa[256];

/* functions in query.c */

void query_open(const char * filename);

int query_getnext(char ** head, long * head_len,
                  char ** seq, long * seq_len, long * qno,
		  long * qsize);

void query_close();

long query_getfilesize();

long query_getfilepos();

/* functions in util.c */

long gcd(long a, long b);
void fatal(const char * format, ...);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
char * xstrchrnul(char *s, int c);
unsigned long hash_fnv_1a_64(unsigned char * s, unsigned long n);
long getusec(void);
void show_rusage();

/* functions in score.c */

long * score_read_file_aa (const char * filename, unsigned char * chrmap);
long score (long * score_matrix, int d, int q);
long score_chrs (long * score_matrix, char d, char q, unsigned char * chrmap);
void score_matrix_put(long * score_matrix, unsigned char * chrmap);
