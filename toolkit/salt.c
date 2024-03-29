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

#include "salt.h"
#include <assert.h>
#include <time.h>

static char * progname;
char * opt_list_reads;
char * opt_overlap_file;
char * opt_algorithm;

int    opt_run_test;
int    opt_runs;
int    opt_reads_min_len;
int    opt_reads_max_len;
int    opt_min_overlap;
int    opt_verbose;
int    opt_seed;

char * infilename;

long opt_help;
long opt_version;

static char progheader[80];
static char * cmdline;

#define SCORE_MATRIX_SIZE 32

void args_init(int argc, char **argv)
{
  /* Set defaults */
  progname = argv[0];

  opt_help          = 0;
  opt_version       = 0;
  opt_list_reads    = 0;
  opt_overlap_file  = 0;

  opt_algorithm     = xstrdup_aligned("CPU",8);
  opt_run_test      = 0;
  opt_runs          = 10;
  opt_reads_min_len = 150;
  opt_reads_max_len = 300;
  opt_min_overlap   = 20;
  opt_verbose       = 0;
  opt_seed          = time(NULL);

  static struct option long_options[] =
  {
    {"help",          no_argument,       0, 0 },
    {"version",       no_argument,       0, 0 },
    {"list-reads",    required_argument, 0, 0 },
    {"overlap",       required_argument, 0, 0 },
    {"algorithm",     required_argument, 0, 0 },
    {"test",          no_argument,       0, 0 },
    {"runs",          required_argument, 0, 0 },
    {"reads_min_len", required_argument, 0, 0 },
    {"reads_max_len", required_argument, 0, 0 },
    {"min_overlap",   required_argument, 0, 0 },
    {"verbose",       no_argument,       0, 0 },
    {"seed",          required_argument, 0, 0 },
    { 0, 0, 0, 0 }
  };

  int option_index = 0;
  int c;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    switch (option_index)
     {
       case 0:
         /* help */
         opt_help = 1;
         break;

       case 1:
         /* version */
         opt_version = 1;
         break;

       case 2:
         /* list-reads */
         opt_list_reads = optarg;
         break;

       case 3:
         /* overlap */
         opt_overlap_file = optarg;
         break;

       case 4:
         /* select algorithm */
         free (opt_algorithm);
         opt_algorithm = optarg;
         break;

       case 5:
         /* test */
         opt_run_test = 1;
         break;

       case 6:
         /* number of runs */
         opt_runs = atoi(optarg);
         break;

       case 7:
         /* reads_min_len */
         opt_reads_min_len = atoi(optarg);
         break;

       case 8:
         /* reads_max_len */
         opt_reads_max_len = atoi(optarg);
         break;

       case 9:
         /* min_overlap */
         opt_min_overlap = atoi(optarg);
         break;

       case 10:
         /* verbose */
         opt_verbose = 1;
         break;

       case 11:
         /* seeed - dancehall caballeros*/
         opt_seed = atoi(optarg);
         break;

       default:
         fatal("Internal error in option parsing");
     }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands = 0;
  if (opt_list_reads)
    commands++;
  if (opt_overlap_file)
    commands++;
  if (opt_run_test)
    commands++;
  if (opt_help)
    commands++;
  if (opt_version)
    commands++;

  if (commands == 0)
    opt_version = 1;

  if (commands > 1)
    fatal("More than one command specified");
}

void cmd_help()
{
  fprintf (stderr,
           "Usage: %s [OPTIONS]\n", progname);
  fprintf (stderr,
           "\n"
           "General options:\n"
           "  --help                      display help information\n"
           "  --version                   display version information\n"
           "  --list-reads FILENAME       display reads in input fasta file\n"
          );
}

void convert(char * s)
{
  char c;
  char m;
  char * p = s;
  int i = 0;

  while ((c = *p++))
  {
    if ((m = chrmap_2bit[(int)c]) >= 0)
      *(s + i++) = m;
    else
      fatal("Illegal character in sequence.");
  }
}

// simpler version...
/*
void convert (char* s)
{
    char  c;
    while ((c = *s)) {
        *s = chrmap_2bit[(int)c];
        s++;
    }
}
*/

void init_scoring_matrices (
    long *scorematrix_long,
    WORD *scorematrix_word,
    char *scorematrix_char
)
{
    for (int i = 0; i < SCORE_MATRIX_SIZE; ++i) {
        for (int j = 0; j < SCORE_MATRIX_SIZE; ++j) {
            if (i == j) {
                scorematrix_long[(i<<5) + j] = 1;
                scorematrix_word[(i<<5) + j] = 1;
                scorematrix_char[(i<<5) + j] = 1;
            } else {
                scorematrix_long[(i<<5) + j] = -1;
                scorematrix_word[(i<<5) + j] = -1;
                scorematrix_char[(i<<5) + j] = -1;
            }
        }
    }
}

void cmd_overlap()
{
  char * head;
  long head_len;
  char * seq[2];
  long seq_len[2];
  long qno;
  long qsize;
  long scorematrix_long[SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE] __attribute__((aligned(SALT_ALIGNMENT_MAX)));
  WORD scorematrix_word[SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE] __attribute__((aligned(SALT_ALIGNMENT_MAX)));
  char scorematrix_char[SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE] __attribute__((aligned(SALT_ALIGNMENT_MAX)));
  salt_fasta_t * fd;

  long psmscore = 0, overlaplen = 0, matchcase = 0;

  fd = salt_fasta_open(opt_overlap_file);

  /* get first sequence */
  salt_fasta_getnext(fd, &head, &head_len,
                     &seq[0], &seq_len[0], &qno, &qsize);

  seq[0] = xstrdup_aligned(seq[0], SALT_ALIGNMENT_MAX);

  /* get second sequence */
  salt_fasta_getnext(fd, &head, &head_len,
                     &seq[1], &seq_len[1], &qno, &qsize);

  seq[1] = xstrdup_aligned(seq[1], SALT_ALIGNMENT_MAX);

  /* setup scoring matrix */
  init_scoring_matrices (scorematrix_long, scorematrix_word, scorematrix_char);

  printf ("dbs: %s len: %ld\n", seq[0], seq_len[0]);
  printf ("qry: %s len: %ld\n", seq[1], seq_len[1]);

  /* convert to a number representation */
  convert(seq[0]);
  convert(seq[1]);

  salt_overlap_nuc4(seq[0], seq[0] + seq_len[0],
                    seq[1], seq[1] + seq_len[1],
                    (long *)scorematrix_long,
                    &psmscore,
                    &overlaplen,
                    &matchcase);


  printf("\nCPU       : psmscore: %ld, overlaplen: %ld, matchcase: %ld\n", psmscore, overlaplen, matchcase);

  salt_overlap_nuc4_sse_8((BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                          (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                          scorematrix_char,
                          &psmscore,
                          &overlaplen,
                          &matchcase);

  printf("SSE   8bit: psmscore: %ld, overlaplen: %ld, matchcase: %ld\n", psmscore, overlaplen, matchcase);

  salt_overlap_nuc4_sse_16((BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                           (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                           scorematrix_word,
                           &psmscore,
                           &overlaplen,
                           &matchcase);

  printf("SSE  16bit: psmscore: %ld, overlaplen: %ld, matchcase: %ld\n", psmscore, overlaplen, matchcase);

  salt_overlap_nuc4_avx2_16((BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                            (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                            scorematrix_word,
                            &psmscore,
                            &overlaplen,
                            &matchcase);

  printf("AVX2 16bit: psmscore: %ld, overlaplen: %ld, matchcase: %ld\n", psmscore, overlaplen, matchcase);

  salt_fasta_close(fd);
}

/*
void cmd_run_test ()
{
    // get command line options
    int runs          = opt_runs;
    int reads_min_len = opt_reads_min_len;
    int reads_max_len = opt_reads_max_len;
    int min_overlap   = opt_min_overlap;
    int verbose       = opt_verbose;

    assert (min_overlap < reads_min_len);

    // user output showing options
    printf ("Starting with\n");
    printf ("   algorithm:     %s\n", opt_algorithm);
    printf ("   runs:          %i\n", runs);
    printf ("   reads_min_len: %i\n", reads_min_len);
    printf ("   reads_max_len: %i\n", reads_max_len);
    printf ("   min_overlap:   %i\n", min_overlap);

    // prepare variables
    char* seq[2];
    long  seq_len[2];
    int overlap;
    long psmscore = 0, overlaplen = 0, matchcase = 0;

    reads_max_len++; // functions expect min<max, because the interval is [min,max)
    seq[0] = (char*) xmalloc (reads_max_len * sizeof(char), SALT_ALIGNMENT_MAX);
    seq[1] = (char*) xmalloc (reads_max_len * sizeof(char), SALT_ALIGNMENT_MAX);

    long scorematrix_long[SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE] __attribute__((aligned(SALT_ALIGNMENT_MAX)));
    WORD scorematrix_word[SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE] __attribute__((aligned(SALT_ALIGNMENT_MAX)));
    char scorematrix_char[SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE] __attribute__((aligned(SALT_ALIGNMENT_MAX)));

    init_scoring_matrices (scorematrix_long, scorematrix_word, scorematrix_char);

    // generate two random overlapping sequences
    seq_len[0] = random_int_range (reads_min_len, reads_max_len);
    seq_len[1] = random_int_range (reads_min_len, reads_max_len);
    overlap    = random_int_range (min_overlap,   seq_len[0] + seq_len[1] - min_overlap);
    generate_pair (seq[0], seq_len[0], seq[1], seq_len[1], overlap);
    seq[0][seq_len[0]] = '\0';
    seq[1][seq_len[1]] = '\0';

    // user output the sequences
    if (verbose) {
        printf ("=============================================================\n");
        printf ("dbs: %s len: %ld\n", seq[0], seq_len[0]);
        printf ("qry: %s len: %ld\n", seq[1], seq_len[1]);
        printf ("overlap: %u\n\n", overlap);
    }

    // convert sequences to a number representation
    convert(seq[0]);
    convert(seq[1]);

    clock_t stime = clock();

    // test-analyze many different random sequences
    for (int i = 0; i < runs; i++) {
        // select algorithm and run it
        if (strcmp(opt_algorithm, "CPU") == 0) {
            salt_overlap_nuc4(
                seq[0], seq[0] + seq_len[0],
                seq[1], seq[1] + seq_len[1],
                (long *)scorematrix_long,
                &psmscore,
                &overlaplen,
                &matchcase
            );
        } else if (strcmp(opt_algorithm, "SSE8") == 0) {
            salt_overlap_nuc4_sse_8(
                (BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                scorematrix_char,
                &psmscore,
                &overlaplen,
                &matchcase
            );
        } else if (strcmp(opt_algorithm, "SSE16") == 0) {
            salt_overlap_nuc4_sse_16(
                (BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                scorematrix_word,
                &psmscore,
                &overlaplen,
                &matchcase
            );
        } else if (strcmp(opt_algorithm, "AVX8") == 0) {
            salt_overlap_nuc4_avx2_8(
                (BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                scorematrix_char,
                &psmscore,
                &overlaplen,
                &matchcase
            );
        } else if (strcmp(opt_algorithm, "AVX16") == 0) {
            salt_overlap_nuc4_avx2_16(
                (BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                scorematrix_word,
                &psmscore,
                &overlaplen,
                &matchcase
            );
        } else {
            printf ("Unknown algorithm: '%s'. Aborting.\n", opt_algorithm);
            exit(0);
        }
    }

    // user output the results
    int noverlap;
    if (matchcase) {
        noverlap = seq_len[0] + seq_len[1] - overlaplen;
    } else {
        noverlap = overlaplen;
    }
    if (verbose) {
        printf ("%s: psmscore: %3ld, overlaplen: %3ld, matchcase: %ld, actual overlap: %3d\n", opt_algorithm, psmscore, overlaplen, matchcase, noverlap);
    }

    if (overlap!=noverlap) {
        printf ("Warning: wrong overlaps %i vs %i! Seed: %i", overlap, noverlap, opt_seed);
    }

    // user output timing information
    clock_t etime = clock();
    printf ("\nFinished.\nClock time: %li\n", etime-stime);

    // clean up
    xfree (seq[0]);
    xfree (seq[1]);
}
*/

void getentirecommandline(int argc, char ** argv)
{
  int len = 0;
  for (int i = 0; i < argc; i++)
    len += strlen(argv[i]);

  cmdline = (char*) xmalloc(len + argc + 1, SALT_ALIGNMENT_SSE);
  cmdline[0] = 0;

  for (int i = 0; i < argc; i++)
   {
     strcat(cmdline, argv[i]);
     strcat(cmdline, " ");
   }
}

void fillheader()
{
  snprintf(progheader, 80,
           "%s %s_%s",
           PROG_NAME, PROG_VERSION, PROG_ARCH);
}

void show_header()
{
  fprintf(stdout, "                 ____ \n"
                  "     _________ _/ / /_\n"
                  "    / ___/ __ `/ / __/\n"
                  "   (__  ) /_/ / / /_ \n"
                  "  /____/\\__,_/_/\\__/ \n");
  fprintf (stdout, "%s\n\n", progheader);
}

int main (int argc, char * argv[])
{
  char * head;
  long head_len;
  char * seq;
  long seq_len;
  long qno;
  long qsize;

  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);
  srand(opt_seed);

  show_header();

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_list_reads)
  {
    salt_fasta_t * fd = salt_fasta_open(opt_list_reads);
    while (salt_fasta_getnext(fd, &head, &head_len,
                              &seq, &seq_len, &qno, &qsize))
    {
      fprintf(stdout, "%s\n%s\n\n", head, seq);
    }
    salt_fasta_close(fd);
  }
  else if (opt_overlap_file)
  {
    cmd_overlap();
  } /*else if (opt_run_test) {
      cmd_run_test();
  }*/

  return (EXIT_SUCCESS);
}
