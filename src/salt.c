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

static char * progname;
char * opt_list_reads;
char * opt_overlap_file;
char * infilename;

long opt_help;
long opt_version;

static char progheader[80];
static char * cmdline;

void args_init(int argc, char **argv)
{
  /* Set defaults */
  progname = argv[0];

  opt_help         = 0;
  opt_version      = 0;
  opt_list_reads   = 0;
  opt_overlap_file = 0;


  static struct option long_options[] =
  {
    {"help",       no_argument,       0, 0 },
    {"version",    no_argument,       0, 0 },
    {"list-reads", required_argument, 0, 0 },
    {"overlap",    required_argument, 0, 0 },
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

void cmd_overlap()
{
  char * head;
  long head_len;
  char * seq[2];
  long seq_len[2];
  long qno;
  long qsize;
  long scorematrix[32*32];
  WORD scorematrix_word[32*32];
  int id;

  long psmscore = 0, overlaplen = 0, matchcase = 0;

  id = salt_fasta_open(opt_overlap_file);

  /* get first sequence */
  salt_fasta_getnext(id, &head, &head_len,
                     &seq[0], &seq_len[0], &qno, &qsize);

  seq[0] = xstrdup(seq[0]);

  /* get second sequence */
  salt_fasta_getnext(id, &head, &head_len,
                     &seq[1], &seq_len[1], &qno, &qsize);

  seq[1] = xstrdup(seq[1]);

  /* setup scoring matrix */
  for (int i = 0; i < 32; ++i)
    for (int j = 0; j < 32; ++j)
      if (i == j)
      {
        scorematrix_word[(i<<5) + j] = 1;
        scorematrix[(i<<5) + j] = 1;
      }
      else
      {
        scorematrix_word[(i<<5) + j] = 0;
        scorematrix[(i<<5) + j] = 0;
      }

  printf (" d : %s len: %ld\n", seq[0], seq_len[0]);
  printf ("qry: %s len: %ld\n", seq[1], seq_len[1]);

  /* convert to a number representation */
  convert(seq[0]);
  convert(seq[1]);
  
  salt_overlap_plain(seq[0], seq[0] + seq_len[0],
                     seq[1], seq[1] + seq_len[1],
                     (long *)scorematrix,
                     &psmscore,
                     &overlaplen,
                     &matchcase);


  printf("\n CPU: psmscore: %ld, overlaplen: %ld, matchcase: %ld\n", psmscore, overlaplen, matchcase);

  salt_overlap_plain16_sse((BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                           (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                           scorematrix_word,
                           &psmscore,
                           &overlaplen,
                           &matchcase);

  printf("SSE3: psmscore: %ld, overlaplen: %ld, matchcase: %ld\n", psmscore, overlaplen, matchcase);

  salt_overlap_plain16_avx2((BYTE *)seq[0], (BYTE *)seq[0] + seq_len[0],
                            (BYTE *)seq[1], (BYTE *)seq[1] + seq_len[1],
                            scorematrix_word,
                            &psmscore,
                            &overlaplen,
                            &matchcase);

  printf("AVX2: psmscore: %ld, overlaplen: %ld, matchcase: %ld\n", psmscore, overlaplen, matchcase);

  salt_fasta_close(id);
}

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

  show_header();

<<<<<<< HEAD
=======
  // print stuff for mapping
//   int i;
//   for (i=0; i<256; i++) {
//     if (i%16==0) fprintf(stdout, "\n    ");
//     if(i<33 || i>126) {
//       fprintf(stdout, "0,  ");
//     } else {
//       fprintf(stdout, "%c,  ", i);
//     }
//   }
//   fprintf(stdout, "\n\n");
//

  //score_chrmap_set(chrmap_5bit_aa);
  //score_matrix_read_aa("../data/score_matrix");
  //score_matrix_put();
  //fprintf(stdout, "%lu\n", score_chr('A', 'C'));


>>>>>>> minor changes
  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_list_reads)
  {
    int id = salt_fasta_open(opt_list_reads);
    while (salt_fasta_getnext(id, &head, &head_len,
                              &seq, &seq_len, &qno, &qsize))
    {
      fprintf(stdout, "%s\n%s\n\n", head, seq);
    }
    salt_fasta_close(id);
  }
  else if (opt_overlap_file)
  {
    cmd_overlap();
  }

  return (EXIT_SUCCESS);
}
