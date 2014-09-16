#include "salt.h"

static char * progname;
char * opt_list_reads;
char * infilename;

long opt_help;
long opt_version;

static char progheader[80];
static char * cmdline;

void args_init(int argc, char **argv)
{
  /* Set defaults */
  progname = argv[0];

  opt_help       = 0;
  opt_version    = 0;
  opt_list_reads = 0;

  static struct option long_options[] =
  {
    {"help",       no_argument,       0, 0 },
    {"version",    no_argument,       0, 0 },
    {"list-reads", required_argument, 0, 0 },
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

       default:
         fatal("Internal error in option parsing");
     }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands = 0;
  if (opt_list_reads)
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


void getentirecommandline(int argc, char ** argv)
{
  int len = 0;
  for (int i = 0; i < argc; i++)
    len += strlen(argv[i]);

  cmdline = (char*) xmalloc(len + argc + 1);

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

//   query_open("file.fa");

  show_header();

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

  score_chrmap_set(&chrmap_5bit_aa);
  score_matrix_read_aa("../data/score_matrix");
  score_matrix_put();
  fprintf(stdout, "%lu\n", score_chr('A', 'C'));


  if (opt_help)
   {
     cmd_help();
   }
  else if (opt_list_reads)
   {
     while (query_getnext(&head, &head_len,
                           &seq, &seq_len, &qno, &qsize))

      {

        fprintf(stdout, "%s\n%s\n\n", head, seq);
      }
     query_close();
   }

  return (EXIT_SUCCESS);
}
