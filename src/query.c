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

/* please note that these functions will return a pointer to a buffer
   allocated here for the query header and sequence. This buffers will
   be overwritten on the next call of query_getnext. */

#define MEMCHUNK 4096
#define LINEALLOC LINE_MAX

typedef struct
{
  FILE * query_fp;
  char query_line[LINEALLOC];

  long query_no;

  char * query_head;
  char * query_seq;

  long query_head_len;
  long query_seq_len;

  long query_head_alloc;
  long query_seq_alloc;

  long query_filesize;

  long query_lineno;

  long query_stripped_count;
  long query_stripped[256];

  regex_t q_regexp;
} salt_FASTA;

static salt_FASTA * open_files = NULL;
static int open_files_count = 0;
/*
extern unsigned int chrstatus[256];

static FILE * query_fp;
static char query_line[LINEALLOC];

static long query_no = -1;

static char * query_head = 0;
static char * query_seq = 0;

static long query_head_len = 0;
static long query_seq_len = 0;

static long query_head_alloc = 0;
static long query_seq_alloc = 0;

static long query_filesize = 0;

static long query_lineno;

static long query_stripped_count;
static long query_stripped[256];

regex_t q_regexp;
*/

static void init_open_files(int index)
{
  open_files = xrealloc((void *)open_files, index * sizeof(salt_FASTA));

  open_files[index-1].query_no   = -1;
  open_files[index-1].query_head =  0;
  open_files[index-1].query_seq  =  0;

  open_files[index-1].query_head_len =  0;
  open_files[index-1].query_seq_len  =  0;

  open_files[index-1].query_head_alloc =  0;
  open_files[index-1].query_seq_alloc  =  0;

  open_files[index-1].query_filesize =  0;
}

long salt_fasta_getfilesize(int id)
{
  return open_files[id].query_filesize;
}

long salt_fasta_getfilepos(int id)
{
  return ftell(open_files[id].query_fp);
}

int salt_fasta_open(const char * filename)
{
  init_open_files(++open_files_count);

  if (regcomp(&(open_files[open_files_count - 1].q_regexp),
              "(^|;)size=([0-9]+)(;|$)",
              REG_EXTENDED))
    fatal("Regular expression compilation failed");

  //unsigned long query_line_len;
  /* allocate space */

  open_files[open_files_count - 1].query_head = NULL;
  open_files[open_files_count - 1].query_seq = NULL;

  open_files[open_files_count - 1].query_head_len = 0;
  open_files[open_files_count - 1].query_seq_len = 0;

  open_files[open_files_count - 1].query_head_alloc = MEMCHUNK;
  open_files[open_files_count - 1].query_seq_alloc = MEMCHUNK;

  open_files[open_files_count - 1].query_head = (char *) xmalloc((size_t)(open_files[open_files_count - 1].query_head_alloc), SALT_ALIGNMENT_SSE);
  open_files[open_files_count - 1].query_seq = (char *) xmalloc((size_t)(open_files[open_files_count - 1].query_seq_alloc), SALT_ALIGNMENT_SSE);

  open_files[open_files_count - 1].query_no = -1;

  /* open query file */
  open_files[open_files_count - 1].query_fp = NULL;
  open_files[open_files_count - 1].query_fp = fopen(filename, "r");
  if (!open_files[open_files_count - 1].query_fp)
    fatal("Error: Unable to open query file (%s)", filename);

  if (fseek(open_files[open_files_count - 1].query_fp, 0, SEEK_END))
    fatal("Error: Unable to seek in query file (%s)", filename);

  open_files[open_files_count - 1].query_filesize = ftell(open_files[open_files_count - 1].query_fp);

  rewind(open_files[open_files_count - 1].query_fp);

  open_files[open_files_count - 1].query_line[0] = 0;
  fgets(open_files[open_files_count - 1].query_line, LINEALLOC, open_files[open_files_count - 1].query_fp);
  open_files[open_files_count - 1].query_lineno = 1;

  open_files[open_files_count - 1].query_stripped_count = 0;
  for(int i=0; i<256; i++)
    open_files[open_files_count - 1].query_stripped[i] = 0;

  return open_files_count - 1;
}

void salt_fasta_close(int id)
{
  /* Warn about stripped chars */

  if (open_files[id].query_stripped_count)
    {
      fprintf(stderr, "Warning: invalid characters stripped from query:");
      for (int i=0; i<256;i++)
        if (open_files[id].query_stripped[i])
          fprintf(stderr, " %c(%ld)", i, open_files[id].query_stripped[i]);
      fprintf(stderr, "\n");
    }

  fclose(open_files[id].query_fp);

  if (open_files[id].query_seq)
    free(open_files[id].query_seq);
  if (open_files[id].query_head)
    free(open_files[id].query_head);

  open_files[id].query_head = 0;
  open_files[id].query_seq = 0;
}

int salt_fasta_getnext(int id, char ** head, long * head_len,
                       char ** seq, long * seq_len, long * qno,
                       long * qsize)
{
  char msg[200];
  while (open_files[id].query_line[0])
    {
      /* read header */

      if (open_files[id].query_line[0] != '>')
        fatal("Illegal header line in query fasta file");

      long headerlen = xstrchrnul(open_files[id].query_line+1, '\n') - (open_files[id].query_line+1);
      open_files[id].query_head_len = headerlen;

      if (headerlen + 1 > open_files[id].query_head_alloc)
        {
          open_files[id].query_head_alloc = headerlen + 1;
          open_files[id].query_head = (char *) xrealloc(open_files[id].query_head,
                                                        (size_t)(open_files[id].query_head_alloc));
        }

      memcpy(open_files[id].query_head, open_files[id].query_line + 1, (size_t)headerlen);
      open_files[id].query_head[headerlen] = 0;

      /* read size/abundance annotation */

      regmatch_t pmatch[4];

      if (!regexec(&(open_files[id].q_regexp), open_files[id].query_head, 4, pmatch, 0))
        {
          unsigned long size = atol(open_files[id].query_head + pmatch[2].rm_so);
          if (size > 0)
            * qsize = size;
          else
            fatal("Size annotation zero in query sequence");
        }
      else
        *qsize = 1;

      /* get next line */

      open_files[id].query_line[0] = 0;
      fgets(open_files[id].query_line, LINEALLOC, open_files[id].query_fp);
      open_files[id].query_lineno++;

      /* read sequence */

      open_files[id].query_seq_len = 0;

      while (open_files[id].query_line[0] && (open_files[id].query_line[0] != '>'))
        {
          char c;
          char m;
          char * p = open_files[id].query_line;

          while((c = *p++))
            {
              m = chrstatus[(int)c];
              switch(m)
                {
                case 0:
                  /* character to be stripped */
                  open_files[id].query_stripped_count++;
                  open_files[id].query_stripped[(int)c]++;
                  break;

                case 1:
                  /* legal character */
                  if (open_files[id].query_seq_len + 1 > open_files[id].query_seq_alloc)
                    {
                      open_files[id].query_seq_alloc += MEMCHUNK;
                      open_files[id].query_seq = (char *) xrealloc(open_files[id].query_seq, (size_t)(open_files[id].query_seq_alloc));
                    }
                  *(open_files[id].query_seq + open_files[id].query_seq_len) = c;
                  open_files[id].query_seq_len++;

                  break;

                case 2:
                  /* fatal character */
                  if (c>=32)
                    snprintf(msg, 200, "illegal character '%c' on line %ld in the query file", c, open_files[id].query_lineno);
                  else
                    snprintf(msg, 200, "illegal unprintable character %#.2x (hexadecimal) on line %ld in the query file", c, open_files[id].query_lineno);
                  fatal(msg);
                  break;

                case 3:
                  /* silently stripped chars */
                  break;

                }
            }

          open_files[id].query_line[0] = 0;
          fgets(open_files[id].query_line, LINEALLOC, open_files[id].query_fp);
          open_files[id].query_lineno++;
        }

      /* add zero after sequence */

      if (open_files[id].query_seq_len + 1 > open_files[id].query_seq_alloc)
        {
          open_files[id].query_seq_alloc += MEMCHUNK;
          open_files[id].query_seq = (char *) xrealloc(open_files[id].query_seq, (size_t)open_files[id].query_seq_alloc);
        }
      *(open_files[id].query_seq + open_files[id].query_seq_len) = 0;




      open_files[id].query_no++;
      *head = open_files[id].query_head;
      *seq = open_files[id].query_seq;
      *head_len = open_files[id].query_head_len;
      *seq_len = open_files[id].query_seq_len;
      *qno = open_files[id].query_no;

      return 1;
    }

  return 0;
}


