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

static salt_fasta_t ** of = NULL;
static int of_count = 0;
/*
extern unsigned int chrstatus[256];

static FILE * fp;
static char line[LINEALLOC];

static long no = -1;

static char * head = 0;
static char * seq = 0;

static long head_len = 0;
static long seq_len = 0;

static long head_alloc = 0;
static long seq_alloc = 0;

static long filesize = 0;

static long lineno;

static long stripped_count;
static long stripped[256];

regex_t q_regexp;
*/
static salt_fasta_t * init_file_descriptor(int i)
{
  of = (salt_fasta_t **)xrealloc(of, (i+1)*sizeof(salt_fasta_t *));

  of[i] = xmalloc(sizeof(salt_fasta_t),8);

  of[i]->no   = -1;
  of[i]->head =  0;
  of[i]->seq  =  0;

  of[i]->head_len =  0;
  of[i]->seq_len  =  0;

  of[i]->head_alloc =  0;
  of[i]->seq_alloc  =  0;

  of[i]->filesize =  0;

  of[i]->ofid = i;

  return of[i];
}

long salt_fasta_getfilesize(salt_fasta_t * fd)
{
  return fd->filesize;
}

long salt_fasta_getfilepos(salt_fasta_t * fd)
{
  return ftell(fd->fp);
}

salt_fasta_t * salt_fasta_open(const char * filename)
{
  salt_fasta_t * fd = init_file_descriptor(of_count++);

  if (regcomp(&(fd->q_regexp),
              "(^|;)size=([0-9]+)(;|$)",
              REG_EXTENDED))
    fatal("Regular expression compilation failed");

  //unsigned long line_len;
  /* allocate space */

  fd->head = NULL;
  fd->seq = NULL;

  fd->head_len = 0;
  fd->seq_len = 0;

  fd->head_alloc = MEMCHUNK;
  fd->seq_alloc = MEMCHUNK;

  fd->head = (char *) xmalloc((size_t)(fd->head_alloc), 
                                    SALT_ALIGNMENT_SSE);
  fd->seq  = (char *) xmalloc((size_t)(fd->seq_alloc), 
                                    SALT_ALIGNMENT_SSE);

  fd->no = -1;

  /* open queyfile */
  fd->fp = NULL;
  fd->fp = fopen(filename, "r");
  if (!fd->fp)
    fatal("Error: Unable to open query file (%s)", filename);

  if (fseek(fd->fp, 0, SEEK_END))
    fatal("Error: Unable to seek in query file (%s)", filename);

  fd->filesize = ftell(fd->fp);

  rewind(fd->fp);

  fd->line[0] = 0;
  fgets(fd->line, LINEALLOC, fd->fp);
  fd->lineno = 1;

  fd->stripped_count = 0;
  for(int i=0; i<256; i++)
    fd->stripped[i] = 0;

  return fd;
}

void salt_fasta_close(salt_fasta_t * fd)
{
  /* Warn about stripped chars */

  if (fd->stripped_count)
    {
      fprintf(stderr, "Warning: invalid characters stripped from query:");
      for (int i=0; i<256;i++)
        if (fd->stripped[i])
          fprintf(stderr, " %c(%ld)", i, fd->stripped[i]);
      fprintf(stderr, "\n");
    }

  fclose(fd->fp);

  if (fd->seq)
    free(fd->seq);
  if (fd->head)
    free(fd->head);

  fd->head = 0;
  fd->seq = 0;

  if (fd->ofid != of_count-1)
  {
    of[fd->ofid] = of[of_count-1];
  }

  free(fd);
  of_count--;
}

int salt_fasta_getnext(salt_fasta_t * fd, char ** head, long * head_len,
                       char ** seq, long * seq_len, long * qno,
                       long * qsize)
{
  char msg[200];
  while (fd->line[0])
    {
      /* read header */

      if (fd->line[0] != '>')
        fatal("Illegal header line in query fasta file");

      long headerlen = xstrchrnul(fd->line+1, '\n') - (fd->line+1);
      fd->head_len = headerlen;

      if (headerlen + 1 > fd->head_alloc)
        {
          fd->head_alloc = headerlen + 1;
          fd->head = (char *) xrealloc(fd->head,
                                             (size_t)(fd->head_alloc));
        }

      memcpy(fd->head, fd->line + 1, (size_t)headerlen);
      fd->head[headerlen] = 0;

      /* read size/abundance annotation */

      regmatch_t pmatch[4];

      if (!regexec(&(fd->q_regexp), fd->head, 4, pmatch, 0))
        {
          unsigned long size = atol(fd->head + pmatch[2].rm_so);
          if (size > 0)
            * qsize = size;
          else
            fatal("Size annotation zero in query sequence");
        }
      else
        *qsize = 1;

      /* get next line */

      fd->line[0] = 0;
      fgets(fd->line, LINEALLOC, fd->fp);
      fd->lineno++;

      /* read sequence */

      fd->seq_len = 0;

      while (fd->line[0] && (fd->line[0] != '>'))
        {
          char c;
          char m;
          char * p = fd->line;

          while((c = *p++))
            {
              m = chrstatus[(int)c];
              switch(m)
                {
                case 0:
                  /* character to be stripped */
                  fd->stripped_count++;
                  fd->stripped[(int)c]++;
                  break;

                case 1:
                  /* legal character */
                  if (fd->seq_len + 1 > fd->seq_alloc)
                    {
                      fd->seq_alloc += MEMCHUNK;
                      fd->seq = (char *) xrealloc(fd->seq, (size_t)(fd->seq_alloc));
                    }
                  *(fd->seq + fd->seq_len) = c;
                  fd->seq_len++;

                  break;

                case 2:
                  /* fatal character */
                  if (c>=32)
                    snprintf(msg, 200, "illegal character '%c' on line %ld in the query file", c, fd->lineno);
                  else
                    snprintf(msg, 200, "illegal unprintable character %#.2x (hexadecimal) on line %ld in the query file", c, fd->lineno);
                  fatal(msg);
                  break;

                case 3:
                  /* silently stripped chars */
                  break;

                }
            }

          fd->line[0] = 0;
          fgets(fd->line, LINEALLOC, fd->fp);
          fd->lineno++;
        }

      /* add zero after sequence */

      if (fd->seq_len + 1 > fd->seq_alloc)
        {
          fd->seq_alloc += MEMCHUNK;
          fd->seq = (char *) xrealloc(fd->seq, (size_t)fd->seq_alloc);
        }
      *(fd->seq + fd->seq_len) = 0;




      fd->no++;
      *head = fd->head;
      *seq = fd->seq;
      *head_len = fd->head_len;
      *seq_len = fd->seq_len;
      *qno = fd->no;

      return 1;
    }

  return 0;
}


