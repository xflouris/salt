#include "salt.h"

/* please note that these functions will return a pointer to a buffer
   allocated here for the query header and sequence. This buffers will
   be overwritten on the next call of query_getnext. */

#define MEMCHUNK 4096
#define LINEALLOC LINE_MAX

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

long query_getfilesize()
{
  return query_filesize;
}

long query_getfilepos()
{
  return ftell(query_fp);
}

void query_open(const char * filename)
{
  if (regcomp(&q_regexp, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
    fatal("Regular expression compilation failed");
  
  //unsigned long query_line_len;
  /* allocate space */

  query_head = NULL;
  query_seq = NULL;

  query_head_len = 0;
  query_seq_len = 0;

  query_head_alloc = MEMCHUNK;
  query_seq_alloc = MEMCHUNK;

  query_head = (char *) xmalloc((size_t)query_head_alloc);
  query_seq = (char *) xmalloc((size_t)query_seq_alloc);

  query_no = -1;

  /* open query file */
  query_fp = NULL;
  query_fp = fopen(filename, "r");
  if (!query_fp)
    fatal("Error: Unable to open query file (%s)", filename);

  if (fseek(query_fp, 0, SEEK_END))
    fatal("Error: Unable to seek in query file (%s)", filename);

  query_filesize = ftell(query_fp);
  
  rewind(query_fp);

  query_line[0] = 0;
  fgets(query_line, LINEALLOC, query_fp);
  query_lineno = 1;

  query_stripped_count = 0;
  for(int i=0; i<256; i++)
    query_stripped[i] = 0;
}

void query_close()
{
  /* Warn about stripped chars */

  if (query_stripped_count)
    {
      fprintf(stderr, "Warning: invalid characters stripped from query:");
      for (int i=0; i<256;i++)
        if (query_stripped[i])
          fprintf(stderr, " %c(%ld)", i, query_stripped[i]);
      fprintf(stderr, "\n");
    }
  
  fclose(query_fp);
  
  if (query_seq)
    free(query_seq);
  if (query_head)
    free(query_head);

  query_head = 0;
  query_seq = 0;
}

int query_getnext(char ** head, long * head_len,
                  char ** seq, long * seq_len, long * qno,
                  long * qsize)
{
  char msg[200];
  while (query_line[0])
    {
      /* read header */

      if (query_line[0] != '>')
        fatal("Illegal header line in query fasta file");
      
      long headerlen = xstrchrnul(query_line+1, '\n') - (query_line+1);
      query_head_len = headerlen;

      if (headerlen + 1 > query_head_alloc)
        {
          query_head_alloc = headerlen + 1;
          query_head = (char *) xrealloc(query_head, (size_t)query_head_alloc);
        }

      memcpy(query_head, query_line + 1, (size_t)headerlen);
      query_head[headerlen] = 0;

      /* read size/abundance annotation */

      regmatch_t pmatch[4];

      if (!regexec(&q_regexp, query_head, 4, pmatch, 0))
        {
          unsigned long size = atol(query_head + pmatch[2].rm_so);
          if (size > 0)
            * qsize = size;
          else
            fatal("Size annotation zero in query sequence");
        }
      else
        *qsize = 1;

      /* get next line */

      query_line[0] = 0;
      fgets(query_line, LINEALLOC, query_fp);
      query_lineno++;

      /* read sequence */

      query_seq_len = 0;

      while (query_line[0] && (query_line[0] != '>'))
        {
          char c;
          char m;
          char * p = query_line;

          while((c = *p++))
            {
              m = chrstatus[(int)c];
              switch(m)
                {
                case 0:
                  /* character to be stripped */
                  query_stripped_count++;
                  query_stripped[(int)c]++;
                  break;

                case 1:
                  /* legal character */
                  if (query_seq_len + 1 > query_seq_alloc)
                    {
                      query_seq_alloc += MEMCHUNK;
                      query_seq = (char *) xrealloc(query_seq, (size_t)query_seq_alloc);
                    }
                  *(query_seq + query_seq_len) = c;
                  query_seq_len++;

                  break;

                case 2:
                  /* fatal character */
                  if (c>=32)
                    snprintf(msg, 200, "illegal character '%c' on line %ld in the query file", c, query_lineno);
                  else
                    snprintf(msg, 200, "illegal unprintable character %#.2x (hexadecimal) on line %ld in the query file", c, query_lineno);
                  fatal(msg);
                  break;

                case 3:
                  /* silently stripped chars */
                  break;

                }
            }

          query_line[0] = 0;
          fgets(query_line, LINEALLOC, query_fp);
          query_lineno++;
        }

      /* add zero after sequence */

      if (query_seq_len + 1 > query_seq_alloc)
        {
          query_seq_alloc += MEMCHUNK;
          query_seq = (char *) xrealloc(query_seq, (size_t)query_seq_alloc);
        }
      *(query_seq + query_seq_len) = 0;




      query_no++;
      *head = query_head;
      *seq = query_seq;
      *head_len = query_head_len;
      *seq_len = query_seq_len;
      *qno = query_no;

      return 1;
    }
  
  return 0;
}


