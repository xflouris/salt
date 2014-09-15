#include "salt.h"

/*

  Optimal prefix-suffix matching with mismatches only

  finds the best overlap with a minimum cost
  there should be positive costs/penalties for mismatches
  matches should have zero cost (0)

  dseq: the database/horizontal sequence
  qseq: the query/vertical sequence

  typical costs:
  match: 0
  mismatch: 1

  input

  dseq: pointer to start of database sequence
  dend: pointer after database sequence
  qseq: pointer to start of query sequence
  qend: pointer after query sequence
  score_matrix: 32x32 matrix of longs with scores for aligning two symbols

  output

  psmscore: the best possible score of the alignment
  overlaplen: length of the best overlap
  matchcase: 0 if the best score was achieved by aligning a prefix of query with
             a suffix of database, otherwise 1 if a prefix of database was
             aligned with a suffix of query.

*/

unsigned long qarray_alloc = 0;
unsigned long darray_alloc = 0;

unsigned long * qarray;
unsigned long * darray;

void prefsuf_mis_init()
{
  qarray = 0;
  darray = 0;

  qarray_alloc = 0;
  darray_alloc = 0;
}

void prefsuf_mis_exit()
{
  if (qarray)
    free(qarray);
  if (darray)
    free(darray);
}

void prefsuf_mis (char * dseq,
                  char * dend,
                  char * qseq,
                  char * qend,
                  long * score_matrix,
                  long * psmscore,
                  long * overlaplen,
                  long * matchcase)
{
  long h, n, score, len;
  unsigned long *qa, *da;

  long dlen = dend - dseq;
  long qlen = qend - qseq;

  qarray_alloc = qlen * sizeof(long);
  darray_alloc = dlen * sizeof(long);

  qarray = (unsigned long *) xrealloc(qarray, qarray_alloc);
  darray = (unsigned long *) xrealloc(darray, darray_alloc);
  da = darray;

  memset (qarray, 0, qarray_alloc);

  long i, j;

  /* compute the matrix */
  for (j = 0; j < dlen; ++j) 
  {
    h = 0;
    qa = qarray;
    
    for (i = 0; i < qlen; ++i)
    {
      n = *qa;
      h += score_matrix[(dseq[j] << 5) + qseq[i]];

      *qa++ = h;
      h = n;
    }
    
    *da++ = qarray[qlen - 1];
  }

  /* pick the best values */
  *matchcase = 0;
  for (i = 0, score = qarray[0]; i < qlen; ++i)
  {
    if (qarray[i] >= score)
    {
      len = i;
      score = qarray[i];
    }
  }

  for (i = 0; i < dlen; ++i)
  {
    if (darray[i] >= score)
    {
      len = i;
      score = darray[i];
      *matchcase = 1;
    }
  }
  
  *psmscore = score;
  *overlaplen = len;
}

