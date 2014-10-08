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

static unsigned long qarray_alloc = 0;
static unsigned long darray_alloc = 0;

static unsigned long * qarray;
static unsigned long * darray;

void salt_overlap_plain(char * dseq,
                        char * dend,
                        char * qseq,
                        char * qend,
                        long * score_matrix,
                        long * psmscore,
                        long * overlaplen,
                        long * matchcase)
{
  long i, j, h, n, score, len = 0;
  unsigned long *qa, *da;

  long dlen = dend - dseq;
  long qlen = qend - qseq;

  qarray_alloc = qlen * sizeof(long);
  darray_alloc = dlen * sizeof(long);

  qarray = (unsigned long *) xrealloc(qarray, qarray_alloc);
  darray = (unsigned long *) xrealloc(darray, darray_alloc);
  da = darray;

  memset (qarray, 0, qarray_alloc);

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

  /* pick the best overlap in non run-through case*/
  *matchcase = 0;
  for (i = 0, score = qarray[0]; i < qlen; ++i)
  {
    if (qarray[i] >= score)
    {
      len = i+1;
      score = qarray[i];
    }
  }

  /* check the run-through case */
  for (i = 0; i < dlen; ++i)
  {
    if (darray[i] >= score)
    {
      len = i+1;
      score = darray[i];
      *matchcase = 1;
    }
  }
  
  *psmscore = score;
  *overlaplen = len;
}

