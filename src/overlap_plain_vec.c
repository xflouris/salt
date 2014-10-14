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
#define shft 3

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

static WORD * qprofile   = NULL;
static WORD * hh         = NULL;
static WORD * ee         = NULL;

static long qprofile_len = 0;
static long ee_len       = 0;
static long hh_len       = 0;

static void qprofile_fill16_avx(WORD * score_matrix_word,
                                  BYTE * qseq,
                                  BYTE * qend)
{
  WORD * offset;
  long qlen = qend - qseq;
  long padded_len = roundup(qlen, 16);
  long i;

  if (padded_len > qprofile_len)
  {
    free(qprofile);
    qprofile = xmalloc(4*padded_len*sizeof(WORD), SALT_ALIGNMENT_SSE);
    qprofile_len = padded_len;
  }

  /* currently only for DNA with A,C,G,T as 0,1,2,3 */
  for (i = 0, offset = qprofile; i < 4; offset += padded_len, i++)
  {
    for (long j = 0; j < qlen; ++j)
    {
      offset[j] = score_matrix_word[(i << 5) + qseq[j]];
    }
    for (long j = qlen; j < padded_len; ++j)
    {
      offset[j] = 0;
    }
  }
}

//#ifdef DEBUG
void pprint_avx(__m256i x)
{
  short * p = (short *) & x;

  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d ", *p++);
  printf("%04d", *p++);
}

void pshow_avx(char * name, __m256i x)
{
  printf("%s: ", name);
  pprint_avx(x);
  printf("\n");
}
//#endif

void salt_overlap_plain16_avx2(BYTE * dseq,
                               BYTE * dend,
                               BYTE * qseq,
                               BYTE * qend,
                               WORD * score_matrix,
                               long * psmscore,
                               long * overlaplen,
                               long * matchcase)
{
  long len = 0;
  WORD score = 0;
  long dlen = dend - dseq;
  long qlen = qend - qseq;
  long qlen_padded = roundup(qlen,16);

  char c;

  if (qlen_padded > hh_len)
  {
    free(hh);
    hh = xmalloc(qlen_padded*sizeof(WORD), SALT_ALIGNMENT_AVX);
    hh_len = qlen_padded;
  }
  if (dlen > ee_len)
  {
    free(ee);
    ee = xmalloc(roundup(dlen,8)*sizeof(WORD), SALT_ALIGNMENT_AVX);
    ee_len = dlen;
  }

  qprofile_fill16_avx(score_matrix, qseq, qend);

  __m256i X, H, T1, xmm0, xmm1, xmm2, xmm3, xmm4;

  xmm2 = _mm256_set_epi16(0xffff, 0xffff, 0xffff, 0xffff,
                          0xffff, 0xffff, 0xffff, 0xffff,
                          0x0000, 0x0000, 0x0000, 0x0000,
                          0x0000, 0x0000, 0x0000, 0x0000);

  xmm3 = _mm256_set_epi16(0x0000, 0x0000, 0x0000, 0x0000,
                          0x0000, 0x0000, 0x0000, 0x0000,
                          0xffff, 0xffff, 0xffff, 0xffff,
                          0xffff, 0xffff, 0xffff, 0xffff);

  xmm0 = _mm256_setzero_si256();

  for (long i = 0; i < qlen_padded; i += 16)
  {
    _mm256_store_si256((__m256i *)(hh + i), xmm0);
  }

  WORD * lastbyte= hh+qlen-1;

  for (long j = 0; j < dlen; j++)
  {
    X = xmm0;
    c = dseq[j];
    for (long i = 0; i < qlen_padded; i += 16)
     {
       H  = _mm256_load_si256((__m256i *)(hh+i));

       xmm1 = _mm256_permute2x128_si256(H,H, _MM_SHUFFLE(0,0,0,3));
       xmm4 = _mm256_and_si256(xmm1, xmm3);
       T1 = _mm256_alignr_epi8(xmm4,xmm0,0x1e);

       xmm4 = _mm256_and_si256(xmm1, xmm2);
       H    = _mm256_alignr_epi8(H,xmm4,14);

       H  = _mm256_or_si256(H,X);
       X  = T1;

       xmm1 = _mm256_load_si256((__m256i *)(qprofile + c*qlen_padded + i));
       H = _mm256_add_epi16(H,xmm1);

       _mm256_store_si256((__m256i *)(hh+i),H);
     }

     *(ee+j) = *lastbyte;
  }

  // pick best value
  *matchcase = 0;
  score = hh[0];

  // normal case
  for (long i = 0; i < qlen; ++i)
  {
    if (hh[i] >= score)
    {
      len = i+1;
      score = hh[i];
    }
  }

  // run through case
  for (long i = 0; i < dlen; ++i)
  {
    if (ee[i] >= score)
    {
      len = i+1;
      score = ee[i];
      *matchcase = 1;
    }
  }

  *psmscore = score;
  *overlaplen = len;
}
