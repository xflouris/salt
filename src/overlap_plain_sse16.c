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

static WORD * qprofile = NULL;
static WORD * hh = NULL;
static WORD * ee = NULL;

static long qprofile_len = 0;
static long ee_len = 0;
static long hh_len = 0;

static void pprint_sse(__m128i x)
{
  short * p = (short *) & x;

  printf("%02d ", *p++);
  printf("%02d ", *p++);
  printf("%02d ", *p++);
  printf("%02d ", *p++);
  printf("%02d ", *p++);
  printf("%02d ", *p++);
  printf("%02d ", *p++);
  printf("%02d", *p++);
}

static void pshow_sse(char * name, __m128i x)
{
  printf("%s: ", name);
  pprint_sse(x);
  printf("\n");
}

#if NONVEC
static void qprofile_fill16_sse(WORD * score_matrix_word,
                                BYTE * qseq,
                                BYTE * qend)
{
  WORD * offset;
  long qlen = qend - qseq;
  long padded_len = roundup(qlen, 8);
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
#else
static void qprofile_fill16_sse(WORD * score_matrix_word,
                                BYTE * qseq,
                                BYTE * qend)
{
  long qlen = qend - qseq;
  long padded_len = roundup(qlen, 16);

  __m128i xmm0, xmm1, xmm2,  xmm3,  xmm4,  xmm5,  xmm6,   xmm7;
  __m128i xmm8, xmm9, xmm10, xmm11, xmm12, xmm13, xmm14;

  if (padded_len > qprofile_len)
  {
    free(qprofile);
    qprofile = xmalloc(4*padded_len*sizeof(WORD), SALT_ALIGNMENT_SSE);
    qprofile_len = padded_len;
  }

  xmm0 = _mm_set_epi8(0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,
                      0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02);

  xmm1 = _mm_load_si128((__m128i *)(score_matrix_word+0));  /* A */
  xmm2 = _mm_load_si128((__m128i *)(score_matrix_word+32)); /* C */
  xmm3 = _mm_load_si128((__m128i *)(score_matrix_word+64)); /* G */
  xmm4 = _mm_load_si128((__m128i *)(score_matrix_word+96)); /* T */

  /* pick first two elements from each of the four nucleotides */
  xmm1 = _mm_unpacklo_epi32(xmm1,xmm2);
  xmm2 = _mm_unpacklo_epi32(xmm3,xmm4);

  xmm3 = _mm_unpacklo_epi64(xmm1,xmm2);  /* contains AC for A,C,G,T */
  xmm4 = _mm_unpackhi_epi64(xmm1,xmm2);  /* contains GT for A,C,G,T */
  
  xmm5 = _mm_set_epi8(0x80, 0x07, 0x80, 0x06, 0x80, 0x05, 0x80, 0x04,
                      0x80, 0x03, 0x80, 0x02, 0x80, 0x01, 0x80, 0x00);

  xmm6 = _mm_set_epi8(0x80, 0x0f, 0x80, 0x0e, 0x80, 0x0d, 0x80, 0x0c,
                      0x80, 0x0b, 0x80, 0x0a, 0x80, 0x09, 0x80, 0x08);

  xmm7 = _mm_set_epi8(0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00,
                      0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00);

  for (long i = 0; i < padded_len; i += 16)
  {
    /* load sequence */
    xmm8 = _mm_load_si128((__m128i *)(qseq+i));

    /* left shift values by 2 */
    xmm8 = _mm_add_epi8(xmm8,xmm8);
    xmm8 = _mm_add_epi8(xmm8,xmm8);

    /* get first 8 values and expand them to 16-bit */
    xmm9  = _mm_shuffle_epi8(xmm8,xmm5);
    xmm10 = _mm_slli_si128(xmm9, 1);
    xmm9 = _mm_or_si128(xmm9,xmm10);
    xmm9 = _mm_add_epi8(xmm9,xmm7);


    /* get next 8 values and expand them to 16-bit */
    xmm10 = _mm_shuffle_epi8(xmm8,xmm6);
    xmm11 = _mm_slli_si128(xmm10, 1);
    xmm10 = _mm_or_si128(xmm10,xmm11);
    xmm10 = _mm_add_epi8(xmm10,xmm7);
    
    /* get A component */
    xmm11 = _mm_shuffle_epi8(xmm3,xmm9); /* store it */

    /* get G component */
    xmm12 = _mm_shuffle_epi8(xmm4,xmm9); /* store it */
    
    xmm9 = _mm_add_epi8(xmm9,xmm0);

    /* get C component */
    xmm13 = _mm_shuffle_epi8(xmm3,xmm9);

    /* get T component */
    xmm14 = _mm_shuffle_epi8(xmm4,xmm9);

    _mm_store_si128((__m128i *)(qprofile+0*padded_len+i),xmm11);
    _mm_store_si128((__m128i *)(qprofile+1*padded_len+i),xmm13);
    _mm_store_si128((__m128i *)(qprofile+2*padded_len+i),xmm12);
    _mm_store_si128((__m128i *)(qprofile+3*padded_len+i),xmm14);

    /* get A component */
    xmm11 = _mm_shuffle_epi8(xmm3,xmm10); /* store it */

    /* get G component */
    xmm12 = _mm_shuffle_epi8(xmm4,xmm10); /* store it */
    
    xmm10 = _mm_add_epi8(xmm10,xmm0);

    /* get C component */
    xmm13 = _mm_shuffle_epi8(xmm3,xmm10);

    /* get T component */
    xmm14 = _mm_shuffle_epi8(xmm4,xmm10);

    _mm_store_si128((__m128i *)(qprofile+0*padded_len+i+8),xmm11);
    _mm_store_si128((__m128i *)(qprofile+1*padded_len+i+8),xmm13);
    _mm_store_si128((__m128i *)(qprofile+2*padded_len+i+8),xmm12);
    _mm_store_si128((__m128i *)(qprofile+3*padded_len+i+8),xmm14);
  }
}
#endif

void salt_overlap_plain16_sse (BYTE * dseq,
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
  long qlen_padded = roundup(qlen,8);
  long qlen_padded16 = roundup(qlen,16); /* needed for query profile */

  char c;

  if (qlen_padded > hh_len)
  {
    free(hh);
    hh = xmalloc(qlen_padded*sizeof(WORD), SALT_ALIGNMENT_SSE);
    hh_len = qlen_padded;
  }
  if (dlen > ee_len)
  {
    free(ee);
    ee = xmalloc(roundup(dlen,8)*sizeof(WORD), SALT_ALIGNMENT_SSE);
    ee_len = dlen;
  }

  qprofile_fill16_sse(score_matrix,
                      qseq,
                      qend);

  __m128i X, H, T1, xmm0, xmm1;

  xmm0 = _mm_setzero_si128();

  for (long i = 0; i < qlen_padded; i += 8)
  {
    _mm_store_si128((__m128i *)(hh + i), xmm0);
  }

  WORD * lastbyte= hh+qlen-1;

  for (long j = 0; j < dlen; ++j)
  {
    X = xmm0;
    c = dseq[j];
    for (long i = 0; i < qlen_padded; i += 8)
     {
       H  = _mm_load_si128((__m128i *)(hh+i));
       
       T1 = _mm_srli_si128(H,14);
       H  = _mm_slli_si128(H,2);
       H  = _mm_or_si128(H,X);
       X  = T1;

       /* TODO: Note that the query profile is padded at 16! if
          we use the vectorized version. Otherwise use qlen_padded
          instead of qlen_padded16 */
#if NONVEC
       xmm1 = _mm_load_si128((__m128i *)(qprofile+c*qlen_padded+i));
#else
       xmm1 = _mm_load_si128((__m128i *)(qprofile+c*qlen_padded16+i));
#endif
       H = _mm_add_epi16(H,xmm1);

       _mm_store_si128((__m128i *)(hh+i),H);
     }
    *(ee+j) = *lastbyte;
  }

  /* pick the best values 
     TODO: vectorize it */
  *matchcase = 0;
  score = hh[0];
  for (long i = 0; i < qlen; ++i)
  {
    if (hh[i] >= score)
    {
      len = i+1;
      score = hh[i];
    }
  }

  /* check the run-through case */
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
