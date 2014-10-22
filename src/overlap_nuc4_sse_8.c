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
static char * qprofile = NULL;
static char * hh = NULL;
static char * ee = NULL;

static long qprofile_len = 0;
static long ee_len = 0;
static long hh_len = 0;

#if 0
/* original non-vectorized version that does not require aligned memory */
static void qprofile_fill8_sse(char * score_matrix_byte,
                                 BYTE * qseq,
                                 BYTE * qend)
{
  char * offset;
  long qlen = qend - qseq;
  long padded_len = roundup(qlen, 16);
  long i;

  if (padded_len > qprofile_len)
  {
    free(qprofile);
    qprofile = xmalloc(4*padded_len*sizeof(char), SALT_ALIGNMENT_SSE);
    qprofile_len = padded_len;
  }

  /* currently only for DNA with A,C,G,T as 0,1,2,3 */
  for (i = 0, offset = qprofile; i < 4; offset += padded_len, i++)
  {
    for (long j = 0; j < qlen; ++j)
    {
      offset[j] = score_matrix_byte[(i << 5) + qseq[j]];
    }
    for (long j = qlen; j < padded_len; ++j)
    {
      offset[j] = 0;
    }
  }
  return qprofile;
}
#endif

/* TODO: Note requires aligned memory for the read */
static void qprofile_fill8_sse(char * score_matrix_byte,
                               BYTE * qseq,
                               BYTE * qend)
{
  long qlen = qend - qseq;
  long padded_len = roundup(qlen, 16);

  __m128i xmm0, xmm1, xmm2,  xmm3,  xmm4,  xmm5,  xmm6,   xmm7;
  __m128i xmm8, xmm9, xmm10, xmm11, xmm12;

  if (padded_len > qprofile_len)
  {
    free(qprofile);
    qprofile = xmalloc(4*padded_len*sizeof(char), SALT_ALIGNMENT_SSE);
    qprofile_len = padded_len;
  }

  /* mask */
  xmm0 = _mm_set_epi8(0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
                      0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01);

  /* scoring matrix */
  xmm1 = _mm_load_si128((__m128i *)(score_matrix_byte+0));  /* A */
  xmm2 = _mm_load_si128((__m128i *)(score_matrix_byte+32)); /* C */
  xmm3 = _mm_load_si128((__m128i *)(score_matrix_byte+64)); /* G */
  xmm4 = _mm_load_si128((__m128i *)(score_matrix_byte+96)); /* T */

  /* pack 16 8-bit values of scoring matrix into one register */
  xmm1 = _mm_unpacklo_epi32(xmm1,xmm2);
  xmm2 = _mm_unpacklo_epi32(xmm3,xmm4);
  xmm3 = _mm_unpacklo_epi64(xmm1,xmm2);

  for (long i = 0; i < padded_len; i += 16)
  {
    /* load sequence */
    xmm5 = _mm_load_si128((__m128i *)(qseq+i));

    /* left shift 16 values by 2 */
    xmm5 = _mm_add_epi8(xmm5,xmm5);
    xmm5 = _mm_add_epi8(xmm5,xmm5);

    xmm6 = _mm_add_epi8(xmm5,xmm0);
    xmm7 = _mm_add_epi8(xmm6,xmm0);
    xmm8 = _mm_add_epi8(xmm7,xmm0);

    /* shuffle */
    xmm9  = _mm_shuffle_epi8(xmm3, xmm5);
    xmm10 = _mm_shuffle_epi8(xmm3, xmm6);
    xmm11 = _mm_shuffle_epi8(xmm3, xmm7);
    xmm12 = _mm_shuffle_epi8(xmm3, xmm8);

    /* store qprofile vectors */
    _mm_store_si128((__m128i *)(qprofile+0*padded_len+i),xmm9);
    _mm_store_si128((__m128i *)(qprofile+1*padded_len+i),xmm10);
    _mm_store_si128((__m128i *)(qprofile+2*padded_len+i),xmm11);
    _mm_store_si128((__m128i *)(qprofile+3*padded_len+i),xmm12);
  }
}

void pprint_sse8(__m128i x)
{
  char * p = (char *) & x;

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

void pshow_sse8(char * name, __m128i x)
{
  printf("%s: ", name);
  pprint_sse8(x);
  printf("\n");
}

#define ONESTEP(H,j)                                                 \
  xmm3 = xmm0;                                                       \
  c = dseq[j];                                                       \
  for (long i = 0; i < qlen_padded; i += 16)                         \
    {                                                                \
      xmm1  = _mm_load_si128((__m128i *)(hh+i));                     \
      xmm2  = _mm_srli_si128(xmm1,15);                               \
      xmm1  = _mm_slli_si128(xmm1,1);                                \
      xmm1  = _mm_or_si128(xmm1,xmm3);                               \
      xmm3  = xmm2;                                                  \
      xmm2  = _mm_load_si128((__m128i *)(qprofile+c*qlen_padded+i)); \
      H     = _mm_add_epi8(xmm1,xmm2);                               \
      _mm_store_si128((__m128i *)(hh+i),H);                          \
    }

#define PROCESS16_NA                                                 \
  for (long j=0; j<dlen16; j+=16)                                    \
  {                                                                  \
    ONESTEP(xmm4, j);                                                \
    ONESTEP(xmm5, j+1);                                              \
    ONESTEP(xmm6, j+2);                                              \
    ONESTEP(xmm7, j+3);                                              \
    ONESTEP(xmm8, j+4);                                              \
    ONESTEP(xmm9, j+5);                                              \
    ONESTEP(xmm10,j+6);                                              \
    ONESTEP(xmm11,j+7);                                              \
                                                                     \
    xmm4  = _mm_unpacklo_epi8(xmm4,xmm5);                            \
    xmm5  = _mm_unpacklo_epi8(xmm6,xmm7);                            \
    xmm6  = _mm_unpacklo_epi8(xmm8,xmm9);                            \
    xmm7  = _mm_unpacklo_epi8(xmm10,xmm11);                          \
                                                                     \
    ONESTEP(xmm8, j+8);                                              \
    ONESTEP(xmm9, j+9);                                              \
    ONESTEP(xmm10,j+10);                                             \
    ONESTEP(xmm11,j+11);                                             \
    ONESTEP(xmm12,j+12);                                             \
    ONESTEP(xmm13,j+13);                                             \
    ONESTEP(xmm14,j+14);                                             \
    ONESTEP(xmm15,j+15);                                             \
                                                                     \
    /* interleave bytes */                                           \
    xmm8  = _mm_unpacklo_epi8(xmm8, xmm9);                           \
    xmm9  = _mm_unpacklo_epi8(xmm10,xmm11);                          \
    xmm10 = _mm_unpacklo_epi8(xmm12,xmm13);                          \
    xmm11 = _mm_unpacklo_epi8(xmm14,xmm15);                          \
                                                                     \
    /* interleave pairs */                                           \
    xmm12 = _mm_unpacklo_epi16(xmm4, xmm5);                          \
    xmm13 = _mm_unpacklo_epi16(xmm6, xmm7);                          \
    xmm14 = _mm_unpacklo_epi16(xmm8, xmm9);                          \
    xmm15 = _mm_unpacklo_epi16(xmm10,xmm11);                         \
                                                                     \
    /* interleave quadruplets */                                     \
    xmm5 = _mm_unpacklo_epi32(xmm12,xmm13);                          \
    xmm6 = _mm_unpacklo_epi32(xmm14,xmm15);                          \
                                                                     \
    /* interleave octets */                                          \
    xmm7 = _mm_unpacklo_epi64(xmm5,xmm6);                            \
                                                                     \
    /* store horizontal */                                           \
    _mm_store_si128((__m128i *)(ee+j),xmm7);                         \
  }

#define PROCESS16(N)                                                 \
  for (long j=0; j<dlen16; j+=16)                                    \
  {                                                                  \
    ONESTEP(xmm4, j);                                                \
    ONESTEP(xmm5, j+1);                                              \
    ONESTEP(xmm6, j+2);                                              \
    ONESTEP(xmm7, j+3);                                              \
    ONESTEP(xmm8, j+4);                                              \
    ONESTEP(xmm9, j+5);                                              \
    ONESTEP(xmm10,j+6);                                              \
    ONESTEP(xmm11,j+7);                                              \
                                                                     \
    /* align at lsb */                                               \
    xmm4  = _mm_alignr_epi8(xmm4, xmm4, N);                          \
    xmm5  = _mm_alignr_epi8(xmm5, xmm5, N);                          \
    xmm6  = _mm_alignr_epi8(xmm6, xmm6, N);                          \
    xmm7  = _mm_alignr_epi8(xmm7, xmm7, N);                          \
    xmm8  = _mm_alignr_epi8(xmm8, xmm8, N);                          \
    xmm9  = _mm_alignr_epi8(xmm9, xmm9, N);                          \
    xmm10 = _mm_alignr_epi8(xmm10,xmm10,N);                          \
    xmm11 = _mm_alignr_epi8(xmm11,xmm11,N);                          \
                                                                     \
    /* interleave bytes */                                           \
    xmm4  = _mm_unpacklo_epi8(xmm4,xmm5);                            \
    xmm5  = _mm_unpacklo_epi8(xmm6,xmm7);                            \
    xmm6  = _mm_unpacklo_epi8(xmm8,xmm9);                            \
    xmm7  = _mm_unpacklo_epi8(xmm10,xmm11);                          \
                                                                     \
    ONESTEP(xmm8, j+8);                                              \
    ONESTEP(xmm9, j+9);                                              \
    ONESTEP(xmm10,j+10);                                             \
    ONESTEP(xmm11,j+11);                                             \
    ONESTEP(xmm12,j+12);                                             \
    ONESTEP(xmm13,j+13);                                             \
    ONESTEP(xmm14,j+14);                                             \
    ONESTEP(xmm15,j+15);                                             \
                                                                     \
    /* align at lsb */                                               \
    xmm8  = _mm_alignr_epi8(xmm8, xmm8, N);                          \
    xmm9  = _mm_alignr_epi8(xmm9, xmm9, N);                          \
    xmm10 = _mm_alignr_epi8(xmm10,xmm10,N);                          \
    xmm11 = _mm_alignr_epi8(xmm11,xmm11,N);                          \
    xmm12 = _mm_alignr_epi8(xmm12,xmm12,N);                          \
    xmm13 = _mm_alignr_epi8(xmm13,xmm13,N);                          \
    xmm14 = _mm_alignr_epi8(xmm14,xmm14,N);                          \
    xmm15 = _mm_alignr_epi8(xmm15,xmm15,N);                          \
                                                                     \
    /* interleave bytes */                                           \
    xmm8  = _mm_unpacklo_epi8(xmm8, xmm9);                           \
    xmm9  = _mm_unpacklo_epi8(xmm10,xmm11);                          \
    xmm10 = _mm_unpacklo_epi8(xmm12,xmm13);                          \
    xmm11 = _mm_unpacklo_epi8(xmm14,xmm15);                          \
                                                                     \
    /* interleave pairs */                                           \
    xmm12 = _mm_unpacklo_epi16(xmm4, xmm5);                          \
    xmm13 = _mm_unpacklo_epi16(xmm6, xmm7);                          \
    xmm14 = _mm_unpacklo_epi16(xmm8, xmm9);                          \
    xmm15 = _mm_unpacklo_epi16(xmm10,xmm11);                         \
                                                                     \
    /* interleave quadruplets */                                     \
    xmm5 = _mm_unpacklo_epi32(xmm12,xmm13);                          \
    xmm6 = _mm_unpacklo_epi32(xmm14,xmm15);                          \
                                                                     \
    /* interleave octets */                                          \
    xmm7 = _mm_unpacklo_epi64(xmm5,xmm6);                            \
                                                                     \
    /* store horizontal */                                           \
    _mm_store_si128((__m128i *)(ee+j),xmm7);                         \
  }

static void donormal8(BYTE * dseq, BYTE * qseq,
                      long dlen,long qlen)
{
  long qlen_padded = roundup(qlen,16);
  long dlen16 = (dlen >> 4) << 4;
  char c;

  __m128i xmm0, xmm1, xmm2,  xmm3,  xmm4,  xmm5,  xmm6,   xmm7;
  __m128i xmm8, xmm9, xmm10, xmm11, xmm12, xmm13, xmm14, xmm15;

  xmm0 = _mm_setzero_si128();

  for (long i = 0; i < qlen_padded; i += 16)
  {
    _mm_store_si128((__m128i *)(hh + i), xmm0);
  }

  switch (15 - (qlen_padded - qlen))
  {
    case 0:
      PROCESS16_NA;
      break;
    case 1:
      PROCESS16(1);;
      break;
    case 2:
      PROCESS16(2);
      break;
    case 3:
      PROCESS16(3);
      break;
    case 4:
      PROCESS16(4);
      break;
    case 5:
      PROCESS16(5);
      break;
    case 6:
      PROCESS16(6);
      break;
    case 7:
      PROCESS16(7);
      break;
    case 8:
      PROCESS16(8);
      break;
    case 9:
      PROCESS16(9);
      break;
    case 10:
      PROCESS16(10);
      break;
    case 11:
      PROCESS16(11);
      break;
    case 12:
      PROCESS16(12);
      break;
    case 13:
      PROCESS16(13);
      break;
    case 14:
      PROCESS16(14);
      break;
    case 15:
      PROCESS16(15);
      break;
    default:
      exit(0);
  }
  for (long j = dlen16; j < dlen; ++j)
  {
    ONESTEP(xmm4,j);
    *(ee+j) = *(hh+qlen-1);
  }

}


void salt_overlap_nuc4_sse_8(BYTE * dseq, BYTE * dend,
                             BYTE * qseq, BYTE * qend,
                             char * score_matrix,
                             long * psmscore,
                             long * overlaplen,
                             long * matchcase)
{
  long len = 0;
  char score = 0;
  long dlen = dend - dseq;
  long qlen = qend - qseq;
  long qlen_padded = roundup(qlen,16);
  char c;

  __m128i xmm0, X, H, T1, xmm1;

  xmm0 = _mm_setzero_si128();

  if (qlen_padded > hh_len)
  {
    free(hh);
    hh = xmalloc(qlen_padded*sizeof(char), SALT_ALIGNMENT_SSE);
    hh_len = qlen_padded;
  }
  if (dlen > ee_len)
  {
    free(ee);
    ee = xmalloc(roundup(dlen,16)*sizeof(char), SALT_ALIGNMENT_SSE);
    ee_len = dlen;
  }

  for (long i = 0; i < qlen_padded; i += 16)
  {
    _mm_store_si128((__m128i *)(hh + i), xmm0);
  }

  char * lastbyte= hh+qlen-1;

  qprofile_fill8_sse(score_matrix,
                     qseq,
                     qend);

  for (long j = 0; j < dlen; ++j)
  {
    X = xmm0;
    c = dseq[j];
    for (long i = 0; i < qlen_padded; i += 16)
     {
       H  = _mm_load_si128((__m128i *)(hh+i));

       T1 = _mm_srli_si128(H,15);
       H  = _mm_slli_si128(H,1);
       H  = _mm_or_si128(H,X);
       X  = T1;

       xmm1 = _mm_load_si128((__m128i *)(qprofile+c*qlen_padded+i));
       H = _mm_add_epi8(H,xmm1);

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

void salt_overlap_nuc4_sse2_8(BYTE * dseq, BYTE * dend,
                              BYTE * qseq, BYTE * qend,
                              char * score_matrix,
                              long * psmscore,
                              long * overlaplen,
                              long * matchcase)
{
  long len = 0;
  char score = 0;
  long dlen = dend - dseq;
  long qlen = qend - qseq;
  long qlen_padded = roundup(qlen,16);
  long dlen16 = (dlen >> 4) << 4;

  if (qlen_padded > hh_len)
  {
    free(hh);
    hh = xmalloc(qlen_padded*sizeof(char), SALT_ALIGNMENT_SSE);
    hh_len = qlen_padded;
  }
  if (dlen > ee_len)
  {
    free(ee);
    ee = xmalloc(roundup(dlen,16)*sizeof(char), SALT_ALIGNMENT_SSE);
    ee_len = dlen;
  }

  qprofile_fill8_sse(score_matrix,
                     qseq,
                     qend);

  donormal8(dseq, qseq,
            dlen, qlen);


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
  for (long i = 0; i < dlen16; ++i)
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
