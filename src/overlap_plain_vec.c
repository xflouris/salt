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

static WORD * qprofile_fill16_avx(WORD * score_matrix_word,
                                  BYTE * qseq,
                                  BYTE * qend)
{
  WORD * qprofile;
  WORD * offset;
  long qlen = qend - qseq;
  long padded_len = roundup(qlen, 16);
  long i;

  qprofile = xmalloc(4*padded_len*sizeof(WORD), SALT_ALIGNMENT_AVX);

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
  return qprofile;
}

static WORD * qprofile_fill16_sse(WORD * score_matrix_word,
                                  BYTE * qseq,
                                  BYTE * qend)
{
  WORD * qprofile;
  WORD * offset;
  long qlen = qend - qseq;
  long padded_len = roundup(qlen, 8);
  long i;

  qprofile = xmalloc(4*padded_len*sizeof(WORD), SALT_ALIGNMENT_SSE);

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
  return qprofile;
}

#ifdef DEBUG
void pprint_avx(__m256i x)
{
  unsigned short * p = (unsigned short *) & x;

  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x", *p++);
}

void pprint_sse(__m128i x)
{
  unsigned short * p = (unsigned short *) & x;

  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x ", *p++);
  printf("%04x", *p++);
}

void pshow_sse(char * name, __m128i x)
{
  printf("%s: ", name);
  pprint_sse(x);
  printf("\n");
}

void pshow_avx(char * name, __m256i x)
{
  printf("%s: ", name);
  pprint_avx(x);
  printf("\n");
}
#endif

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
  long  offset;
  WORD * qprofile;

  char c;

  WORD * hh = xmalloc(qlen_padded*sizeof(WORD), SALT_ALIGNMENT_AVX);

  qprofile = qprofile_fill16_avx(score_matrix,
                                 qseq,
                                 qend);

  __m256i X, H, T1, xmm0, xmm1;

  xmm0 = _mm256_setzero_si256();

  for (long i = 0; i < qlen_padded; i += 16)
  {
    _mm256_store_si256((__m256i *)(hh + i), xmm0);
  }

  for (long j = 0; j < dlen; j++)
  {
    X = xmm0;
    c = dseq[j];
    for (offset = 0; offset < qlen_padded; offset += 16)
     {
       H  = _mm256_load_si256((__m256i *)(hh+offset));
       
       xmm1 = _mm256_permute2x128_si256(H,H,_MM_SHUFFLE(3,0,1,1));
       T1 = _mm256_alignr_epi8(xmm1,xmm0,0x1e);

       xmm1 = _mm256_permute2x128_si256(H, H, _MM_SHUFFLE(0,0,3,0));
       H    = _mm256_alignr_epi8(H,xmm1,16-2);

       H  = _mm256_or_si256(H,X);
       X  = T1;

       xmm1 = _mm256_load_si256((__m256i *)(qprofile+c*qlen_padded+offset));
       H = _mm256_add_epi16(H,xmm1);

       _mm256_store_si256((__m256i *)(hh+offset),H);
     }
  }

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

  *psmscore = score;
  *overlaplen = len;
  
}

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
  long  offset;
  WORD * qprofile;

  char c;

  WORD * hh = xmalloc(qlen_padded * sizeof(WORD), SALT_ALIGNMENT_SSE);

  qprofile = qprofile_fill16_sse(score_matrix,
                                 qseq,
                                 qend);

  __m128i X, H, T1, xmm0, xmm1;

  xmm0 = _mm_setzero_si128();

  for (long i = 0; i < qlen_padded; i += 8)
  {
    _mm_store_si128((__m128i *)(hh + i), xmm0);
  }

  for (long j = 0; j < dlen; ++j)
  {
    X = xmm0;
    c = dseq[j];
    for (offset = 0; offset < qlen_padded; offset += 8)
     {
       H  = _mm_load_si128((__m128i *)(hh+offset));
       
       T1 = _mm_srli_si128(H,14);
       H  = _mm_slli_si128(H,2);
       H  = _mm_or_si128(H,X);
       X  = T1;

       xmm1 = _mm_load_si128((__m128i *)(qprofile+c*qlen_padded+offset));
       H = _mm_add_epi16(H,xmm1);

       _mm_store_si128((__m128i *)(hh+offset),H);
     }
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

  *psmscore = score;
  *overlaplen = len;
  
}
