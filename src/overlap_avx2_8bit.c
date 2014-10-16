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

  dseq: the database/horizontal sequence
  qseq: the query/vertical sequence

  typical scores:
  match:     1
  mismatch: -1

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

static char * qprofile   = NULL;
static char * hh         = NULL;
static char * ee         = NULL;

static long qprofile_len = 0;
static long ee_len       = 0;
static long hh_len       = 0;

//#ifdef DEBUG
void pprint_avx8(__m256i x)
{
    char * p = (char *) & x;
    for (int i = 0; i < 32; i++) {
        printf("%+d ", *p++);
    }
}

void pshow_avx8(char * name, __m256i x)
{
  printf("%s: ", name);
  pprint_avx8(x);
  printf("\n");
}

//#endif

static void qprofile_fill8_avx (char * score_matrix,
                                BYTE * qseq,
                                BYTE * qend)
{
    // get the sizes needed for storage
    long qlen       = qend - qseq;
    long padded_len = roundup(qlen, SALT_ALIGNMENT_AVX);

    // make sure qprofile is big enough
    if (padded_len > qprofile_len) {
        free (qprofile);
        qprofile     = xmalloc(4*padded_len*sizeof(char), SALT_ALIGNMENT_AVX);
        qprofile_len = padded_len;
    }

    // declare all needed register vars
    __m256i xmm0, xmm1, xmm2,  xmm3,  xmm4,  xmm5;

    // load scoring values for each letter
    // (only [31:0] are interesting, rest is garbage)
    xmm1 = _mm256_load_si256((__m256i *)(score_matrix+0));  // A
    xmm2 = _mm256_load_si256((__m256i *)(score_matrix+32)); // C
    xmm3 = _mm256_load_si256((__m256i *)(score_matrix+64)); // G
    xmm4 = _mm256_load_si256((__m256i *)(score_matrix+96)); // T

    // copy [127:0] to [255:128]
    xmm1 = _mm256_permute4x64_epi64 (xmm1, 0x44);
    xmm2 = _mm256_permute4x64_epi64 (xmm2, 0x44);
    xmm3 = _mm256_permute4x64_epi64 (xmm3, 0x44);
    xmm4 = _mm256_permute4x64_epi64 (xmm4, 0x44);

    // loop over qseq to process it, jumping a vector size per iteration
    for (long i = 0; i < padded_len; i += 32) {
        // load data of one vector size from qseq
        xmm0 = _mm256_load_si256 ((__m256i*) (qseq+i));

        // A
        xmm5 = _mm256_shuffle_epi8 (xmm1, xmm0);
        _mm256_store_si256 ((__m256i*)(qprofile+0*padded_len+i), xmm5);

        // C
        xmm5 = _mm256_shuffle_epi8 (xmm2, xmm0);
        _mm256_store_si256 ((__m256i*)(qprofile+1*padded_len+i), xmm5);

        // G
        xmm5 = _mm256_shuffle_epi8 (xmm3, xmm0);
        _mm256_store_si256 ((__m256i*)(qprofile+2*padded_len+i), xmm5);

        // T
        xmm5 = _mm256_shuffle_epi8 (xmm4, xmm0);
        _mm256_store_si256 ((__m256i*)(qprofile+3*padded_len+i), xmm5);
    }
}

void salt_overlap_avx2_8bit (BYTE * dseq,
                             BYTE * dend,
                             BYTE * qseq,
                             BYTE * qend,
                             char * score_matrix,
                             long * psmscore,
                             long * overlaplen,
                             long * matchcase)
{
    long dlen = dend - dseq;
    long qlen = qend - qseq;
    long qlen_padded = roundup(qlen,16);



    if (qlen_padded > hh_len) {
        free(hh);
        hh = xmalloc(qlen_padded*sizeof(WORD), SALT_ALIGNMENT_AVX);
        hh_len = qlen_padded;
    }
    if (dlen > ee_len) {
        free(ee);
        ee = xmalloc(roundup(dlen,8)*sizeof(WORD), SALT_ALIGNMENT_AVX);
        ee_len = dlen;
    }

    qprofile_fill8_avx(score_matrix, qseq, qend);

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

    for (long i = 0; i < qlen_padded; i += 16) {
        _mm256_store_si256((__m256i *)(hh + i), xmm0);
    }

    char c;
    WORD * lastbyte= hh+qlen-1;

    for (long j = 0; j < dlen; j++) {
        X = xmm0;
        c = dseq[j];
        for (long i = 0; i < qlen_padded; i += 16) {
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

    // prepare to pick best value
    *matchcase = 0;
    WORD score = hh[0];
    long len   = 0;

    // find best value in normal case...
    for (long i = 0; i < qlen; ++i) {
        if (hh[i] >= score) {
            len = i+1;
            score = hh[i];
        }
    }

    // ... and run through case
    for (long i = 0; i < dlen; ++i) {
        if (ee[i] >= score) {
            len = i+1;
            score = ee[i];
            *matchcase = 1;
        }
    }

    // hand over results
    *psmscore = score;
    *overlaplen = len;
}
