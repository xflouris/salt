/*
 * Generate test data
 */

#include <time.h>
#include <assert.h>
#include "salt.h"

// characters and their probs to use for generating sequences
char  cmap[] = {'A',  'C',  'G',  'T'};
float cprb[] = {0.25, 0.25, 0.25, 0.25}; // sum must be 1
#define CLEN 4 // num of elements in above lists

/*
 * Returns min of the ints
 */
inline int min (int a, int b)
{
    return a < b ? a : b;
}

/*
 * Returns a random int within range [min, max).
 *
 * It makes sure that numbers are equally distributed in the range.
 */
int random_int_range (int min, int max)
{
    assert (min < max);
    int ret;
    int n = max-min;
    int end = RAND_MAX / n;
    assert (end>0); // make sure compiler does not optimize too much
    end *= n;
    while ((ret = rand()) >= end);
    return min + (ret % n);
}

/*
 * Returns a random float in range (0, 1).
 */
inline float random_float ()
{
    return (float)rand()/(float)RAND_MAX;
}

/*
 * Returns a random char from the cmap, with probs given by cprb.
 */
char random_char ()
{
    float s = 0.0;
    float r = random_float();

    // add the probs to s until they are higher than r
    for (int j = 0; j < CLEN; j++) {
        s += cprb[j];
        if (r<=s) {
            return cmap[j];
        }
    }

    // make compiler happy.
    // this case should not occur as long as the sum of cprb is 1
    return cmap[CLEN-1];
}

/*
 * Fill a string with random chars.
 * Does not make sure the string is null-terminated!
 */
void generate_sequence (char* seq, int len)
{
    for (int i = 0; i < len; i++) {
        seq[i] = random_char();
    }
}

/*
 * Fill two strings with random chars, with an overlap between them.
 * Does not make sure the strings are null-terminated!
 *
 * overlap marks the number of overlapping chars, starting from the
 * end of seq1 to the beginning of seq2, so that
 * seq2 = seq1 + len1 - overlap
 */
void generate_pair (char* seq1, int len1, char* seq2, int len2, int overlap)
{
    // there has to be an actual overlap
    assert (overlap > 0 && overlap < len1 + len2);

    // fill seq1 completely with random chars
    generate_sequence (seq1, len1);
    generate_sequence (seq2, len2);

    char* p1; // where to start copying chars from seq1
    char* p2; // where to start storing the copied chars in seq2
    int   lc; // how many chars to copy

    if (overlap < len1) { // normal case
        p1 = seq1 + len1 - overlap;
        p2 = seq2;
        lc = min (overlap, len2);
    } else { // run through case
        p1 = seq1;
        p2 = seq2 + overlap - len1;
        lc = len1 + min (0, len2 - overlap);
    }

    // copy the overlap from seq1 to seq2
    strncpy (p2, p1, lc);
}

/*
 *
 */
void generate_reads (int total_len,
                     int reads_number,
                     int reads_min_len,
                     int reads_max_len,
                     char*  total_seq,
                     char** reads
                    )
{
    // it needs to somehow be possible to cover the total seq with reads
    assert (reads_number * reads_max_len > total_len);

    generate_sequence (total_seq, total_len);

    for (int i = 0; i < reads_number; i++) {

    }
}

/*
 * Overwrite chars with new ones with a certain proability.
 */
void induce_errors (char* seq, int len, float prob)
{
    char c;
    for (int i = 0; i < len; i++) {
        if (random_float() <= prob) {
            // make sure it's actually a different char
            do {
                c = random_char();
            } while (c == seq[i]);
            seq[i] = c;
        }
    }
}
