/*
 * Generate test data
 */

#include <time.h>
#include <assert.h>
#include "salt.h"

char  cmap[] = {'A',  'C',  'G',  'T'};
float cprb[] = {0.25, 0.25, 0.25, 0.25};
#define CLEN 4

/*
 * Returns a random float.
 */
inline float random_float ()
{
    return (float)rand()/(float)RAND_MAX;
}

/*
 * Returns a random char from the cmap, with probs given by cprb.
 */
inline char random_char ()
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

    // get rid of compiler warning.
    // this case should not occur as long as the sum of cprb is 1
    return cmap[CLEN-1];
}

/*
 * Fill a string with random chars.
 * Does not make sure the string is null-terminated!
 */
void generate_sequence (char* seq, int len)
{
    srand(time(NULL));
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

    char* p1; // where to start copying from seq1 to seq2
    char* p2; // where to start storing the copied chars
    int   lc; // how many chars to copy

    char* pf; // where to start filling seq2 with randomness
    int   lf; // how many chars to fill

    if (overlap < len1) { // normal case
        p1 = seq1 + len1 - overlap;
        p2 = seq2;
        lc = overlap;
        pf = seq2 + overlap;
        lf = len2 - lc;
    } else { // run through case
        p1 = seq1;
        p2 = seq2 - len1 + overlap;
        lc = len1 + len2 - overlap;
        pf = seq2;
        lf = len2 - lc;
    }

    // copy the overlap from seq1 to seq2, then fill the rest with random chars
    strncpy (p2, p1, lc);
    generate_sequence (pf, lf);
}

/*
 * Overwrite chars with new ones with a certain proability.
 */
void induce_errors (char* seq, int len, float prob)
{
    char c;
    for (int i = 0; i < len; i++) {
        if (random_float() < prob) {
            // make sure it's actually a different char
            do {
                c = random_char();
            } while (c == seq[i]);
            seq[i] = c;
        }
    }
}
