#include "salt.h"

static unsigned char * chrmap       = NULL;
static          long * score_matrix = NULL;

/**
 *
 */
void score_chrmap_set(unsigned char * val)
{
    chrmap = val;
}

/**
 *
 */
void score_matrix_read_aa (const char * filename)
{
    // alloc mem for result array
    if (score_matrix) xfree(score_matrix);
    score_matrix = (long *) xmalloc(32*32*sizeof(long));

    // init array
    int i;
    for(i=0; i<32*32; i++) {
        score_matrix[i] = 0;
    }

    // open file with scores
    FILE * fp;
    if(!(fp = fopen(filename, "r"))) {
        fatal("Cannot open file %s for reading.", filename);
    }

    // init vars for line reading
    char * line = (char *) xmalloc(LINE_MAX * sizeof(char));
    char del[] = " \t\n";
    char * in1, * in2, * in3;
    int pos;

    // iterate over lines of file
    while (fgets(line, LINE_MAX, fp)) {
        // get tokens
        in1 = strtok(line, del);
        in2 = strtok(NULL, del);
        in3 = strtok(NULL, del);

        // sanity check
        if (!in1 || !in2 || !in3)                 continue;
        if (strlen(in1) != 1 || strlen(in2) != 1) continue;

        // put value in3 into array position corresponding to indices in1 and in2
        pos = (chrmap[(unsigned char)*in1] << 5) + chrmap[(unsigned char)*in2];
        score_matrix[pos] = atol(in3);
    }

    fclose(fp);
}

/**
 *
 */
void score_matrix_put()
{
    int i,j;
    printf("\t\t");
    for(i=1; i<32; i++) {
        for(j=32; j<127; j++) {
            if(chrmap[j] == i) {
                printf("%c", j);
                break;
            }
        }
        printf("\t");
    }

    for(i=0; i<32*32; i++) {
        if (i%32==0) {
            printf("\n");
            for(j=32; j<127; j++) {
                if(chrmap[j] == i/32) {
                    printf("%c", j);
                    break;
                }
            }
            printf("\t");
        }

        printf("%lu\t", score_matrix[i]);
    }
    printf("\n\n");
}

/**
 *
 */
inline long score_int (int d, int q)
{
    return score_matrix[(d << 5) + q];
}

/**
 *
 */
inline long score_chr (char d, char q)
{
    return score_matrix[(chrmap[(unsigned char)d] << 5) + chrmap[(unsigned char)q]];
}
