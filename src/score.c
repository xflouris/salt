#include "salt.h"

/**
 *
 */
long * score_read_file_aa (const char * filename, unsigned char * chrmap)
{
    // alloc mem for result array
    long * score_matrix = (long *) xmalloc(32*32*sizeof(long));

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
    return score_matrix;
}

/**
 *
 */
inline long score (long * score_matrix, int d, int q)
{
    return score_matrix[(d << 5) + q];
}

/**
 *
 */
long score_chrs (long * score_matrix, char d, char q, unsigned char * chrmap)
{
    return score_matrix[(chrmap[(unsigned char)d] << 5) + chrmap[(unsigned char)q]];
}

/**
 *
 */
void score_matrix_put(long * score_matrix, unsigned char * chrmap)
{
    int i,j,f;
    fprintf(stdout, "\t\t");
    for(i=1; i<32; i++) {
        f=0;
        for(j=32; j<127; j++) {
            if(chrmap[j] == i) {
                fprintf(stdout, "%c\t", j);
                f=1;
                break;
            }
        }
        if(!f) fprintf(stdout, "\t");
    }

    for(i=0; i<32*32; i++) {
        if (i%32==0) {
            fprintf(stdout, "\n");
            if (i==0) {
                fprintf(stdout, "\t");
                continue;
            }
            f=0;
            for(j=32; j<127; j++) {
                if(chrmap[j] == i/32) {
                    fprintf(stdout, "%c\t", j);
                    f=1;
                    break;
                }
            }
             if(!f) fprintf(stdout, "\t");
        }

        fprintf(stdout, "%lu\t", score_matrix[i]);
    }
    fprintf(stdout, "\n\n");
}
