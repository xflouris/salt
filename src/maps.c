#include "salt.h"

/*
  legal symbols: *abcdefghiklmnpqrstuvxyz (all except j and o), also upper case
  fatal symbols: .-
  fatal: ascii 0-26 except tab (9), newline (10 and 13), vt (11), formfeed (12)
  stripped: !"#$&'()+,/0123456789:;<=>?@JO^_`jo�����ŧ�� as well as chrs 9-13
  
  includes both amino acid and nucleotide sequences, adapt to nt only
*/

unsigned int chrstatus[256] =
  {
    /*
    0=stripped, 1=legal, 2=fatal, 3=silently stripped
    @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
    P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

    2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,
    0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,
    0,  1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,
    0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
  };

unsigned int chrmap_2bit[256] =
  {
    /* 
       map from ascii to 2-bit code
       A and all others: 0
       C: 1
       G: 2
       TU: 3

    @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
    P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _              
    */

    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
  };


unsigned int chrmap_4bit[256] =
  {
    /*
      map from ascii to 4-bit code

      A: 1
      C: 2
      G: 3
      TU: 4
      R: 5
      Y: 6
      S: 7
      W: 8
      K: 9
      M: 10
      B: 11
      D: 12
      H: 13
      V: 14
      N: 15
      other: 0

     @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
     P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */
    
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  1, 11,  2, 12,  0,  0,  3, 13,  0,  0,  9,  0, 10, 15,  0,
     0,  0,  5,  7,  4,  4, 14,  8,  0,  6,  0,  0,  0,  0,  0,  0,
     0,  1, 11,  2, 12,  0,  0,  3, 13,  0,  0,  9,  0, 10, 15,  0,
     0,  0,  5,  7,  4,  4, 14,  8,  0,  6,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
  };


char chrmap_complement[256] =
  {
    /*
     map from ascii to ascii, complementary nucleotide
     @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
     P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','T','V','G','H','N','N','C','D','N','N','M','N','K','N','N',
    'N','N','Y','S','A','A','B','W','N','R','N','N','N','N','N','N',
    'N','t','v','g','h','N','N','c','d','N','N','m','N','k','n','N',
    'N','N','y','s','a','a','b','w','N','r','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
  };
