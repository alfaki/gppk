$ontext
    gppE4 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:10 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*35 /;
    set s(i) / 1*10 /;
    set t(i) /21*35 /;
    set k    / 1*12 /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
         11      12      13      14      15      16      17      18      19      20      21      22      23      24      25      26      27      28      29      30      31      32      33      34      35
  1    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -6.00    0.00    0.00    0.00    0.00  -14.00   -6.00    0.00    0.00    0.00    0.00   -9.00    0.00    0.00
  2    0.00    1.00    0.00    0.00    1.00    1.00    1.00    0.00    1.00    0.00    0.00   -5.00    0.00    0.00    0.00  -11.00    0.00    0.00    0.00    0.00  -12.00    0.00    0.00    0.00    0.00
  3    4.00    4.00    0.00    0.00    4.00    0.00    0.00    0.00    0.00    4.00    0.00    0.00    0.00    0.00   -2.00    0.00    0.00   -2.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
  4    0.00    2.00    2.00    0.00    0.00    2.00    2.00    2.00    0.00    2.00    0.00    0.00    0.00    0.00   -4.00    0.00    0.00    0.00    0.00    0.00    0.00   -4.00    0.00   -5.00    0.00
  5    0.00    0.00    0.00    0.00    0.00    0.00    4.00    0.00    4.00    0.00    0.00   -2.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -3.00    0.00
  6    0.00    0.00    5.00    0.00    0.00    5.00    0.00    0.00    0.00    0.00    0.00   -1.00    0.00   -4.00   -1.00    0.00   -9.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
  7    0.00    0.00    2.00    0.00    0.00    2.00    0.00    2.00    2.00    2.00   -7.00   -4.00   -4.00   -7.00    0.00    0.00    0.00   -4.00    0.00   -3.00    0.00    0.00    0.00    0.00    0.00
  8    0.00    0.00    4.00    0.00    0.00    0.00    4.00    0.00    0.00    0.00   -5.00    0.00    0.00    0.00   -2.00    0.00    0.00   -2.00    0.00    0.00    0.00   -2.00   -5.00    0.00    0.00
  9    0.00    1.00    0.00    0.00    0.00    1.00    1.00    1.00    0.00    1.00   -8.00   -5.00    0.00    0.00    0.00    0.00  -13.00    0.00    0.00    0.00    0.00    0.00   -8.00    0.00   -8.00
 10    0.00    0.00    5.00    5.00    0.00    5.00    0.00    5.00    5.00    0.00    0.00   -1.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -1.00    0.00    0.00   -4.00
 11    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -9.00   -6.00    0.00   -9.00   -6.00  -12.00  -14.00   -6.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
 12    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -6.00    0.00   -9.00    0.00  -12.00    0.00    0.00   -5.00   -5.00    0.00   -6.00   -9.00    0.00    0.00
 13    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -6.00   -6.00   -9.00   -6.00    0.00  -14.00    0.00   -5.00    0.00  -13.00    0.00   -9.00   -7.00    0.00
 14    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -9.00   -6.00   -6.00   -9.00   -6.00  -12.00    0.00    0.00   -5.00   -5.00  -13.00   -6.00   -9.00   -7.00    0.00
 15    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -9.00   -6.00  -12.00    0.00    0.00    0.00   -5.00  -13.00   -6.00   -9.00   -7.00    0.00
 16    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -9.00    0.00    0.00    0.00   -6.00  -12.00  -14.00    0.00   -5.00   -5.00  -13.00   -6.00   -9.00   -7.00   -9.00
 17    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -9.00   -6.00   -6.00    0.00    0.00    0.00  -14.00   -6.00   -5.00    0.00    0.00    0.00    0.00    0.00   -9.00
 18    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -9.00    0.00   -6.00    0.00   -6.00    0.00    0.00   -6.00   -5.00   -5.00  -13.00    0.00   -9.00   -7.00    0.00
 19    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -6.00    0.00    0.00    0.00    0.00  -14.00   -6.00   -5.00   -5.00  -13.00   -6.00   -9.00    0.00   -9.00
 20    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -9.00   -6.00    0.00   -9.00   -6.00    0.00  -14.00   -6.00   -5.00   -5.00    0.00    0.00    0.00   -7.00    0.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
     11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35
  1   0   1   0   1   1   0   0   0   0   1   0   1   0   0   0   0   1   1   0   0   0   0   1   0   0
  2   0   1   0   0   1   1   1   0   1   0   0   1   0   0   0   1   0   0   0   0   1   0   0   0   0
  3   1   1   0   0   1   0   0   0   0   1   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0
  4   0   1   1   0   0   1   1   1   0   1   0   0   0   0   1   0   0   0   0   0   0   1   0   1   0
  5   0   0   0   0   0   0   1   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0
  6   0   0   1   0   0   1   0   0   0   0   0   1   0   1   1   0   1   0   0   0   0   0   0   0   0
  7   0   0   1   0   0   1   0   1   1   1   1   1   1   1   0   0   0   1   0   1   0   0   0   0   0
  8   0   0   1   0   0   0   1   0   0   0   1   0   0   0   1   0   0   1   0   0   0   1   1   0   0
  9   0   1   0   0   0   1   1   1   0   1   1   1   0   0   0   0   1   0   0   0   0   0   1   0   1
 10   0   0   1   1   0   1   0   1   1   0   0   1   0   0   0   0   0   0   1   0   0   1   0   0   1
 11   0   0   1   1   0   1   0   1   1   0   1   1   0   1   1   1   1   1   0   0   0   0   0   0   0
 12   0   0   0   0   0   1   1   1   1   1   0   1   0   1   0   1   0   0   1   1   0   1   1   0   0
 13   1   1   0   1   1   0   1   1   0   1   0   1   1   1   1   0   1   0   1   0   1   0   1   1   0
 14   0   1   1   0   0   1   0   1   1   1   1   1   1   1   1   1   0   0   1   1   1   1   1   1   0
 15   1   1   0   1   0   1   1   1   1   0   0   0   0   1   1   1   0   0   0   1   1   1   1   1   0
 16   0   0   1   1   0   0   1   0   0   1   1   0   0   0   1   1   1   0   1   1   1   1   1   1   1
 17   1   1   1   1   1   1   0   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   1
 18   1   1   0   1   0   0   1   0   1   0   1   0   1   0   1   0   0   1   1   1   1   0   1   1   0
 19   1   1   0   1   0   1   1   1   0   1   0   1   0   0   0   0   1   1   1   1   1   1   1   0   1
 20   0   0   0   0   0   0   0   0   0   0   1   1   0   1   1   0   1   1   1   1   0   0   0   1   0 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3       4       5       6       7       8       9      10      11      12
  1    8.00    7.00    0.00    5.00    7.00    8.00    1.00    6.00    1.00    2.00    2.00    5.00
  2    1.00    7.00    1.00    5.00    5.00    9.00    0.00    5.00    0.00    5.00    1.00    2.00
  3    2.00    2.00    2.00    4.00    5.00    5.00    3.00    5.00    4.00    2.00    5.00    0.00
  4    8.00    9.00    9.00    4.00    1.00    5.00    9.00    4.00    5.00    7.00    4.00    0.00
  5    6.00    7.00    1.00    5.00    4.00    1.00    3.00    0.00    7.00    3.00    3.00    1.00
  6    7.00    8.00    5.00    4.00    2.00    0.00    6.00    6.00    4.00    1.00    8.00    3.00
  7    3.00    4.00    4.00    3.00    3.00    5.00    1.00    3.00    7.00    0.00    0.00    2.00
  8    1.00    4.00    9.00    1.00    6.00    9.00    0.00    7.00    8.00    3.00    3.00    5.00
  9    6.00    0.00    1.00    0.00    3.00    1.00    8.00    5.00    6.00    6.00    5.00    9.00
 10    5.00    2.00    4.00    5.00    9.00    7.00    8.00    7.00    9.00    3.00    8.00    5.00
 21    4.00    3.00    2.00    6.00    6.00    4.00    2.00    2.00    4.00    5.00    2.00    3.00
 22    5.00    3.00    5.00    4.00    4.00    5.00    4.00    4.00    3.00    3.00    5.00    4.00
 23    2.00    2.00    6.00    6.00    2.00    3.00    5.00    6.00    5.00    4.00    5.00    3.00
 24    5.00    5.00    5.00    6.00    4.00    4.00    4.00    5.00    4.00    4.00    3.00    5.00
 25    3.00    6.00    4.00    5.00    4.00    4.00    4.00    6.00    2.00    2.00    3.00    6.00
 26    3.00    2.00    3.00    6.00    5.00    2.00    3.00    5.00    4.00    4.00    6.00    6.00
 27    3.00    2.00    6.00    3.00    6.00    3.00    2.00    5.00    5.00    4.00    4.00    6.00
 28    6.00    5.00    3.00    5.00    6.00    4.00    2.00    2.00    5.00    4.00    2.00    2.00
 29    6.00    4.00    6.00    3.00    4.00    4.00    6.00    5.00    4.00    3.00    4.00    2.00
 30    2.00    3.00    2.00    3.00    5.00    5.00    5.00    4.00    5.00    3.00    6.00    3.00
 31    4.00    5.00    3.00    4.00    6.00    6.00    2.00    5.00    4.00    4.00    2.00    4.00
 32    4.00    3.00    3.00    3.00    2.00    2.00    4.00    3.00    3.00    3.00    5.00    3.00
 33    3.00    3.00    5.00    2.00    5.00    6.00    6.00    2.00    3.00    6.00    3.00    5.00
 34    4.00    5.00    4.00    2.00    6.00    6.00    5.00    5.00    6.00    6.00    5.00    6.00
 35    5.00    3.00    4.00    5.00    3.00    4.00    4.00    3.00    6.00    6.00    5.00    5.00 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   4     0.00
                   5     0.00
                   6     0.00
                   7     0.00
                   8     0.00
                   9     0.00
                  10     0.00
                  21     0.00
                  22     0.00
                  23     0.00
                  24     0.00
                  25     0.00
                  26     0.00
                  27     0.00
                  28     0.00
                  29     0.00
                  30     0.00
                  31     0.00
                  32     0.00
                  33     0.00
                  34     0.00
                  35     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    48.00
                   2    49.00
                   3    30.00
                   4    27.00
                   5    21.00
                   6    45.00
                   7    31.00
                   8    29.00
                   9    29.00
                  10    20.00
                  11    20.00
                  12    28.00
                  13    29.00
                  14    23.00
                  15    37.00
                  16    51.00
                  17    20.00
                  18    41.00
                  19    52.00
                  20    40.00
                  21    26.00
                  22    35.00
                  23    38.00
                  24    45.00
                  25    26.00
                  26    22.00
                  27    56.00
                  28    30.00
                  29    36.00
                  30    26.00
                  31    57.00
                  32    33.00
                  33    58.00
                  34    37.00
                  35    24.00 / ;

$include xmodel.gms