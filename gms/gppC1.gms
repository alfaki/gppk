$ontext
    gppC1 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:09 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*20 /;
    set s(i) / 1*8  /;
    set t(i) /15*20 /;
    set k    / 1*4  /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
          9      10      11      12      13      14      15      16      17      18      19      20
  1    5.00    0.00    0.00    0.00    5.00    0.00    0.00    0.00    0.00    0.00    0.00   -7.00
  2    5.00    5.00    5.00    0.00    5.00    0.00    0.00    0.00   -5.00    0.00    0.00   -7.00
  3    0.00    0.00    5.00    5.00    5.00    5.00    0.00   -7.00    0.00    0.00    0.00    0.00
  4    0.00    2.00    0.00    2.00    0.00    2.00   -8.00    0.00    0.00    0.00    0.00  -10.00
  5    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -12.00  -10.00    0.00    0.00    0.00
  6    0.00    0.00    0.00    0.00    5.00    5.00   -5.00    0.00    0.00   -4.00    0.00    0.00
  7    0.00    0.00    0.00    0.00    4.00    0.00    0.00   -8.00   -6.00   -5.00    0.00    0.00
  8    3.00    0.00    0.00    3.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
  9    0.00    0.00    0.00    0.00    0.00    0.00  -10.00  -12.00  -10.00   -9.00    0.00  -12.00
 10    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -10.00   -9.00  -12.00    0.00
 11    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -10.00    0.00  -12.00    0.00
 12    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -12.00  -10.00    0.00  -12.00    0.00
 13    0.00    0.00    0.00    0.00    0.00    0.00  -10.00  -12.00    0.00    0.00    0.00    0.00
 14    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      9  10  11  12  13  14  15  16  17  18  19  20
  1   1   0   0   0   1   0   0   0   0   0   0   1
  2   1   1   1   0   1   0   0   0   1   0   0   1
  3   0   0   1   1   1   1   0   1   0   0   0   0
  4   0   1   0   1   0   1   1   0   0   0   0   1
  5   1   1   0   1   1   1   0   1   1   0   0   0
  6   0   0   0   0   1   1   1   0   0   1   0   0
  7   0   0   0   0   1   0   0   1   1   1   0   0
  8   1   0   0   1   0   0   0   0   0   0   0   0
  9   0   1   0   1   1   1   1   1   1   1   0   1
 10   0   0   1   0   0   1   0   0   1   1   1   0
 11   0   0   0   0   0   0   0   0   1   0   1   0
 12   0   0   0   0   0   0   0   1   1   0   1   0
 13   0   0   0   0   0   0   1   1   0   0   0   0
 14   0   0   0   0   0   0   0   0   0   0   0   0 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3       4
  1    3.00    6.00    1.00    0.00
  2    3.00    2.00    4.00    2.00
  3    1.00    0.00    4.00    6.00
  4    1.00    0.00    9.00    2.00
  5    6.00    4.00    8.00    4.00
  6    2.00    9.00    6.00    0.00
  7    3.00    1.00    1.00    3.00
  8    8.00    8.00    6.00    3.00
 15    2.00    4.00    4.00    6.00
 16    2.00    2.00    3.00    5.00
 17    3.00    5.00    5.00    4.00
 18    2.00    3.00    4.00    2.00
 19    6.00    2.00    5.00    2.00
 20    4.00    2.00    2.00    4.00 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   4     0.00
                   5     0.00
                   6     0.00
                   7     0.00
                   8     0.00
                  15     0.00
                  16     0.00
                  17     0.00
                  18     0.00
                  19     0.00
                  20     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    39.00
                   2    35.00
                   3    51.00
                   4    57.00
                   5    35.00
                   6    37.00
                   7    23.00
                   8    34.00
                   9    49.00
                  10    41.00
                  11    22.00
                  12    46.00
                  13    46.00
                  14    35.00
                  15    57.00
                  16    26.00
                  17    55.00
                  18    59.00
                  19    54.00
                  20    33.00 / ;

$include xmodel.gms