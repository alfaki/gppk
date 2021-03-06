$ontext
    gppC4 pooling problem data.
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
  1    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -6.00    0.00    0.00
  2    0.00    5.00    5.00    0.00    5.00    0.00    0.00    0.00    0.00    0.00    0.00   -6.00
  3    5.00    5.00    5.00    5.00    5.00    5.00    0.00    0.00    0.00   -1.00    0.00   -6.00
  4    0.00    2.00    2.00    2.00    0.00    2.00   -9.00    0.00   -9.00    0.00    0.00    0.00
  5    0.00    0.00    1.00    1.00    1.00    1.00    0.00    0.00  -10.00    0.00    0.00    0.00
  6    3.00    0.00    0.00    0.00    3.00    3.00    0.00    0.00    0.00   -3.00    0.00    0.00
  7    0.00    0.00    1.00    0.00    1.00    0.00    0.00   -4.00  -10.00   -5.00    0.00  -10.00
  8    1.00    0.00    1.00    0.00    1.00    1.00    0.00    0.00  -10.00   -5.00    0.00    0.00
  9    0.00    0.00    0.00    0.00    0.00    0.00  -11.00   -5.00  -11.00   -6.00  -11.00    0.00
 10    0.00    0.00    0.00    0.00    0.00    0.00  -11.00   -5.00  -11.00   -6.00  -11.00  -11.00
 11    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -5.00    0.00   -6.00  -11.00    0.00
 12    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -5.00    0.00   -6.00  -11.00  -11.00
 13    0.00    0.00    0.00    0.00    0.00    0.00  -11.00    0.00    0.00   -6.00  -11.00    0.00
 14    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -5.00  -11.00   -6.00    0.00  -11.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      9  10  11  12  13  14  15  16  17  18  19  20
  1   1   1   1   1   1   1   0   0   0   1   0   0
  2   0   1   1   0   1   0   0   0   0   0   0   1
  3   1   1   1   1   1   1   0   1   0   1   0   1
  4   0   1   1   1   0   1   1   0   1   0   0   0
  5   0   0   1   1   1   1   0   0   1   0   0   0
  6   1   0   0   0   1   1   0   0   0   1   0   0
  7   0   0   1   0   1   0   0   1   1   1   0   1
  8   1   0   1   0   1   1   0   0   1   1   0   0
  9   0   0   1   1   0   1   1   1   1   1   1   0
 10   0   0   1   1   1   1   1   1   1   1   1   1
 11   0   0   0   1   0   1   0   1   0   1   1   0
 12   0   0   0   0   1   0   0   1   0   1   1   1
 13   0   0   0   0   0   0   1   0   0   1   1   0
 14   0   0   0   0   0   0   0   1   1   1   0   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3       4
  1    1.00    0.00    9.00    4.00
  2    3.00    1.00    9.00    8.00
  3    0.00    0.00    2.00    9.00
  4    7.00    6.00    2.00    1.00
  5    9.00    3.00    6.00    6.00
  6    0.00    0.00    1.00    1.00
  7    9.00    9.00    5.00    3.00
  8    4.00    7.00    3.00    0.00
 15    5.00    2.00    3.00    6.00
 16    2.00    6.00    3.00    2.00
 17    5.00    2.00    5.00    3.00
 18    6.00    6.00    2.00    2.00
 19    5.00    6.00    4.00    5.00
 20    3.00    4.00    4.00    4.00 ;

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
parameter bu(i) /  1    36.00
                   2    32.00
                   3    27.00
                   4    35.00
                   5    28.00
                   6    50.00
                   7    41.00
                   8    40.00
                   9    43.00
                  10    50.00
                  11    50.00
                  12    39.00
                  13    31.00
                  14    46.00
                  15    57.00
                  16    35.00
                  17    27.00
                  18    37.00
                  19    26.00
                  20    29.00 / ;

$include xmodel.gms