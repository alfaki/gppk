$ontext
    gppL4 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:09 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*15 /;
    set s(i) / 1*8  /;
    set t(i) /12*15 /;
    set k    / 1*6  /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
          9      10      11      12      13      14      15
  1    7.00    0.00    0.00    0.00    0.00    0.00    0.00
  2    3.00    0.00    0.00    0.00    0.00    0.00    0.00
  3    0.00    2.00    0.00    0.00    0.00    0.00    0.00
  4    0.00   10.00    0.00    0.00    0.00    0.00    0.00
  5    0.00    5.00    0.00    0.00    0.00    0.00    0.00
  6    0.00    0.00    5.00    0.00    0.00    0.00    0.00
  7    0.00    0.00    9.00    0.00    0.00    0.00    0.00
  8    0.00    0.00   11.00    0.00    0.00    0.00    0.00
  9    0.00    0.00    0.00  -16.00  -25.00  -15.00  -10.00
 10    0.00    0.00    0.00  -16.00  -25.00  -15.00  -10.00
 11    0.00    0.00    0.00  -16.00  -25.00  -15.00  -10.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      9  10  11  12  13  14  15
  1   1   0   0   0   0   0   0
  2   1   0   0   0   0   0   0
  3   0   1   0   0   0   0   0
  4   0   1   0   0   0   0   0
  5   0   1   0   0   0   0   0
  6   0   0   1   0   0   0   0
  7   0   0   1   0   0   0   0
  8   0   0   1   0   0   0   0
  9   0   1   1   1   1   1   1
 10   1   0   1   1   1   1   1
 11   1   1   0   1   1   1   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3       4       5       6
  1    1.00    6.00    4.00    0.50    5.00    9.00
  2    4.00    1.00    3.00    2.00    4.00    4.00
  3    4.00    5.50    3.00    0.90    7.00   10.00
  4    3.00    3.00    3.00    1.00    3.00    4.00
  5    1.00    2.70    4.00    1.60    3.00    7.00
  6    1.80    2.70    4.00    3.50    6.10    3.00
  7    5.00    1.00    1.70    2.90    3.50    2.90
  8    3.00    3.00    3.00    1.00    5.00    2.00
 12    3.00    3.00    3.25    0.75    6.00    5.00
 13    4.00    2.50    3.50    1.50    7.00    6.00
 14    1.50    5.50    3.90    0.80    7.00    6.00
 15    3.00    4.00    4.00    1.80    6.00    6.00 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   4     0.00
                   5     0.00
                   6     0.00
                   7     0.00
                   8     0.00
                  12     0.00
                  13     0.00
                  14     0.00
                  15     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    75.00
                   2    75.00
                   3    75.00
                   4    75.00
                   5    75.00
                   6    75.00
                   7    75.00
                   8    75.00
                   9    75.00
                  10    75.00
                  11    75.00
                  12    10.00
                  13    25.00
                  14    30.00
                  15    10.00 / ;

$include xmodel.gms