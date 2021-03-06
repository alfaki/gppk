$ontext
    gppB4 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:09 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*12 /;
    set s(i) / 1*5  /;
    set t(i) /10*12 /;
    set k    / 1*3  /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
          6       7       8       9      10      11      12
  1    1.00    1.00    1.00    1.00    0.00    0.00   -9.00
  2    2.00    2.00    2.00    2.00    0.00   -4.00   -8.00
  3    5.00    5.00    5.00    5.00    0.00   -1.00    0.00
  4    5.00    5.00    5.00    5.00   -8.00    0.00    0.00
  5    0.00    0.00    0.00    0.00  -13.00    0.00    0.00
  6    0.00    0.00    0.00    0.00  -13.00    0.00  -10.00
  7    0.00    0.00    0.00    0.00  -13.00   -6.00  -10.00
  8    0.00    0.00    0.00    0.00  -13.00   -6.00  -10.00
  9    0.00    0.00    0.00    0.00    0.00    0.00  -10.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      6   7   8   9  10  11  12
  1   1   1   1   1   0   0   1
  2   1   1   1   1   0   1   1
  3   1   1   1   1   0   1   0
  4   1   1   1   1   1   0   0
  5   1   0   1   1   1   0   0
  6   0   1   1   1   1   0   1
  7   0   0   1   1   1   1   1
  8   0   0   0   1   1   1   1
  9   0   0   0   0   0   0   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3
  1    4.00    1.00    4.00
  2    0.00    2.00    1.00
  3    5.00    8.00    5.00
  4    7.00    7.00    6.00
  5    1.00    0.00    8.00
 10    3.00    2.00    3.00
 11    6.00    6.00    3.00
 12    2.00    4.00    3.00 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   4     0.00
                   5     0.00
                  10     0.00
                  11     0.00
                  12     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    41.00
                   2    40.00
                   3    36.00
                   4    40.00
                   5    52.00
                   6    20.00
                   7    21.00
                   8    57.00
                   9    20.00
                  10    59.00
                  11    52.00
                  12    22.00 / ;

$include xmodel.gms