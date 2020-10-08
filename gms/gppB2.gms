$ontext
    gppB2 pooling problem data.
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
  1    0.00    4.00    4.00    4.00    0.00   -2.00    0.00
  2    0.00    2.00    2.00    2.00   -3.00   -4.00   -5.00
  3    5.00    5.00    5.00    5.00    0.00   -1.00   -2.00
  4    0.00    0.00    3.00    3.00   -2.00    0.00   -4.00
  5    0.00    0.00    0.00    0.00    0.00   -6.00    0.00
  6    0.00    0.00    0.00    0.00   -5.00   -6.00   -7.00
  7    0.00    0.00    0.00    0.00    0.00   -6.00    0.00
  8    0.00    0.00    0.00    0.00   -5.00    0.00   -7.00
  9    0.00    0.00    0.00    0.00   -5.00   -6.00   -7.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      6   7   8   9  10  11  12
  1   0   1   1   1   0   1   0
  2   0   1   1   1   1   1   1
  3   1   1   1   1   0   1   1
  4   0   0   1   1   1   0   1
  5   1   0   1   1   0   1   0
  6   0   0   0   1   1   1   1
  7   0   0   1   1   0   1   0
  8   0   0   0   1   1   0   1
  9   0   0   0   0   1   1   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3
  1    5.00    2.00    3.00
  2    5.00    8.00    8.00
  3    7.00    9.00    0.00
  4    8.00    1.00    3.00
  5    6.00    5.00    1.00
 10    2.00    5.00    3.00
 11    3.00    5.00    2.00
 12    6.00    5.00    2.00 ;

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
parameter bu(i) /  1    22.00
                   2    48.00
                   3    34.00
                   4    41.00
                   5    51.00
                   6    41.00
                   7    26.00
                   8    51.00
                   9    38.00
                  10    46.00
                  11    24.00
                  12    30.00 / ;

$include xmodel.gms