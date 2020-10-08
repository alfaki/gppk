$ontext
    gppB1 pooling problem data.
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
  1    1.00    1.00    0.00    1.00    0.00    0.00  -11.00
  2    0.00    0.00    0.00    0.00    0.00   -9.00    0.00
  3    4.00    4.00    4.00    4.00    0.00    0.00    0.00
  4    0.00    0.00    0.00    5.00    0.00    0.00    0.00
  5    0.00    3.00    3.00    3.00   -3.00    0.00    0.00
  6    0.00    0.00    0.00    0.00   -6.00    0.00  -12.00
  7    0.00    0.00    0.00    0.00    0.00    0.00    0.00
  8    0.00    0.00    0.00    0.00   -6.00   -9.00  -12.00
  9    0.00    0.00    0.00    0.00   -6.00   -9.00    0.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      6   7   8   9  10  11  12
  1   1   1   0   1   0   0   1
  2   1   0   1   1   0   1   0
  3   1   1   1   1   0   0   0
  4   0   0   0   1   0   0   0
  5   0   1   1   1   1   0   0
  6   0   1   0   1   1   0   1
  7   0   0   0   1   0   0   0
  8   0   0   0   0   1   1   1
  9   0   0   0   0   1   1   0 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3
  1    8.00    9.00    4.00
  2    4.00    4.00    7.00
  3    7.00    8.00    2.00
  4    6.00    6.00    0.00
  5    3.00    6.00    2.00
 10    2.00    6.00    6.00
 11    5.00    6.00    4.00
 12    5.00    3.00    2.00 ;

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
parameter bu(i) /  1    48.00
                   2    45.00
                   3    53.00
                   4    45.00
                   5    58.00
                   6    45.00
                   7    21.00
                   8    45.00
                   9    54.00
                  10    38.00
                  11    58.00
                  12    52.00 / ;

$include xmodel.gms