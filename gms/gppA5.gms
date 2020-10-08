$ontext
    gppA5 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:09 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*8  /;
    set s(i) / 1*3  /;
    set t(i) / 6*8  /;
    set k    / 1*2  /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
          4       5       6       7       8
  1    0.00    4.00   -8.00   -8.00   -8.00
  2    3.00    3.00   -9.00    0.00    0.00
  3    3.00    3.00    0.00   -9.00   -9.00
  4    0.00    0.00  -12.00  -12.00  -12.00
  5    0.00    0.00  -12.00  -12.00  -12.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      4   5   6   7   8
  1   0   1   1   1   1
  2   1   1   1   0   0
  3   1   1   0   1   1
  4   0   1   1   1   1
  5   0   0   1   1   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2
  1    1.00    6.00
  2    5.00    3.00
  3    3.00    7.00
  6    3.00    4.00
  7    2.00    2.00
  8    4.00    5.00 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   6     0.00
                   7     0.00
                   8     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    34.00
                   2    40.00
                   3    21.00
                   4    51.00
                   5    55.00
                   6    37.00
                   7    28.00
                   8    22.00 / ;

$include xmodel.gms