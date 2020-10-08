$ontext
    gppA2 pooling problem data.
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
  1    2.00    2.00    0.00  -11.00    0.00
  2    0.00    0.00  -12.00  -13.00    0.00
  3    0.00    3.00   -9.00    0.00    0.00
  4    0.00    0.00  -12.00    0.00    0.00
  5    0.00    0.00  -12.00  -13.00    0.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      4   5   6   7   8
  1   1   1   0   1   0
  2   1   1   1   1   0
  3   0   1   1   0   0
  4   0   1   1   0   0
  5   0   0   1   1   0 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2
  1    3.00    1.00
  2    6.00    3.00
  3    8.00    1.00
  6    4.00    2.00
  7    4.00    4.00
  8    5.00    5.00 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   6     0.00
                   7     0.00
                   8     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    37.00
                   2    48.00
                   3    20.00
                   4    35.00
                   5    20.00
                   6    41.00
                   7    49.00
                   8    47.00 / ;

$include xmodel.gms