$ontext
    gppA1 pooling problem data.
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
  1    5.00    5.00    0.00    0.00   -8.00
  2    0.00    0.00    0.00  -10.00    0.00
  3    3.00    3.00    0.00   -7.00    0.00
  4    0.00    0.00    0.00  -10.00  -13.00
  5    0.00    0.00  -13.00  -10.00  -13.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      4   5   6   7   8
  1   1   1   0   0   1
  2   1   1   0   1   0
  3   1   1   0   1   0
  4   0   1   0   1   1
  5   0   0   1   1   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2
  1    6.00    8.00
  2    1.00    2.00
  3    3.00    0.00
  6    5.00    6.00
  7    5.00    6.00
  8    6.00    5.00 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   6     0.00
                   7     0.00
                   8     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    22.00
                   2    50.00
                   3    45.00
                   4    35.00
                   5    59.00
                   6    50.00
                   7    31.00
                   8    35.00 / ;

$include xmodel.gms