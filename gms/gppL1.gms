$ontext
    gppL1 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:09 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*8  /;
    set s(i) / 1*3  /;
    set t(i) / 6*8  /;
    set k    / 1*1  /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
          4       5       6       7       8
  1    6.00    0.00   -3.00    0.00    0.00
  2   16.00    0.00    0.00    0.00    0.00
  3    0.00   10.00    0.00   -3.00    0.00
  4    0.00    0.00    0.00  -13.00    0.00
  5    0.00    0.00   -9.00    0.00  -14.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      4   5   6   7   8
  1   1   0   1   0   0
  2   1   0   0   0   0
  3   0   1   0   1   0
  4   0   1   0   1   0
  5   0   0   1   0   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1
  1    3.00
  2    1.00
  3    2.00
  6    2.50
  7    1.75
  8    1.50 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   6     0.00
                   7     0.00
                   8     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    18.00
                   2    18.00
                   3    18.00
                   4    20.00
                   5    20.00
                   6    10.00
                   7    15.00
                   8    20.00 / ;

$include xmodel.gms