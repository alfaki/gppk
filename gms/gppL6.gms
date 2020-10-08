$ontext
    gppL6 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:09 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*8  /;
    set s(i) / 1*4  /;
    set t(i) / 7*8  /;
    set k    / 1*1  /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
          5       6       7       8
  1    6.00    0.00    0.00    0.00
  2   15.00    0.00    0.00    0.00
  3   16.00    0.00    0.00    0.00
  4    0.00   10.00    0.00    0.00
  5    0.00    0.00   -9.00  -15.00
  6    0.00    0.00   -9.00  -15.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      5   6   7   8
  1   1   0   0   0
  2   1   0   0   0
  3   1   0   0   0
  4   0   1   0   0
  5   0   1   1   1
  6   1   0   1   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1
  1    3.00
  2    1.00
  3    1.00
  4    2.00
  7    2.50
  8    1.50 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   4     0.00
                   7     0.00
                   8     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1   300.00
                   2    50.00
                   3   300.00
                   4   300.00
                   5   300.00
                   6   300.00
                   7   100.00
                   8   200.00 / ;

$include xmodel.gms