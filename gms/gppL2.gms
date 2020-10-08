$ontext
    gppL2 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:09 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*11 /;
    set s(i) / 1*5  /;
    set t(i) / 8*11 /;
    set k    / 1*4  /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
          6       7       8       9      10      11
  1    7.00    0.00    0.00    0.00    0.00    0.00
  2    3.00    0.00    0.00    0.00    0.00    0.00
  3    0.00    2.00    0.00    0.00    0.00    0.00
  4    0.00   10.00    0.00    0.00    0.00    0.00
  5    0.00    5.00    0.00    0.00    0.00    0.00
  6    0.00    0.00  -16.00  -25.00  -15.00  -10.00
  7    0.00    0.00  -16.00  -25.00  -15.00  -10.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
      6   7   8   9  10  11
  1   1   0   0   0   0   0
  2   1   0   0   0   0   0
  3   0   1   0   0   0   0
  4   0   1   0   0   0   0
  5   0   1   0   0   0   0
  6   0   1   1   1   1   1
  7   1   0   1   1   1   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3       4
  1    1.00    6.00    4.00    0.50
  2    4.00    1.00    3.00    2.00
  3    4.00    5.50    3.00    0.90
  4    3.00    3.00    3.00    1.00
  5    1.00    2.70    4.00    1.60
  8    3.00    3.00    3.25    0.75
  9    4.00    2.50    3.50    1.50
 10    1.50    5.50    3.90    0.80
 11    3.00    4.00    4.00    1.80 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   4     0.00
                   5     0.00
                   8     0.00
                   9     0.00
                  10     0.00
                  11     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1    75.00
                   2    75.00
                   3    75.00
                   4    75.00
                   5    75.00
                   6    75.00
                   7    75.00
                   8    10.00
                   9    25.00
                  10    30.00
                  11    10.00 / ;

$include xmodel.gms