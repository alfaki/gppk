function [c, ceq] = nlcHaverly1p(x)
c   = [-0.5*x(4)+x(6)*x(1)-2.5*x(6)
       +0.5*x(5)+x(7)*x(1)-1.5*x(7)];

ceq = [+3*x(2)+1*x(3)-x(6)*x(1)-x(7)*x(1)];
