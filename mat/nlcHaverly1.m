function [c, ceq] = nlcHaverly1(x)
c   = [+0.5*x(1)*x(7)-1.5*x(2)*x(7)-0.5*x(5)
       +1.5*x(1)*x(8)-0.5*x(2)*x(8)+0.5*x(6)];

ceq = [+x(1)*x(7)+x(1)*x(8)-x(3)
       +x(2)*x(7)+x(2)*x(8)-x(4)];
