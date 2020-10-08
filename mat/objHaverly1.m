function f = objHaverly1(x)
if strcmp(x,'init')
    f.LB = zeros(1,8);
    f.UB = [1 1 300 300 100 200 100 200];
    i = [1 2 3 3 4 4 5 5 6 6];
    j = [3 4 5 6 3 4 5 7 6 8];
    s = [1 1 1 1 1 1 1 1 1 1];
    f.Aineq = full(sparse(i,j,s,6,8));
    f.bineq = [300 300 300 300 100 200]';
    i = [1 1];
    j = [1 2];
    s = [1 1];
    f.Aeq = full(sparse(i,j,s,1,8));
    f.beq = [1]';
    f.fitnessfcn = @objHaverly1;
    f.nonlcon = 'nlcHaverly1';
    f.nvars = 8;
    f.options.PopulationSize = 50;
    f.options.Generations = 10;
    f.options.ConstrBoundary = 'absorb';
else
    f = +6*x(3)+16*x(4)+1*x(5)-5*x(6)-9*x(7)-15*x(8);
end
