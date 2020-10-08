function f = objHaverly1p(x)
if strcmp(x,'init')
    f.LB = [1 0 0 0 0 0 0];
    f.UB = [3 300 300 100 200 100 200];
    i = [1 2 3 3 4 4 5 5 6 6];
    j = [2 3 4 5 2 3 4 6 5 7];
    s = [1 1 1 1 1 1 1 1 1 1];
    f.Aineq = full(sparse(i,j,s,6,7));
    f.bineq = [300 300 300 300 100 200]';
    i = [1 1 1 1];
    j = [2 3 6 7];
    s = [1 1 -1 -1];
    f.Aeq = full(sparse(i,j,s,1,7));
    f.beq = [0]';
    f.fitnessfcn = @objHaverly1p;
    f.nonlcon = 'nlcHaverly1p';
    f.nvars = 7;
    f.options.PopulationSize = 50;
    f.options.Generations = 10;
    f.options.ConstrBoundary = 'absorb';
else
    f = +6*x(2)+16*x(3)+1*x(4)-5*x(5)-9*x(6)-15*x(7);
end
