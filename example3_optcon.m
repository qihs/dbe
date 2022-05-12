function [c,ceq] = example3_optcon(x)
% nonlinear constraint for the optmization prolem in example 3
c = [];
ceq = x'*x - 1;