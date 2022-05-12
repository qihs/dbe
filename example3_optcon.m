function [c,ceq] = example3_optcon(x)
% Nonlinear constraint for the optmization problem in example 3
c = [];
ceq = x'*x - 1;
