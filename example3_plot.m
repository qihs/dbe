optnum = 100;
sigma = [1.00E-01 6.00E-02 3.00E-02 1.00E-02 6.00E-03 3.00E-03 1.00E-03 ...
         6.00E-04 3.00E-04 1.00E-04 6.00E-05 3.00E-05 1.00E-05 1.00E-06 ...
         1.00E-07 1.00E-08]';
x = 0:length(sigma)-1;
% The following data are from example3_data.txt which is obtained by
% running example3.m
t1 = [72 87 88 97 97 98 99 100 100 100 100 100 100 100 100 100]';
t2 = [66 78 85 96 98 95 99 99 100 100 100 100 100 100 100 100]';
t3 = [61 72 81 91 98 97 99 100 100 100 100 100 100 100 100 100]';
tt = [28 53 55 85 93 90 97 99 100 100 100 100 100 100 100 100]';

plot(x, t1/optnum, 'o-', x, t2/optnum, 's-', x, t3/optnum, 'd-', x, tt/optnum, '*-','LineWidth',1);
title('Rates of Correctness for T=500');
grid on;
ylim([0 1]);
xlabel('\sigma');
ylabel('Rate of Correctness');
hl = legend('Rate of Correctness for Node 1', 'Rate of Correctness for Node 2', ...
    'Rate of Correctness for Node 3', 'Rate of Correctness for All Nodes');
set(hl,'Interpreter','latex','location','best');
xticks(x);
xticklabels({'10^{-1}','6\times10^{-2}','3\times10^{-2}','10^{-2}','6\times10^{-3}',...
    '3\times10^{-3}','10^{-3}','6\times10^{-4}','3\times10^{-4}','10^{-4}',...
    '6\times10^{-5}','3\times10^{-5}','10^{-5}','10^{-6}','10^{-7}','10^{-8}'});
xtickangle(45);


