% This example is from http://www.pnas.org/cgi/doi/10.1073/pnas.0305937101
% with 11 Nodes, which are Cln3; MBF; SBF; Cln1,2; Cdh1; Swi5; Cdc20; Clb5,6; Sic1; Clb1,2; Mcm1

m = 11;  % number of nodes, could be set to less than or equal to 11

d = 2^m; % number of states
p = d+1;
if m > 11
    error('The number of nodes cannot be large than 11.');
end
% the adjacent matrix of the Boolean network
bn_adj = [
    0	0	0	0	0	0	0	0	0	0	0;
    1	0	0	0	0	0	0	0	0	-1	0;
    1	0	0	0	0	0	0	0	0	-1	0;
    0	0	1	0	0	0	0	0	0	0	0;
    0	0	0	-1	0	0	1	-1	0	-1	0;
    0	0	0	0	0	0	1	0	0	-1	1;
    0	0	0	0	0	0	0	0	0	1	1;
    0	1	0	0	0	0	-1	0	-1	0	0;
    0	0	0	-1	0	1	1	-1	0	-1	0;
    0	0	0	0	-1	0	-1	1	-1	0	1;
    0	0	0	0	0	0	0	1	0	1	0];
bn_adj = bn_adj(1:m,1:m);
theta = zeros(m,1); % threshhold for each node
self_degradation_nodes = [1 4 6 7 11]; % the node index set with self degradation
T = zeros(m,d); % initialize the transition matrix of the Boolean network
TT = zeros(m,d); % initialize the transition matrix of the logical equations
  
for j = 1:d
    jj = d-j+1;
    x = dec2bin(j-1,m) - '0';
    x1 = bn_adj*x';
    g_theta = find(x1>theta);
    l_theta = find(x1<theta);
    e_theta = find(x1==theta);
    x1(g_theta) = 1;
    x1(l_theta) = 0;
    x1(e_theta) = x(e_theta);
    T(g_theta,j) = 1;
    T(l_theta,j) = 0;
    T(e_theta,j) = x(e_theta);
    for sdn = self_degradation_nodes
        if sdn<=m && x(sdn)==1 && any(find(e_theta==sdn))
            x1(sdn) = 0;
            T(sdn,j) = 0;
        end
    end
    for i = 1:m
        if x(i) == x1(i)
            TT(i,j) = 0;
        else
            TT(i,j) = 1;
        end
    end
end

HH = zeros(2*m,d);
for i = 1:m
    HH(2*i-1,:) = 1-TT(i,:);
    HH(2*i,:) = TT(i,:);
end
zz = zeros(2*m,1);
zz(1:2:2*m) = 1;

% To accelerate the speed, we computer some matrices in advance
invHH = pinv(HH);
eHH = eye(size(HH,2)) - invHH*HH;

y = zeros(d,p);
for s = 1:p
    for i = 1:m
        x0 = rand(d,1);
        y(:,s) = y(:,s)+(eHH*x0+invHH*zz)/m;
%         y(:,s) = y(:,s)+Proj(HH,zz,x0)/m;
    end
end

fprintf('Conceptual Implementation:\n\n');
sol = BooleanVectorSearch(y,1e-10);
% if any solutions found, convert them from tensor space to binary space
if any(sol)
    num_sols = size(sol,2);
    sol_bin = zeros(num_sols,m);
    ss = find(sol==1)-d*(0:(num_sols-1))';
    for i = 1:num_sols
        sol_bin(i,:) = dec2bin(ss(i)-1,m)-'0';
    end
    fprintf('There are %d solution(s) found (each row is a solution):\n\n',num_sols);
    disp(sol_bin);
else
    fprintf('No solusions found!\n');
end
