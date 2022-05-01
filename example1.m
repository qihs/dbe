m = 3; % number of nodes
iter = 300; % iteration times
d = 2^m;
p = d + 1;

f1 = @(x) x(1) | x(2)| ~x(3) ;
f2 = @(x) x(1) & (~x(1)|x(2)) & (x(1)|~x(2));
f3 = @(x) x(2) & x(3);

H1 = BooleanMatricization(f1,m);
H2 = BooleanMatricization(f2,m);
H3 = BooleanMatricization(f3,m);

z1 = Theta1(1);
z2 = Theta1(0);
z3 = Theta1(0);


H = {H1, H2, H3};
z = {z1, z2, z3};

% Generating an adjacent matrix L of a connected graph for all the nodes
c = 1/(2*m); % a parameter less than 1/m
L = ones(m,m);
L = diag(ones(m,1))+ c*(L - diag(m*ones(m,1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distributed Linear Equation Solver
Y = DistributedLAE(H, z, L, iter);

fprintf('Distributed Linear Equation Solver: \n\n');
for k = 1:m
    sol = BooleanVectorSearch(Y{k});
    % if any solutions found, convert them from tensor space to binary space
    if any(sol)
        num_sols = size(sol,2);
        sol_bin = zeros(num_sols,m);
        ss = find(sol==1)-d*(0:(num_sols-1))';
        for i = 1:num_sols
            sol_bin(i,:) = dec2bin(ss(i)-1,m)-'0';
        end
        fprintf('There are %d solution(s) found for node %d (each row is a solution):\n\n',num_sols,k);
        disp(sol_bin);
    else
        fprintf('No solusions found for node %d!\n',k);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Centralized Linear Equation Solver
HH = vertcat(H{:});
zz = vertcat(z{:});
y = zeros(d,p);
for s = 1:p
    for i = 1:m
        x0 = rand(d,1);
        y(:,s) = y(:,s)+Proj(HH,zz,x0)/m;
    end
end
sol = BooleanVectorSearch(y);
fprintf('\nCentralized Linear Equation Solver: ');
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
