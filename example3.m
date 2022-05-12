rng(1000); % set the random seed

m = 3; % number of nodes
d = 2^m;
p = d + 1;

f1 = @(x) x(1) | x(2) | ~x(3) ;
f2 = @(x) x(1) & (~x(1)|x(2)) & (x(1)|~x(2));
f3 = @(x) x(2) & x(3);

H1 = BooleanMatricization(f1,m);
H2 = BooleanMatricization(f2,m);
H3 = BooleanMatricization(f3,m);

z1 = Theta1(1);
z2 = Theta1(0);
z3 = Theta1(0);

% normalizing these matrices for all nodes
H1norm = norm(H1);
H1 = H1/H1norm;
z1 = z1/H1norm;
H2norm = norm(H2);
H2 = H2/H2norm;
z2 = z2/H2norm;
H3norm = norm(H3);
H3 = H3/H3norm;
z3 = z3/H3norm;

HH = vertcat(H1, H2, H3);
zz = vertcat(z1, z2, z3);
H = {H1, H2, H3};
z = {z1, z2, z3};

% M = d;
c = 0.1;
L = ones(m,m);
L = diag(ones(m,1))+ c*(L - diag(m*ones(m,1)));
T = 500;
batch = 1;
K = p*batch;
sigma = 1e-2;

ExpNum = 100;
OptNum = 100;
RepNum = 3;
E = eye(d);

opts = optimoptions('fmincon','Display','off');

% To accelerate the speed, we computer some matrices in advance
eHH = cell(size(H));
Hz = cell(size(H));
for i = 1:m
    eHH{i} = eye(size(H{i},2)) - H{i}'*H{i};
    Hz{i} = H{i}'*z{i};
end

fid = fopen('example3_data.txt','a');
for sigma = [1e-1 6e-2 3e-2 1e-2 6e-3 3e-3 1e-3 6e-4 3e-4 1e-4 6e-5 3e-5 1e-5 1e-6 1e-7 1e-8]
    
t1 = 0;
t2 = 0;
t3 = 0;
tt = 0;

for expnum = 1:ExpNum
    y = zeros(T, m, K, d);
    for i = 1:T
        for j = 1:m
            for s = 1:K
                if i == 1
                    y(i,j,s,:) = rand(d,1)*2 - 1;
                else
%                     Py = Proj2(L(j,1)*y(i-1,1,s,:)+L(j,2)*y(i-1,2,s,:)+...
%                          L(j,3)*y(i-1,3,s,:),H{j},z{j})+ laprnd(d,1,0,sigma);
                     Py = eHH{j}*squeeze(L(j,1)*y(i-1,1,s,:)+L(j,2)*y(i-1,2,s,:)+...
                          L(j,3)*y(i-1,3,s,:))+ Hz{j} + laprnd(d,1,0,sigma);
                    ind1 = find(abs(Py) > 1);
                    if (~isempty(ind1))
                        Py(ind1) = sign(Py(ind1));
                    end
                    y(i,j,s,:) = Py;
                end
            end
        end
    end

    yy1 = zeros(d, K);
    yy2 = zeros(d, K);
    yy3 = zeros(d, K);
    for i = 1:K
        yy1(:,i) = y(T,1,i,:);
        yy2(:,i) = y(T,2,i,:);
        yy3(:,i) = y(T,3,i,:);
    end
    yyy1 = zeros(d, p);
    yyy2 = zeros(d, p);
    yyy3 = zeros(d, p);
    for i = 1:p
        yyy1(:,i) = sum(yy1(:,((i-1)*batch+1):((i-1)*batch+batch)),2)/batch;
        yyy2(:,i) = sum(yy2(:,((i-1)*batch+1):((i-1)*batch+batch)),2)/batch;
        yyy3(:,i) = sum(yy3(:,((i-1)*batch+1):((i-1)*batch+batch)),2)/batch;
    end
    yyy1(:,2:end) = yyy1(:,2:end) - yyy1(:,1);
    yyy2(:,2:end) = yyy2(:,2:end) - yyy2(:,1);
    yyy3(:,2:end) = yyy3(:,2:end) - yyy3(:,1);

    myfun1 = @(x) (x'*yyy1(:,9))^2 + (x'*yyy1(:,2))^2 + ...
                  (x'*yyy1(:,3))^2 + (x'*yyy1(:,4))^2 + ...
                  (x'*yyy1(:,5))^2 + (x'*yyy1(:,6))^2 + ...
                  (x'*yyy1(:,7))^2 + (x'*yyy1(:,8))^2;
    myfun2 = @(x) (x'*yyy2(:,9))^2 + (x'*yyy2(:,2))^2 + ...
                  (x'*yyy2(:,3))^2 + (x'*yyy2(:,4))^2 + ...
                  (x'*yyy2(:,5))^2 + (x'*yyy2(:,6))^2 + ...
                  (x'*yyy2(:,7))^2 + (x'*yyy2(:,8))^2;
    myfun3 = @(x) (x'*yyy3(:,9))^2 + (x'*yyy3(:,2))^2 + ...
                  (x'*yyy3(:,3))^2 + (x'*yyy3(:,4))^2 + ...
                  (x'*yyy3(:,5))^2 + (x'*yyy3(:,6))^2 + ...
                  (x'*yyy3(:,7))^2 + (x'*yyy3(:,8))^2;

    L1 = zeros(d, 1);
    L2 = zeros(d, 1);
    L3 = zeros(d, 1);

    xopt = zeros(m, RepNum, OptNum, d);
    fval = zeros(m, RepNum, OptNum);
    
    opt_ind = zeros(m,RepNum);
    for r = 1:RepNum
        for i = 1:OptNum
            x0 = rand(d,1)*2-1; 
            x0 = x0/norm(x0);
            [xopt(1,r,i,:),fval(1,r,i)] = fmincon(myfun1,x0,[],[],[],[],[],[],'example3_optcon',opts);
            [xopt(2,r,i,:),fval(2,r,i)] = fmincon(myfun2,x0,[],[],[],[],[],[],'example3_optcon',opts);
            [xopt(3,r,i,:),fval(3,r,i)] = fmincon(myfun3,x0,[],[],[],[],[],[],'example3_optcon',opts);
        end
        for ii = 1:m
            opt_ind(ii,r) = find(fval(ii,r,:)==min(fval(ii,r,:)),1);
        end
    end
    
    for l = 1:d
        ll = E(:,l)-yyy1(:,1);
        ll = ll/norm(ll);
        for r = 1:RepNum
            L1(l) = L1(l) + (ll'*squeeze(xopt(1,r,opt_ind(1,r),:)))^2; 
            L2(l) = L2(l) + (ll'*squeeze(xopt(2,r,opt_ind(2,r),:)))^2; 
            L3(l) = L3(l) + (ll'*squeeze(xopt(3,r,opt_ind(3,r),:)))^2; 
        end
    end
    r1 = find(L1==min(L1));
    r2 = find(L2==min(L2));
    r3 = find(L3==min(L3));
%     disp([r1 r2 r3]);
    if find(r1==[1 3 5 6]), inc1 = 1; else inc1 = 0; end
    if find(r2==[1 3 5 6]), inc2 = 1; else inc2 = 0; end
    if find(r3==[1 3 5 6]), inc3 = 1; else inc3 = 0; end
    t1 = t1 + inc1;
    t2 = t2 + inc2;
    t3 = t3 + inc3;
    if inc1==1 && inc2==1 && inc3==1
        tt = tt+1;
    end
    text = sprintf('%s sigma:%1.0e, T=%d, batch=%d, optnum=%d, t1:%d, t2:%d, t3:%d, tt:%d, expnum:%d\n',datestr(now),sigma,T, batch,OptNum,t1,t2,t3,tt,expnum);
    disp(text);
end % end for expnum
fwrite(fid, text);
end % end for sigma
fclose(fid);
