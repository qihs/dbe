function Y = DistributedLAE(H, z, L, T)

m = size(L,1);
d = 2^m;
p = d + 1;

for i = 1:m
    Hnorm = norm(H{i});
    H{i} = H{i}/Hnorm;
    z{i} = z{i}/Hnorm;
end

Y = mat2cell(rand(d*m,p),d*ones(1,m),p); 
for s = 1:p
    for iter = 1:T
        for j = 1:m
            yy = zeros(d,1);
            for n = 1:m
                yy = yy + L(j,n)*Y{n}(:,s);
            end
            Y{j}(:,s) = Proj2(yy,H{j},z{j});
        end
    end
end
