function M = BooleanMatricization(fun, m)

d = 2^m;
M = zeros(2,d);
for i=1:d
    x = dec2bin(i-1,m) - '0';
    M(2,i) = fun(x);
end
M(1,:) = 1 - M(2,:);
