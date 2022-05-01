function z = Theta1(sigma)

if sigma ~= 0 && sigma ~= 1
    error('The input argument must be 0 or 1');
end

z = [1-sigma; sigma];
