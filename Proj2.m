function P = Proj2(y, H, z)

y=y(:);
m = size(H,2);
P = (eye(m)-H'*H)*y+H'*z;