function P = Proj(G, z, x)

m = size(G,2);
invG = pinv(G);
P = (eye(m)-invG*G)*x+invG*z;