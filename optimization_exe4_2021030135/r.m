function residual = r(A, b, c, t, y)

[p, n] = size(A);

x = y(1:n);
lambda = y(n+1:2*n);
v = y(2*n+1:2*n+p);

r_d = c - lambda + A'*v;

r_c = diag(lambda)*x - (1/t)*ones(n,1);

r_p = A*x - b;

residual = [r_d; r_c; r_p];
