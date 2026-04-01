%% Second Exercise in Optimization
%% B
close all;
clear all;
n = 2;
A = randn(n,n);
[U,S,V] = svd(A);

% they are the same
one = U*U';
two = (U')*U;

lmin = 10;
lmax = 1000;
z = lmin + (lmax - lmin) * rand(n-2,1);
eigP = [lmin;lmax;z];
L = diag(eigP);

condN = lmax/lmin;
K = condN;

P = U*L*U';
%eig(P)
q = randn(n,1);

% iii
gradf = @(x) P*x + q;
x_opt = -q'/P;

x_opt = x_opt'
fstar = (0.5)*(x_opt')*P*x_opt + (q')*x_opt

% iv
cvx_begin
variable x(n)
minimize((0.5)*(x')*P*x + (q')*x)
cvx_end
x_optimal = x;

% v
% exact line search
k = 1;
e = 0.0001;
a = @(x) (gradf(x)')*P*(gradf(x));
b = @(x) norm(gradf(x))^2;
f = @(x) (0.5)*(x')*P*x + (q')*x;

x_k = [1;1];
all_x_k = zeros(2,k);
f_val = zeros(1, k);

while(norm(gradf(x_k)) > e)
    dx_k = -(gradf(x_k));
    t_star = b(x_k)/a(x_k);
    all_x_k(:,k) = x_k;
    f_val(k) = f(x_k);
    x_k = x_k + t_star*dx_k;
    k = k+1;
end
fprintf('Optimal value f = %f \n',f(x_k));
fprintf('Iterations needed for exact line = %d \n',k);

% backtracking
x0 = [1;1];
step = 1;
a = 0.2;
b = 0.5;
e = 0.0001;
[x,fun_val,f_val_b,all_x_kb,iter] = gradient_method_backtracking(f,gradf,x0,step,a,b,e);

% vi
x1 = -3:0.01:3;
x2 = -3:0.01:3;
%x_v = [x1; x2];
fc = zeros(size(x1,2),size(x2,2));
for i=1:size(x1,2)
    for j = 1:size(x2,2)
        fc(i,j) = (0.5)*([x1(i) x2(j)])*P*([x1(i) x2(j)]') + (q')*([x1(i) x2(j)]');
    end
end
contour(x1, x2, fc', f_val);
hold on;
plot(all_x_k(1, :), all_x_k(2, :), '*-');
xlabel('x1');
ylabel('x2');
title('$Contour\: Plot\: of\: x_k\: and\: f(x_k)\: for\: Exact\: Line$', 'Interpreter','latex');
legend('Contour', 'x_k');
hold off;

figure;
contour(x1,x2,fc' , f_val_b);
hold on;
plot(all_x_kb(1, :), all_x_kb(2, :), '*-');
xlabel('x1');
ylabel('x2');
title('$Contour\: Plot\: of\: x_k\: and\: f(x_k)\: for\: Backtracking$', 'Interpreter','latex');
legend('Contour', 'x_k');
hold off;

% vii
k_axis = 1:k-1;
iter_axis = 1:iter-1;
quant_ex_line = log(f_val(1:k-1) - fstar);
quant_back_line = log(f_val_b(1:iter-1) - fstar);

figure;
plot(k_axis,quant_ex_line);
grid on;
title('$log(f(\mathbf{x}_k) - p_*)\: for \:exact \:line\: search $','Interpreter','latex');

figure;
plot(iter_axis,quant_back_line);
grid on;
title('$log(f(\mathbf{x}_k) - p_*)\: for \: backtracking \:line\: search $','Interpreter','latex');

% viii
k_e= ceil(log((f(x0)-f(x_optimal))/e)/-log(1-(condN)^(-1)))

%% The needed functions for the algorithms

% backtracking function
function [x,fun_val,f_val,all_x_kb,iter]=gradient_method_backtracking(f,g,x0,s,alpha,beta,epsilon)
x=x0;
grad=g(x);
fun_val=f(x);
iter=0;
all_x_kb = zeros(2,iter);
f_val = zeros(1, iter);
while (norm(grad)>epsilon)
    iter=iter+1;
    t=s;
    while (fun_val-f(x-t*grad) < (alpha*t*norm(grad)^2))
        t=beta*t;
    end
    all_x_kb(:,iter) = x;
    f_val(iter) = f(x);
    x=x-t*grad;
    fun_val=f(x);
    grad=g(x);
end
fprintf('Optimal value f = %f \n', fun_val);
fprintf('Iterations needed for backtracking = %d \n',iter);
end
