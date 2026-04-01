%% part C
% 1c
clear all;
close all;
clc;
n = 100;
p = 20;

A = randn(n,n);
[U,S,V] = svd(A);

L = 100;
l = 1;

z = l + (L - l)*rand(n-2,1);

eig_P = [l; L; z];

Lambda = diag(eig_P);

P = U*Lambda*U';     

q = randn(n,1)+0.0001;

A = 2*rand(p,n); 
b = 2*rand(p,1);

% 1c
cvx_begin
variable x(n)
minimize (f_x(x,P,q))
subject to
A*x == b ;
cvx_end
x_optimal = x;

% 2c
x0 = 4*ones(n,1);
[x,fun_val,f_val,all_x_kb,iter]=projected_gradient(x0,P,q,n,A,b);

iter_axis = 1:iter-1;
quant_gradient = (f_val(1:iter-1) - f_x(x_optimal,P,q));
figure;
semilogy(iter_axis, quant_gradient);
grid on;
title('$(f(\mathbf{x}_k) - p_*)$','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

% 3c
x0_ = ones(n,1);
y0 = x0_;
[x,fun_val_x,f_val_n,all_x_kb_n,iter_n]=accelerated_projected_gradient(x0_,y0,P,q,A,b,n);

iter_n_axis = 1:iter_n-1;
quant_acc_gradient = (f_val_n(1:iter_n-1) - f_x(x_optimal,P,q));
figure;
semilogy(iter_n_axis, quant_acc_gradient);
grid on;
title('$(f(\mathbf{x}_k) - p_*)$','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

% 4c
figure;
semilogy(iter_axis, quant_gradient,'r-', 'LineWidth', 1.3);
hold on;
semilogy(iter_n_axis, quant_acc_gradient,'b-', 'LineWidth', 1.3);
grid on;
title('Comparison of the semilogies','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');
legend('Projected', 'Accelerated Projected')

%% projection
function x_projS = projection_c(x,A,b,P)
    A1 = A; 
    b1 = b;
    [p, n] = size(A1);
    P1 = eye(n,n);
    F = [P1, A1'; A1, zeros(p,p)];
    u = [x; b1];
    iF = inv(F);
    x_projS = iF(1:n,:)*u;

end


%% function for projected gradient algorithm
function [x,fun_val,f_val,all_x_kb,iter]=projected_gradient(x0,P,q,n,A,b)
x=x0;
gradf = grad_f(x0,P,q);
fun_val=f_x(x0,P,q);
iter=0;
all_x_kb = zeros(n,iter);
f_val = zeros(1, iter);
thresh = 0.000001;
step_size = 1/max(eig(P));
while((norm(projection_c((x-step_size*gradf),A,b,P)-x))/(norm(projection_c((x-step_size*gradf),A,b,P))) > thresh)
    gradf = grad_f(x,P,q);
    fun_val = f_x(x,P,q);
    x = projection_c((x-step_size*gradf),A,b,P); 
    iter=iter+1;
    all_x_kb(:,iter) = x;
    f_val(iter) = f_x(x,P,q);
end
fprintf('\nOptimal value f(x) = %f \n', fun_val);
fprintf('\nIterations needed for projected algorithm = %d \n',iter);
end

%% function for accelerated projected gradient algorithm
function [x,fun_val_x,f_val,all_x_kb,iter]=accelerated_projected_gradient(x0,y0,P,q,A,b,n)
x=x0;
y = y0;
gradf_x = grad_f(x0,P,q);
fun_val_x=f_x(x0,P,q);
gradf_y = grad_f(y0,P,q);
fun_val_y=f_x(y0,P,q);
iter=0;
all_x_kb = zeros(n,iter);
f_val = zeros(1, iter);
thresh = 0.000001;
L = max(eig(P));
l = min(eig(P));
step_size = 1/L;
beta = (sqrt(L)-sqrt(l))/(sqrt(L)+sqrt(l));

while((norm(projection_c((y-step_size*gradf_y),A,b,P)-x))/(norm(projection_c((y-step_size*gradf_y),A,b,P))) > thresh)
    gradf_x = grad_f(x,P,q);
    fun_val_x = f_x(x,P,q);
    gradf_y = grad_f(y,P,q);
    fun_val_y = f_x(y,P,q);
    x_k = x;
    x = projection_c((y-step_size*gradf_y),A,b,P); 
    y = x + beta*(x-x_k);
    iter=iter+1;
    all_x_kb(:,iter) = x;
    f_val(iter) = f_x(x,P,q);
end
fprintf('\nOptimal value f(x) = %f \n', fun_val_x);
fprintf('\nIterations needed for accelerated projected algorithm = %d \n',iter);
end