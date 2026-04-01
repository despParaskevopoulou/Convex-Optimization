%% part D
close all;
clear all;
clc;
n = 20;
U = randn(n,n);
P = U*U';
q = randn(n,1); 

% 1d
c = 2;
cvx_begin
variable x(n)
minimize (f_x(x,P,q))
subject to
norm(x) <= c  ;
cvx_end
x_optimal = x;

% 2d
x0 = 4*ones(n,1);
[x,fun_val,f_val,all_x_kb,iter]=projected_gradient(x0,P,q,n);

iter_axis = 1:iter-1;
quant_gradient = (f_val(1:iter-1) - f_x(x_optimal,P,q));
figure;
semilogy(iter_axis, quant_gradient);
grid on;
title('$(f(\mathbf{x}_k) - p_*)$','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

% 3d
x0_ = ones(n,1);
y0 = x0_;
[x,fun_val_x,f_val_n,all_x_kb_n,iter_n]=accelerated_projected_gradient(x0_,y0,P,q,n);

iter_n_axis = 1:iter_n-1;
quant_acc_gradient = (f_val_n(1:iter_n-1) - f_x(x_optimal,P,q));
figure;
semilogy(iter_n_axis, quant_acc_gradient);
grid on;
title('$(f(\mathbf{x}_k) - p_*)$','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

% 4d
figure;
semilogy(iter_axis, quant_gradient,'r-', 'LineWidth', 1.3);
hold on;
semilogy(iter_n_axis, quant_acc_gradient,'b-', 'LineWidth', 1.3);
grid on;
title('Comparison of the semilogies','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

%% projection
function x_projS = projection_d(x)
x_projS = x;
c = 2;
if(norm(x)>c)
    x_projS = (c/norm(x))*x;
elseif (norm(x)<=c)
    x_projS = x;
end
end

%% function for projected gradient algorithm
function [x,fun_val,f_val,all_x_kb,iter]=projected_gradient(x0,P,q,n)
x=x0;
gradf = grad_f(x0,P,q);
fun_val=f_x(x0,P,q);
iter=0;
all_x_kb = zeros(n,iter);
f_val = zeros(1, iter);
thresh = 0.000001;
step_size = 1/max(eig(P));
while((norm(projection_d(x-step_size*gradf)-x))/(norm(projection_d(x-step_size*gradf))) > thresh)
    gradf = grad_f(x,P,q);
    fun_val = f_x(x,P,q);
    x = projection_d(x-step_size*gradf); 
    iter=iter+1;
    all_x_kb(:,iter) = x;
    f_val(iter) = f_x(x,P,q);
end
fprintf('\nOptimal value f(x) = %f \n', fun_val);
fprintf('\nIterations needed for projected algorithm = %d \n',iter);
end

%% function for accelerated projected gradient algorithm
function [x,fun_val_x,f_val,all_x_kb,iter]=accelerated_projected_gradient(x0,y0,P,q,n)
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
b = (sqrt(L)-sqrt(l))/(sqrt(L)+sqrt(l));

while((norm(projection_d(y-step_size*gradf_y)-x))/(norm(projection_d(y-step_size*gradf_y))) > thresh)
    gradf_x = grad_f(x,P,q);
    fun_val_x = f_x(x,P,q);
    gradf_y = grad_f(y,P,q);
    fun_val_y = f_x(y,P,q);
    x_k = x;
    x = projection_d(y-step_size*gradf_y); 
    y = x + b*(x-x_k);
    iter=iter+1;
    all_x_kb(:,iter) = x;
    f_val(iter) = f_x(x,P,q);
end
fprintf('\nOptimal value f(x) = %f \n', fun_val_x);
fprintf('\nIterations needed for accelerated projected algorithm = %d \n',iter);
end

