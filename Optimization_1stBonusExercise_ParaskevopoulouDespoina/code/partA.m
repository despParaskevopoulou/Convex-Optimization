%% 1st Bonus Exercise in Optimization
%% case a
% first we will define f
clear all;
close all;
clc;
n = 2;
U = randn(n,n);
lmin = 1;
lmax = 10;
z = lmin + (lmax - lmin) * rand(n-2,1);
eigP = [lmin;lmax;z];
L = diag(eigP);

condN = lmax/lmin;
K = condN;

P = U*U';
q = randn(n,1); 

% 1a
cvx_begin
variable x(n)
minimize (f_x(x,P,q))
subject to
x >= 0;
cvx_end
x_optimal = x;

% 2a
x0 = ones(n,1);
[x,fun_val,f_val,all_x_kb,iter]=projected_gradient(x0,P,q,n);

iter_axis = 1:iter-1;
quant_gradient = (f_val(1:iter-1) - f_x(x_optimal,P,q));
figure;
semilogy(iter_axis, quant_gradient);
grid on;
title('$(f(\mathbf{x}_k) - p_*)$','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

% 3
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

% 4
figure;
semilogy(iter_axis, quant_gradient,'r-', 'LineWidth', 2);
hold on;
semilogy(iter_n_axis, quant_acc_gradient,'b-', 'LineWidth', 2);
grid on;
title('Comparison of the semilogies','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

%% In this part we will compute the projection of xo in the set of part a
function x_projS = projection_a(x)
x_projS = max(x,0);
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
step_size = 1;
L = max(eig(P));
    while(1)
        iter = iter+1;
        xk_plus = projection_a((x-(1/L)*gradf));

        if (norm(xk_plus - x)/(norm(xk_plus)+thresh/100) < thresh)
            fun_val = f_x(x,P,q);
            break;
        else
            x = xk_plus;
            fun_val = f_x(x,P,q);
            gradf = grad_f(x,P,q);
            all_x_kb(:,iter) = x;
            f_val(iter) = f_x(x,P,q);
        end
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

while((norm(projection_a(y-step_size*gradf_y)-x))/(norm(projection_a(y-step_size*gradf_y))) > thresh)
    gradf_x = grad_f(x,P,q);
    fun_val_x = f_x(x,P,q);
    gradf_y = grad_f(y,P,q);
    fun_val_y = f_x(y,P,q);
    x_k = x;
    x = projection_a(y-step_size*gradf_y); 
    y = x + b*(x-x_k);
    iter=iter+1;
    all_x_kb(:,iter) = x;
    f_val(iter) = f_x(x,P,q);
end
fprintf('\nOptimal value f(x) = %f \n', fun_val_x);
fprintf('\nIterations needed for accelerated projected algorithm = %d \n',iter);
end