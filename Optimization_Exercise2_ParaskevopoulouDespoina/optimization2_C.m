%% Second Exercise in Optimization
%% C
close all;
clear all;
n = 50;
m = 200;
A = randn(m,n);
c = randn(n,1);
b = abs(randn(m,1));

% we define f as
f = @(x) c' * x - sum(log(b - A * x));

% a
cvx_begin
variable x(n)
minimize((c')*x - sum(log(b-(A)*x)))
subject to
b - A*x >= 0;
cvx_end

x_cvx_optimal = x;

% b
if (n==2)
    x1 = -x_cvx_optimal-2:0.01:x_cvx_optimal+2;
    x2 = -x_cvx_optimal-2:0.01:x_cvx_optimal+2;
    [X1,X2] = meshgrid(x1,x2);
    domf = all((b-A*x) > 0);
    %x_v = [x1; x2];
    fc = zeros(size(x1,2),size(x2,2));
    for i=1:size(x1,2)
        for j = 1:size(x2,2)
            if (b - A * [x1(i) ;x2(j)] > 0)
                fc(i, j) = c' * [x1(i) ; x2(j)] - sum(log(b - A * [x1(i);x2(j)]));
            else
                fc(i,j) = 1e3;
            end
        end
    end

    figure;
    contour(X1, X2, fc');
    xlabel('x1');
    ylabel('x2');
    title('Level Sets of f(x)');
    grid on;
end

% c
gradf = @(x) c + sum(A./(b - (A)*x),1)';
x0 = zeros(n,1);
t = 1;
a = 0.2;
beta = 0.5;
e = 0.001;
[x,fun_val,f_val,all_x_kb,k] = gradient_method_backtracking(f,gradf,x0,t,a,beta,e,A,b,n,c);

% d
%hessf = @(x) sum((-1/(b - A*x).^2)*A*A');
hessf = h(A,b,c,x);
x0 = zeros(n,1);
t = 1;
a = 0.2;
beta = 0.5;
e = 0.001;
[xn,fun_valn,f_valn,all_x_kn,itern] = gradient_method_newton(f,gradf,hessf,x0,t,a,beta,e,A,b,n,c);

% e
k_axis = 1:k-1;
itern_axis = 1:itern-1;
quant_back_line = (f_val(1:k-1) - f(x_cvx_optimal));
quant_newton_line = (f_valn(1:itern-1) - f(x_cvx_optimal));

figure;
semilogy(k_axis, quant_back_line);
grid on;
title('$(f(\mathbf{x}_k) - p_*)$','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

figure;
semilogy(itern_axis, quant_newton_line);
grid on;
title('$(f(\mathbf{x}_k) - p_*)$','Interpreter','latex');
xlabel('Iteration');
ylabel('(f_k - f_{cvx})');

%% Here start the functions we implement for this exercise

% function for the computing of the hessian
function hessian = h(A, b, c, x)

    hf = zeros(size(c,1),size(c,1));

    for i = 1:size(b,1)
        hf = hf + (1/(b(i)-A(i,:)*x)^2)*(A(i,:)'*A(i,:));
    end

    hessian = hf;

end

% function to check if a point belongs to domf
function domain_check = dom_f(A, b, c, x)

domain_check = true;

for i = 1:size(b,1)
    if(b(i) - A(i,:)*x <= 0)
        domain_check = false;
        break;
    end
end

end

% modified backtracking function
function [x,fun_val,f_val,all_x_kb,iter]=gradient_method_backtracking(f,g,x0,s,alpha,beta,epsilon,A,b,n,c)
x=x0;
grad =  g(x);
fun_val= f(x);
iter=0;
all_x_kb = zeros(n,iter);
f_val = zeros(1, iter);
ft = 1e3;
domain_check = true;

while(norm(grad) > epsilon)

    iter=iter+1;
    t = s;

    domain_check = dom_f(A,b,c,x);
    % that was proved unecessary
    while (domain_check == false)
        t = beta * t;
        x = x -t*grad;
        domain_check = dom_f(A,b,c,x);
    end

    % domain_check = dom_f(A,b,c,x-t*grad);

    while (fun_val - ft < alpha*t*norm(grad)^2)
        if (dom_f(A,b,c,x-t*grad) == false)

            ft = 1e3;
        else
            ft = f(x - t*grad);
        end

        t = beta*t;

    end

    x=x-t*grad;
    all_x_kb(:,iter) = x;
    f_val(iter) = f(x);
    fun_val=f(x);
    grad= g(x);
end
%fprintf('Grad = %f %f \n', grad);
fprintf('Optimal value f(x) = %f \n', fun_val);
fprintf('Iterations needed for backtracking = %d \n',iter);
end

% Newton's method
function [x,fun_val,f_val,all_x_kb,iter]=gradient_method_newton(f,g,hf,x0,s,alpha,beta,epsilon,A,b,n,c)
x=x0;
hessian = h(A,b,c,x);
grad = g(x);
fun_val= f(x);
iter=0;
dxn = -inv(hessian)*grad;
lambda = grad'*inv(hessian)*grad;
ft = 1e3;
domain_check = true;

all_x_kb = zeros(n,iter);
f_val = zeros(1, iter);
t=s;

while (lambda > 2*epsilon)
    iter=iter+1;

    %dxn = -(inv(hessian))*grad;
    domain_check = dom_f(A,b,c,x);

    while (domain_check == false)
        t = beta * t;
        x = x + t*dxn;
        domain_check = dom_f(A,b,c,x+ t*dxn);
    end

    while (fun_val - ft >  norm(dxn)*t*grad' *grad)
        if (dom_f(A,b,c,x +t*dxn) == false)
            ft = 1e3;
        else
            ft = f(x + t*dxn);

        end
        t = beta*t;
    end

    x=x+t*dxn;

    all_x_kb(:,iter) = x;
    f_val(iter) = f(x);

    fun_val=f(x);
    grad=g(x);
    hessian = h(A,b,c,x);
    lambda=  grad'*inv(hessian)*grad;
    dxn = -inv(hessian)*grad;
end

fprintf('Optimal value f(x) = %f \n', fun_val);
fprintf('Iterations needed for Newton = %d \n',iter);
end