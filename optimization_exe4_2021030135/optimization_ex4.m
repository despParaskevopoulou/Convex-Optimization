%% 4th exercise in optimization
clear all;
close all;
clf;
p=1;
n=2;

A = randn(p,n);
x = abs(randn(n,1));
c = randn(n,1);
b = A*x;
d = A*x-b; 
D = diag(1./d);

%% 1
% cvx solution
cvx_begin
    variable x(n)
    minimize ((c')*x)
    subject to
        A*x==b
        x>=0
cvx_end

disp('Optimal value of x:');
x_opt = x;
disp(x_opt);
disp('Optimal value of the objective function:');
disp((c') * x_opt);

% count the number of the nonzero elements
nonzero_elements = 0;
for i = 1:length(x_opt)
    if(x_opt(i) > 10^-5)
        nonzero_elements = nonzero_elements+1;
    end
end
disp(nonzero_elements);

if(nonzero_elements < p)
    fprintf('The nonzero elements of the solution < p');
elseif (nonzero_elements == p)
    fprintf('The nonzero elements of the solution = p');
else
    fprintf('The nonzero elements of the solution > p');
end

%% 2
% cvx solution for the feasible point
cvx_begin
    variable x(n)
    minimize (0)
    subject to
        A*x==b
        x>=0
cvx_end
disp('Feasible point x:');
xf = x;
disp(xf);


%% Interior point method
x0 = xf;
t1 = 10^(0);
t2 = 10^(-7);
alpha = 0.1; 
beta = 0.7; 

fprintf('\nInterior Point Method\n');

[x_int, fun_value_int, int_rec] = interiorPointMethod(A, b, c, n, x0, x_opt, t1, t2, alpha, beta);
fprintf("**************************************************************************************\n\n")

%% 3
%% Primal Dual Algorithm
mu = 10; 
m = 1;
t1 = 10^(-5);
t2 = 10^(-7);

fprintf('\nPrimal Dual Algorithm\n');
[x_int, fun_value_int, int_rec] = primalDual(A, b, c, n, x0+randn(n,1), x_opt, t1, t2, alpha, beta,m,mu);
fprintf("*************************************************************************************\n\n")


% function to check if a point belongs to domf
function domain_check = dom_f(A, b, c, x)
domain_check = true;
for i = 1:size(b,1)
    if(b(i) - A(i,:)*x == 0)
        domain_check = false;
        break;
    end
end

end

