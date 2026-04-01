%% Hessian of function f
function hessf = hess_f(x,n)
hessf = diag(1./(x.^2));
end
%{
function [hessf] = hess_f(x,n)
    [n,~] = size(x);
    hessf = zeros(n,n);
    for i=1:n
        temp = zeros(n,1);
        temp(i) = -1;
        hessf = hessf + (1/(-x(i)^2))*temp*temp';
    end
end
%}