%% function for the computation of the gradient
function gradf = grad_f(x,c,t,n)
gradf = t*c - 1./x;
end