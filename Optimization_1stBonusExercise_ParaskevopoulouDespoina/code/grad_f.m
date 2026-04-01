%% function for gradient of f
function gf = grad_f(x,P,q)
gf = P*x + q;
end