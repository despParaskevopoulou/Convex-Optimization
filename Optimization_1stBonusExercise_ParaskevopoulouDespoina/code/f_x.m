%% function for computing f(x)
function f = f_x(x,P,q)
f = (1/2)*x'*P*x + (q')*x; 
end
