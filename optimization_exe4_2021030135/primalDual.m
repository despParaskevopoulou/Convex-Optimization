%% Primal Dual Algorithm
function [x, fun_value, rec] = primalDual(A, b, c, n, x0, x_opt, t1, t2, alpha, beta, m, mu)
[p, n] = size(A);

x = x0;

t = 1;

lambda = (1/t)*ones(n,1)./x;
v = A'\(lambda-c);

y = [x; lambda; v];

k = 0;
h = x'*lambda;

x_rec = x0;
f_rec = f(c,x);
k_rec = k;

figure;
while(1)

    t = mu*m/h;

    F = [zeros(n,n), -eye(n), A'; diag(lambda), diag(x), zeros(n,p); A, zeros(p,n), zeros(p,p)];

    delta_y = -F\r(A, b, c, t, y);
    delta_x = delta_y(1:n);
    delta_l = delta_y(n+1:2*n);
    delta_v = delta_y(2*n+1:2*n+p);

    %disp('55');
   
    temp = 1;
    for i = 1:n
        if(-lambda(i)/delta_l(i)<temp && delta_l(i) < 0)
            temp = -lambda(i)/delta_l(i);
        end
    end

    s_max = min(1,temp);

    s = 0.99*s_max;

    x_new = x + s*delta_x;
    
    while(not(all(x_new > -10^-15)))
        %disp('88');
        s = beta*s;
        x_new = x + s*delta_x;
    end

    lambda_new = lambda + s*delta_l;
    v_new = v + s*delta_v;
    y_new = [x_new; lambda_new; v_new];

    while(norm(r(A, b, c, t, y_new)) > (1-alpha*s)*norm(r(A, b, c, t, y)))
        %disp('99');
        s = beta*s;
        x_new = x + s*delta_x;
        lambda_new = lambda + s*delta_l;
        v_new = v + s*delta_v;
        y_new = [x_new; lambda_new; v_new];
    end

    x = x + s*delta_x;
    lambda = lambda + s*delta_l;
    v = v + s*delta_v;

    y = [x; lambda; v];

    k = k + 1;
    fun_value = f(c,x);
    
    x_rec = [x_rec, x];
    k_rec = [k_rec, k];
    f_rec = [f_rec, fun_value];
   

    if(n == 2)
        plot(x_rec(1,:),x_rec(2,:),'-ro',LineWidth=0.85)
        hold on
        grid on;
    end

    h = x'*lambda;
    
    residual = r(A, b, c, t, y);
    if((norm(residual(2*n+1:2*n+p))<=t1) && norm(residual(1:n))<=t1 && (h < t2))
        break
    end

    if(k == 1000) break; end

end

fprintf("Total iterations: %d | f = %f\n", k, fun_value);
rec = [k_rec; f_rec];

x_ = -2:0.1:max(x_opt)+5;
y_ = -(A(1)/A(2))*x_ + b/A(2);
plot(x_, y_, 'c', LineWidth=1)
axis([-1 max(x_opt)+4 -1 max(x_opt)+4])
hold on;

plot(x0(1), x0(2), 'b*', 'MarkerSize', 10, 'LineWidth', 1); 
hold on;
plot(x_opt(1), x_opt(2), 'r*', 'MarkerSize', 10, 'LineWidth', 1);
title('Primal Dual Algorithm');
grid on;
hold off;

end