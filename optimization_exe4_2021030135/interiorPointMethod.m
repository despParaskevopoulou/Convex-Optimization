function [x, fun_value, rec] = interiorPointMethod(A, b, c, n, x0,x_opt, t1, t2, alpha, beta)
% Initialization 
threshold_1 = t1;
threshold_2 = t2; 
h = 1;             
mu = 10;   

% Linear constraint parameters: a_1^T x -b_1 <= 0 
a_1 = A; 
b_1 = b;

% Parameter of linear Cost Function f(x) = c^T x 
x_init(:,1) = x0;                           % start from a feasible point
t = 1;                                      % initialize parameter t of interior point methods
fun_value = f(c,x_init(:,1));
f_rec = fun_value;
outer_iter = 1;
X(:,outer_iter) = x_init;
k_rec = outer_iter;
x_rec = x0;

l_x = 1000;
while (1)                                     % Outer loop 
     
     x(:,1) = X(:,outer_iter);       
     inner_iter = 1; 
     while ( 1 ) 
         grad = h*gradf(x_init,c,t,n);
         Hess = h*hessf(x_init,n);
         w = -inv(A*inv(Hess)*A')*A*h*inv(Hess)*grad;
         Dx_Nt = -inv(Hess)*(grad+A'*w);          % Newton step
         l_x_ = l_x;
         l_x = sqrt(Dx_Nt'*Hess*Dx_Nt);           % Newton decrement
         if (l_x^2/2 <= threshold_1 || abs(l_x - l_x_)<=10^-7) break; end    % Newton iterations termination condition
         
         tau = 1;
         x_new = x(:,inner_iter) + tau * Dx_Nt;
         while ( not(all(x_new > 0)) )        
             tau = beta * tau; 
             x_new = x(:,inner_iter) + tau * Dx_Nt;
         end 
         % !!!! x_new is FEASIBLE here !!!

         % Backtracking
         while ( t*f(c, x_new) + logBarrier(x_new) > t*f(c, x(:,inner_iter)) + logBarrier(x(:,inner_iter)) + alpha*tau*grad'*Dx_Nt)
              tau = beta * tau; 
              x_new = x(:,inner_iter) + tau * Dx_Nt;
         end

         % Update x
         x(:,inner_iter+1) = x(:,inner_iter) + tau * Dx_Nt;    % Update x
         inner_iter = inner_iter + 1;
     end

     X(:,outer_iter+1) = x(:,inner_iter);           % X: solution of optimization problem      
     if (1/t < threshold_2) 
         break; end              % Algorithm termination condition
     
     
     if(size(x(:,inner_iter),1) == 2)
        plot(x_rec(1,:),x_rec(2,:),'-ro',LineWidth=0.85)
        hold on
        grid on;
     end

     outer_iter = outer_iter + 1;
     t = t * mu;

     fun_value = f(c,x(:,inner_iter));
        
     x_rec = [x_rec, x(:,inner_iter)];
     k_rec = [k_rec, outer_iter];
     f_rec = [f_rec, fun_value];

     
end

fprintf("Total iterations: %d | f = %f\n", outer_iter, fun_value);
rec = [k_rec; f_rec];

x_ = -2:0.1:max(x_opt)+5;
y_ = -(A(1)/A(2))*x_ + b/A(2);
plot(x_, y_, 'c', LineWidth=1)
axis([-1 max(x_opt)+4 -1 max(x_opt)+4])
hold on;

plot(x0(1), x0(2), 'b*', 'MarkerSize', 10, 'LineWidth', 1); 
hold on;
plot(x_opt(1), x_opt(2), 'r*', 'MarkerSize', 10, 'LineWidth', 1);
hold on;
plot(X(1,:), X(2,:))
title('Interior point Method');
grid on;
hold off;
