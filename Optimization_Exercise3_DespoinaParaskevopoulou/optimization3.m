%% Third Exercise in Optimization
%% scheme for exercise 1
clear all;
close all;
n = 2;
r = 1;
x0 = [-1.5;1];
figure;
hold on;

% the ball
theta = linspace(0, 2*pi, 100);
ball_x = r * cos(theta);
ball_y = r * sin(theta);
plot(ball_x, ball_y, 'b');

% the point x0
quiver(0,0,x0(1), x0(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.2);

% the line from origin to x0
%plot([0, x0(1)], [0, x0(2)], 'k--'); 

% the projection
projection = r * x0 / norm(x0);
quiver(0,0,projection(1), projection(2), 'g-','LineWidth',2,'MaxHeadSize',0.2); 

axis equal;
xlabel('$x-axis$', 'Interpreter', 'latex')
ylabel('$y-axis$', 'Interpreter', 'latex')
legend('Ball B(0,r)','$x_0$','Line From Origin to $x_0$', 'Interpreter', 'latex')
title('Scheme for the first problem','Interpreter','latex');
hold off;

%% scheme for exercise 2
clear all;
close all;
x0 = [-1.5; 1];
r = 1;
y = [1; 1];
x = r/(norm(x0-y))*(x0-y)+y;

%  the projection onto B(y, r)
v = x0 - y;
projection = y + r * v / norm(v);

figure
quiver(0, 0, x0(1), x0(2), 'r', 'LineWidth', 1, 'MaxHeadSize', 0.2)
hold on 
quiver(0,0,projection(1), projection(2), 'g--','LineWidth',2,'MaxHeadSize',0.2);
hold on
quiver(0, 0, x(1), x(2), 'k-', 'LineWidth', 1, 'MaxHeadSize', 0.2)
theta = linspace(0, 2*pi, 100);
x_circle = y(1) + r*cos(theta);
y_circle = y(2) + r*sin(theta);
plot(x_circle, y_circle, 'b', 'LineWidth', 1.2)
axis([-2 2.5 -.5 2.5])
axis equal
title('Scheme for the second problem','Interpreter','latex');
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
legend('$\mathbf{x}_0$','$\mathbf{x}$', 'Interpreter', 'latex')
