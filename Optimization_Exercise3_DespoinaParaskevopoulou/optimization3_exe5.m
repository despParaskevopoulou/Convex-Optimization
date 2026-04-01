%% Fifth Exercise in Optimization
clear all;
close all;

theta = 0:0.01:2*pi;
x1_circle = cos(theta);
x2_circle = sin(theta);
hold on;
plot(x1_circle,x2_circle,'k','LineWidth',1.5,'DisplayName','norm(x) <= 1');

x1 = 0:0.01:1;
x2 = zeros(size(x1));
plot(x1,x2,'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'x1 >= 0');
hold on;
plot(x2,x1,'b-', 'LineWidth', 1.5, ...
    'DisplayName', 'x2 >= 0');
plot(x1, -x1+1, 'g-', 'LineWidth', 1.5, ...
    'DisplayName', 'x1 + x2 <= 1');
axis equal;
axis([-1.5 2 -1.5 2]);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
title('Set $S$', 'Interpreter', 'latex');
legend({'$||x||_2 \leq 1$', '$x_1 \geq 0$' '$x_2 \geq 0$', ...
    '$x_1 + x_2 \leq 1$'}, 'Interpreter', 'latex');
% we will construct the axis in order to ge a better visualization
x = -2:0.01:2;
y = zeros(size(x));
plot(x,y,'k--','LineWidth',1.5);
hold on;
plot(y,x,'k--','LineWidth',1.5);
hold off;