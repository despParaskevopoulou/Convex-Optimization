%% Sixth Exercise
% first, we generate the P as indicated
A = randn(2,2);
P = A*A';

q = [2; 3];
r = 3;

% we define the range for the coordinates x and y
x1_range = -5:0.1:5;
x2_range = -5:0.1:5;

[x1,x2] = meshgrid(x1_range, x2_range);

f= 0.5 * (x1.^2 * P(1, 1) + x2.^2 * P(2, 2) + ...
    (x1 .* x2) * (P(1, 2) + ...
    P(2, 1))) + q(1) * x1 + q(2) * x2 + r;
P
q
x_optimal = -q\P

figure;
mesh(x1, x2, f);
xlabel('x-axis');
ylabel('y-axis');
zlabel('f(x)');

title('Mesh and Contour Plot of f(x)');
legend('Function Surface', 'Contour Lines', 'Optimal Point');
grid on;

figure;
contour(x1, x2, f, 20);
hold on;
plot(x_optimal(1), x_optimal(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('Contour Plot of f');
legend('Optimal Point x*');