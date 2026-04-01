%% Second Exercise
% 1
% define function f(x1,x2)
f12 = @(x1, x2) 1./(1 + x1 + x2);

x1_range = 0:0.2:3;
x2_range = 0:0.2:3;

[x1,x2] = meshgrid(x1_range, x2_range);

f = f12(x1,x2);

figure;
mesh(x1,x2,f);
title('F(x) = 1/1+x1+x2')
xlabel('x')
ylabel('y')
zlabel('z')

grid on;

% 2
contour_levels = linspace(0.1,0.9,10);

figure;
contour(x1,x2,f,contour_levels);
xlabel('x');
ylabel('y');
title('Level sets of f')
grid on;

% 3
% we need the partial derivatives
part_x1 = @(x1, x2) -1/((1+x1+x2).^2);
part_x2 = @(x1, x2) -1/((1+x1+x2).^2);

part2_x1 = @(x1, x2) 2/((x1+x2+1).^3);
part2_x2 = @(x1, x2) 2/((x1+x2+1).^3);
part2_x1_x2 = @(x1, x2) 2/((x1+x2+1).^3);

x01 = 1;
x02 = 1;

x0 = [x01 ; x02];
x = [x1_range ; x2_range];

hf = [2/((x01+x02+1)^3),2/((x01+x02+1)^3) ; 2/((x01+x02+1)^3), 2/((x01+x02+1)^3)];

% first order approximation
f1 = f12(x01,x02) + part_x1(x01,x02) * (x1-x01) + part_x2(x01,x02) * (x2-x02);

% second order approximation
subs = x-x0;
f2 = f12(x01,x02) + ...
    part_x1(x01,x02) * (x1-x01) + part_x2(x01,x02) * (x2-x02) + ...
    (0.5).*(subs')*hf*(x-x0);

% 4
figure;
mesh(x1_range, x2_range, f);
hold on;
mesh(x1_range, x2_range, f1,'EdgeColor','g');
xlabel('x1');
ylabel('x2');
zlabel('f');
title('Mesh Plot of f and Its First-Order Taylor Approximation');
legend('f(x1,x2)', 'First-Order Taylor Approximation');
hold off;

% 5
figure;
mesh(x1_range, x2_range, f);
hold on;
mesh(x1_range, x2_range, f2,'EdgeColor','r');
xlabel('x1');
ylabel('x2');
zlabel('f');
title('Mesh Plot of f and Its Second-Order Taylor Approximation');
legend('f(x1,x2)','Second-Order Taylor Approximation');
hold off;

