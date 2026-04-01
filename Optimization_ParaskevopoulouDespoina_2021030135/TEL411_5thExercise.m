%% 5th exercise
%% c
% first we plot for a >= 1 and a <= 0 
x = linspace(0.1, 2, 5);
a_values_i = [1, -0.5, -1, 3, 0];
figure;
hold on;
for a = a_values_i
    y = x.^a;
    plot(x, y, 'LineWidth', 2, 'DisplayName', ['a = ', num2str(a)]);
end
hold off;
xlabel('x');
ylabel('f(x)');
title('f(x) = x^a for a \geq 1 and a \leq 0');
legend;
grid on;

% now for a>=0 and a <=1
a_values_ii = [0.5 , 0.8, 1];
figure;
hold on;
for a = a_values_ii
    y = x.^a;
    plot(x, y, 'LineWidth', 2, 'DisplayName', ['a = ', num2str(a)]);
end
hold off;
xlabel('x');
ylabel('f(x)');
title('f(x) = x^a for 0 \leq a \leq 1');
legend;
grid on


%% d
x1_range = -2:0.1:2;
x2_range = -2:0.1:2;

[x1,x2] = meshgrid(x1_range, x2_range);
f1 = sqrt(x1.^2 + x2.^2);
f2 = f1.^2;

figure;
mesh(x1,x2,f1)
title('Euclidean Norm f1');

figure
mesh(x1,x2,f2)
title('Square of Euclidean Norm')

xlabel('x1')
ylabel('x2')
