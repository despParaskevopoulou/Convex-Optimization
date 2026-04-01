%% First exercise in Optimization 
%% 1
syms x

% Define x values
x = linspace(0.1, 5, 1000);
f =  @(x) 1./(1 + x);

% 1st and 2nd derivatives of f(x)
f_der1 = @(x) -1./((1 + x).^2);
f_der2 = @(x) 2./((1 + x).^3);

syms x0

% Define x0 values
x0_values1 = [2:5];

% we define different styles for lines
% to make the diagrams more clear
line_styles = {'-', '--', '-.'};

figure;
hold on;

for i = 1:length(x0_values1)
    x0 = x0_values1(i);
    f1 = f(x0) + f_der1(x0)*(x-x0);
    f2 = f(x0) + f_der1(x0)*(x-x0) + 0.5 * f_der2(x0) * (x-x0).^2; 
    subplot(2,2,i);
    plot(x,f(x), line_styles{1}, 'LineWidth', 1, 'DisplayName',sprintf('f(x) for x_0=%d', x0));
    hold on;
    plot(x,f1, line_styles{2}, 'LineWidth', 1, 'DisplayName',sprintf('f1(x) for x_0=%d', x0));
    plot(x,f2, line_styles{3}, 'LineWidth', 1, 'DisplayName',sprintf('f2(x) for x_0=%d', x0));
    title(['x_0 = ' num2str(x0)]);
    xlabel('x');
    ylabel('y');
    legend('show');
    grid on;
    hold off;
end



