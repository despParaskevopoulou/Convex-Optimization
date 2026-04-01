%% Third exercise in Optimization (scheme)
x0 = [-1.5;1];
a=[1,1];

projection = max(x0,a);

figure;
quiver(0,0,x0(1),x0(2), 'b--', 'LineWidth', 1.5,'MaxHeadSize',0.2);
hold on;
quiver(0,0,projection(1),projection(2), 'r', 'LineWidth', 1.5,'MaxHeadSize',0.2);
title('Projection onto Set S','Interpreter','latex');
xlabel('x axis','Interpreter','latex');
ylabel('y axis','Interpreter','latex');

legend('$x_0$','Projection', 'Interpreter','latex');