figure(1);
pp = plot(solution.time,solution.state(:,1),'-bo',solution.time,solution.state(:,2),'-rd',solution.time,solution.state(:,3),'-gs');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$');
yl = ylabel('$(x_1(t),x_2(t),x_3(t))$');
ll = legend('$x_1(t)$','$x_2(t)$','$x_3(t)$','Location','NorthWest');
set(xl,'FontSize',18,'Interpreter','latex');
set(yl,'FontSize',18,'Interpreter','latex');
set(ll,'FontSize',18,'Interpreter','latex');
grid on;
print -depsc2 brysonDenhamState.eps

figure(2);
pp = plot(solution.time,solution.control,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$');
yl = ylabel('$u(t)$','Interpreter','latex');
set(xl,'FontSize',18,'Interpreter','latex');
set(yl,'FontSize',18,'Interpreter','latex');
grid on;
print -depsc2 brysonDenhamControl.eps

figure(3);
pp = plot(solution.time,solution.costate(:,1),'-bo',solution.time,solution.costate(:,2),'-rd',solution.time,solution.costate(:,3),'-gs');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$');
yl = ylabel('$(\lambda_1(t),\lambda_2(t),\lambda_3(t))$');
ll = legend('$\lambda_1(t)$','$\lambda_2(t)$','$\lambda_3(t)$');
set(xl,'FontSize',18,'Interpreter','latex');
set(yl,'FontSize',18,'Interpreter','latex');
set(ll,'FontSize',18,'Interpreter','latex');
grid on;
print -depsc2 brysonDenhamCostate.eps

