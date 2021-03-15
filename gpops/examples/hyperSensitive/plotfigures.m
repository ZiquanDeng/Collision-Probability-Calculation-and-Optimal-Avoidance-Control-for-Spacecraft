figure(1);
pp = plot(solution.time,solution.state,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$x(t)$','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 hypersensitiveState.eps

figure(2);
pp = plot(solution.time,solution.control,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$u(t)$','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 hypersensitiveControl.eps

figure(3);
pp = plot(solution.time,solution.costate,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$\lambda(t)$','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 hypersensitiveCostate.eps
