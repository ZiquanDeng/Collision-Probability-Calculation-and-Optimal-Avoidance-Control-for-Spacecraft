figure(1);
pp = plot(solutionPlot.time(1:5:end),solutionPlot.state(1:5:end,1)/1000,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$h(t)$ (km)','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 brysonMinClimbAltitude.eps

figure(2);
pp = plot(solutionPlot.time(1:5:end),solutionPlot.state(1:5:end,2),'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$v(t)$ (m/s)','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 brysonMinClimbSpeed.eps

figure(3);
pp = plot(solutionPlot.time(1:5:end),solutionPlot.state(1:5:end,3)*180/pi,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$\gamma(t)$ (deg)','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 brysonMinClimbFlightPathAngle.eps

figure(4);
pp = plot(solutionPlot.time(1:5:end),solutionPlot.state(1:5:end,4)/1000,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$m(t)$ (kg$\times 1000$)','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 brysonMinClimbMass.eps

figure(5);
pp = plot(solutionPlot.time(1:5:end),solutionPlot.control(1:5:end,:)*180/pi,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$\alpha(t)$ (deg)','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
print -depsc2 brysonMinClimbAttackAngle.eps

