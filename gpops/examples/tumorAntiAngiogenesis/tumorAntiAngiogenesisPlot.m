t = [solution(1).time;solution(2).time];
p = [solution(1).state(:,1);solution(2).state(:,1)];
q = [solution(1).state(:,2);solution(2).state(:,2)];
y = [solution(1).state(:,3);solution(2).state(:,3)];
u = [solution(1).control;solution(2).control];
lambdap = [solution(1).costate(:,1);solution(2).costate(:,1)];
lambdaq = [solution(1).costate(:,2);solution(2).costate(:,2)];

interpt = [solutionPlot(1).time;solutionPlot(2).time];
interpp = [solutionPlot(1).state(:,1);solutionPlot(2).state(:,1)];
interpq = [solutionPlot(1).state(:,2);solutionPlot(2).state(:,2)];
interpy = [solutionPlot(1).state(:,3);solutionPlot(2).state(:,3)];
interplambdap = [solutionPlot(1).costate(:,1);solutionPlot(2).costate(:,1)];
interplambdaq = [solutionPlot(1).costate(:,2);solutionPlot(2).costate(:,2)];

MarkerSize = 8;
LegendSize = 22;
axisSize = 18;

figure(1);
a1 = plot(interpt,interpp,'-r');
hold on
b1 = plot(t,p,'rd','MarkerSize',MarkerSize,'LineWidth',1.25,'MarkerFaceColor','r');
xl = xlabel('t');
yl = ylabel('p(t)');
l1 = legend('Lagrange interpolation','RPM Solution','Location','NorthEast');
set(xl,'FontSize',axisSize);
set(yl,'FontSize',axisSize);
set(l1,'FontSize',LegendSize);
set(a1,'LineWidth',1.25);
grid on;
print -depsc2 pState.eps

figure(2);
a2 = plot(interpt,interpq,'-b');
hold on
b2 = plot(t,q,'bs','MarkerSize',MarkerSize,'LineWidth',1.25,'MarkerFaceColor','b');
x2 = xlabel('t');
y2 = ylabel('q(t)');
l2 = legend('Lagrange interpolation','RPM Solution','Location','SouthEast');
set(x2,'FontSize',axisSize);
set(y2,'FontSize',axisSize);
set(l2,'FontSize',LegendSize);
set(a2,'LineWidth',1.25);
grid on;
print -depsc2 qState.eps
 
figure(3);
a3 = plot(interpt,interpy,'-g');
hold on
b3 = plot(t,y,'go','MarkerSize',MarkerSize,'LineWidth',1.25,'MarkerFaceColor','g');
x3 = xlabel('t');
y3 = ylabel('y(t)');
l3 = legend('Lagrange interpolation','RPM Solution','Location','SouthEast');
set(x3,'FontSize',axisSize);
set(y3,'FontSize',axisSize);
set(l3,'FontSize',LegendSize);
set(a3,'LineWidth',1.25);
grid on;
print -depsc2 yState.eps

figure(4);
b1 = plot(t,u,'mv','MarkerSize',MarkerSize,'MarkerFaceColor','m');
xl = xlabel('t');
yl = ylabel('u(t)');
l1 = legend('RPM Solution','Location','NorthEast');
set(xl,'FontSize',axisSize);
set(yl,'FontSize',axisSize);
set(l1,'FontSize',LegendSize);
grid on;
print -depsc2 control.eps

figure(5);
a1 = plot(interpt,interplambdap,'-r');
hold on
b1 = plot(t,lambdap,'rd','MarkerSize',MarkerSize,'LineWidth',1.25,'MarkerFaceColor','r');
xl = xlabel('t');
yl = ylabel('Costate for p(t)');
l1 = legend('Lagrange interpolation','RPM Solution','Location','NorthEast');
set(xl,'FontSize',axisSize);
set(yl,'FontSize',axisSize);
set(l1,'FontSize',LegendSize);
set(a1,'LineWidth',1.25);
grid on;
print -depsc2 pCostate.eps

figure(6);
a2 = plot(interpt,interplambdaq,'-b');
hold on
b2 = plot(t,lambdaq,'bs','MarkerSize',MarkerSize,'LineWidth',1.25,'MarkerFaceColor','b');
x2 = xlabel('t');
y2 = ylabel('Costate for q(t)');
l2 = legend('Lagrange interpolation','RPM Solution','Location','SouthEast');
set(x2,'FontSize',axisSize);
set(y2,'FontSize',axisSize);
set(l2,'FontSize',LegendSize);
set(a2,'LineWidth',1.25);
grid on;
print -depsc2 qCostate.eps
