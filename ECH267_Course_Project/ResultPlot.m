clc
close all

%% After GPOPS Plot X,Y,Z,dX,dY,dZ,Ux,Uy,Uz
time = solutionPlot.time;
X = solutionPlot.state(:,1);
Y = solutionPlot.state(:,2);
Z = solutionPlot.state(:,3);
dX = solutionPlot.state(:,4);
dY = solutionPlot.state(:,5);
dZ = solutionPlot.state(:,6);
Rs = (X.*X+Y.*Y+Z.*Z).^0.5;
State = solutionPlot.state;
save('aa.mat','State')

Ux = solutionPlot.control(:,1);
Uy = solutionPlot.control(:,2);
Uz = solutionPlot.control(:,3);
Us = (Ux.*Ux+Uy.*Uy+Uz.*Uz).^0.5;

load('No_Strategy.mat');
t_N = 0:1:24;
X_N = XX(1,1:25);
Y_N = XX(2,1:25);
Z_N = XX(3,1:25);
dX_N = XX(4,1:25);
dY_N = XX(5,1:25);
dZ_N = XX(6,1:25);
Rs_N = (X_N.*X_N+Y_N.*Y_N+Z_N.*Z_N).^0.5;

figure
% plot(time,X,'b-',time,Y,'r-',time,Z,'g-',time,Rs,'k-')
% legend('X','Y','Z','Rs')
plot(t_N,X_N,'-k',time,X,'b*')
legend('No Strategy','Take Strategy')
set(gca, 'FontSize', 15)
xlabel('Time/s')
ylabel('State X/km')

figure
plot(t_N,Y_N,'-k',time,Y,'-ro')
legend('No Strategy','Take Strategy')
set(gca, 'FontSize', 15)
xlabel('Time/s')
ylabel('State Y/km')

figure
plot(t_N,Z_N,'-k',time,Z,'-.g+')
legend('No Strategy','Take Strategy')
set(gca, 'FontSize', 15)
xlabel('Time/s')
ylabel('State Z/km')

figure
plot(t_N,Rs_N,'-k',time,Rs,':b*')
legend('No Strategy','Take Strategy')
set(gca, 'FontSize', 15)
xlabel('Time/s')
ylabel('State Rs/km')

figure
plot(time,Ux,'b*',time,Uy,'-ro',time,Uz,'-.g+',time,Us,'-k')
legend('Ux','Uy','Uz','Us')
set(gca, 'FontSize', 15)
xlabel('Time/s')
ylabel('Control Value')
