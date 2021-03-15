%-------------------------------------------------------%
% Infinte Horizon Optimal Control Problem               %
%-------------------------------------------------------%
% This example is a modification of the problem found   %
% in the following reference:                           %
%-------------------------------------------------------%
% Garg, D., Patterson, M. A., Hager, W. W., Rao, A. V., %
% "Pseudospectral Methods for Solving Infinite-Horizon  %
% Optimal Control Problems," Automatica, Vol. 47, No. 4,%
% April 2011, pp. 829-837.                              %
%-------------------------------------------------------%
% The problem statement is given as follows.  Minimize  % 
%      0.5 * int(ln(y)^2 + u^2)dt                       %
% subject to the dynamic constraint                     %
%      dy/dt = y * ln(y) + y * u                        %
% with the initial condition                            %
%      y(0) = 2                                         %
%                                                       %
%-------------------------------------------------------%
clear all
clc

y0 = 2;

iphase = 1;

limits(iphase).meshPoints = [-1 +1];
limits(iphase).nodesPerInterval = [10];
limits(iphase).time.min = [-1 1];
limits(iphase).time.max = [-1 1];
limits(iphase).state.min(1,:) = [y0 1 1];
limits(iphase).state.max(1,:) = [y0 10 10];
limits(iphase).control.min    = -100;
limits(iphase).control.max    =  100;
limits(iphase).parameter.min  = [];
limits(iphase).parameter.max  = [];
limits(iphase).path.min       = [];
limits(iphase).path.max       = [];
limits(iphase).event.min      = [];
limits(iphase).event.max      = [];
limits(iphase).duration.min   = [];
limits(iphase).duration.max   = [];

guess(iphase).time            = [-1; 1];
guess(iphase).state           = [y0; 1];
guess(iphase).control         = [0; 0];
guess(iphase).parameter       = [];

setup.name  = 'Infinite-Horizon-Problem';
setup.funcs.cost = 'infHorizonCost';
setup.funcs.dae = 'infHorizonDae';
setup.limits = limits;
setup.guess = guess;
setup.linkages = [];
setup.derivatives = 'finite-difference';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;

[output] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;

tau = solution.time;
t = (1+tau)./(1-tau);

% true solution
x = log(2)*exp(-t*sqrt(2));
y = exp(x);
u = -(1+sqrt(2))*x;
lam = (1+sqrt(2))*exp(-x).*x;

figure(1);
pp = plot(tau,solution.state,'o-',tau,y,'-s');
set(pp,'LineWidth',1.5);
xl = xlabel('$\tau$','Interpreter','LaTeX');
yl = ylabel('$y(\tau)$','Interpreter','LaTeX');
ll = legend('Pseudospectral','Exact','Location','Best');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(ll,'FontSize',16,'FontName','Times');
set(gca,'FontName','Times','FontSize',14);
grid on

figure(2);
pp = plot(tau,solution.control,'o-',tau,u,'-s');
set(pp,'LineWidth',1.5);
xl = xlabel('$\tau$','Interpreter','LaTeX');
yl = ylabel('$u(\tau)$','Interpreter','LaTeX');
ll = legend('Pseudospectral','Exact','Location','Best');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(ll,'FontSize',16,'FontName','Times');
set(gca,'FontName','Times','FontSize',14);
grid on

figure(3);
pp = plot(tau,solution.costate,'o-',tau,lam,'-s');
set(pp,'LineWidth',1.5);
xl = xlabel('$\tau$','Interpreter','LaTeX');
yl = ylabel('$\lambda(\tau)$','Interpreter','LaTeX');
ll = legend('Pseudospectral','Exact','Location','Best');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(ll,'FontSize',16,'FontName','Times');
set(gca,'FontName','Times','FontSize',14);
grid on

