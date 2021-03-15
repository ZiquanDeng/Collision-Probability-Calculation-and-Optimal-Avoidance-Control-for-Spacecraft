%-----------------------------------------------------%
% Bryson-Denham Example Problem.                      %
% This example is taken from the following reference: %
% Bryson, A. E., Denham, W. F., and Dreyfus, S. E.,   %
% "Optimal Programming Problems with Inequality       %
% Constraints.  I: Necessary Conditions for Extremal  %
% Solutions, AIAA Journal, Vol. 1, No. 11, November,  %
% 1963, pp. 2544-2550.                                %
%-----------------------------------------------------%
clear all
clc

x10 = 0;
x20 = 1;
x30 = 0;
x1f = 0;
x2f = -1;
x1min = 0;
x1max = 1/9;
x2min = -10;
x2max = 10;
x3min = -10;
x3max =  10;

param_min = [];
param_max = [];
path_min = [];
path_max = [];
event_min = [];
event_max = [];
duration_min = [];
duration_max = [];

iphase = 1;
limits(iphase).time.min = [0 0];
limits(iphase).time.max = [0 50];
limits(iphase).state.min(1,:) = [x10 x1min x1f];
limits(iphase).state.max(1,:) = [x10 x1max x1f];
limits(iphase).state.min(2,:) = [x20 x2min x2f];
limits(iphase).state.max(2,:) = [x20 x2max x2f];
limits(iphase).state.min(3,:) = [x30 x3min x3min];
limits(iphase).state.max(3,:) = [x30 x3max x3max];
limits(iphase).control.min    = -5;
limits(iphase).control.max    =  10;
limits(iphase).parameter.min  = param_min;
limits(iphase).parameter.max  = param_max;
limits(iphase).path.min       = path_min;
limits(iphase).path.max       = path_max;
limits(iphase).event.min      = event_min;
limits(iphase).event.max      = event_max;
limits(iphase).duration.min    = [];
limits(iphase).duration.max    = [];
guess(iphase).time             = [0; 1];
guess(iphase).state(:,1)      = [x10; x1f];
guess(iphase).state(:,2)      = [x20; x2f];
guess(iphase).state(:,3)      = [x30; x30];
guess(iphase).control         = [0; 0];
guess(iphase).parameter       = [];

setup.name  = 'Bryson-Denham-Problem';
setup.funcs.cost = 'brysonDenhamCost';
setup.funcs.dae = 'brysonDenhamDae';
setup.limits = limits;
setup.guess = guess;
setup.linkages = [];
setup.derivatives = 'finite-difference';
setup.checkDerivatives = 0;
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 10;
[output,gpopsHistory] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;

%--------------------------------%
% END: script brysonDenhamMain.m %
%--------------------------------%

