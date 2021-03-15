%---------------------------------------------------------------------%
% BEGIN script tumorAntiAngiogenesisMain.m                            %
%---------------------------------------------------------------------%
% This example is taken from the following reference:                 %
% Ledzewicz, U. and Schattler, H, "Analysis of Optimal Controls for a %
% Mathematical Model of Tumour Anti-angiogenesis," Optimal Control    %
% Applications and Methods, Vol. 29, 2008, pp. 41-57.                 %
%---------------------------------------------------------------------%
clear all
clc

global CONSTANTS

% Parameters:
CONSTANTS.zeta = 0.084;       % per day
CONSTANTS.b = 5.85;           % per day
CONSTANTS.d = 0.00873;        % per mm^2 per day
CONSTANTS.G = 0.15;           % per mg of dose per day
CONSTANTS.mu = 0.02;          % per day
a = 75;           
A = 15;

% Initial Conditions & boundaries:
pMax = ((CONSTANTS.b-CONSTANTS.mu)/CONSTANTS.d)^(3/2);
pMin = 0.1;
qMax = pMax;
qMin = pMin;
yMax = A;
yMin = 0;
uMax = a;
uMin = 0;
t0Max = 0;
t0Min = 0;
tfMax = 5;
tfMin = 0.1;
p0 = pMax/2;
q0 = qMax/4;
y0 = 0;

iphase = 1;
limits(iphase).meshPoints = linspace(-1,1,2);
limits(iphase).nodesPerInterval = [4];
limits(iphase).time.min = [t0Min tfMin];
limits(iphase).time.max = [t0Max tfMax];
limits(iphase).state.min(1,:) = [p0 pMin pMin];
limits(iphase).state.max(1,:) = [p0 pMax pMax];
limits(iphase).state.min(2,:) = [q0 qMin qMin];
limits(iphase).state.max(2,:) = [q0 qMax qMax];
limits(iphase).state.min(3,:) = [y0 yMin yMin];
limits(iphase).state.max(3,:) = [y0 yMax yMax];
limits(iphase).control.min    = uMin;
limits(iphase).control.max    = uMax;
limits(iphase).parameter.min  = [];
limits(iphase).parameter.max  = [];
limits(iphase).path.min       = []; 
limits(iphase).path.max       = []; 
limits(iphase).event.min      = []; 
limits(iphase).event.max      = []; 
limits(iphase).duration.min    = [];
limits(iphase).duration.max    = [];
guess(iphase).time             = [0; 1];
guess(iphase).state(:,1)      = [p0; pMax];
guess(iphase).state(:,2)      = [q0; qMax];
guess(iphase).state(:,3)      = [y0; yMax];
guess(iphase).control         = [uMax; uMax];
guess(iphase).parameter       = [];

setup.name  = 'Tumor-Anti-angiogenesis-Problem';
setup.method = 'radau';
setup.funcs.cost = 'tumorAntiAngiogenesisCost';
setup.funcs.dae = 'tumorAntiAngiogenesisDae';
setup.funcs.event = [];
setup.limits = limits;
setup.guess = guess;
setup.derivatives = 'automatic-intlab';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
output = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;
