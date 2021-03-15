clear all
clc
global CONSTANTS events;
t0 = 0;
tf = 24;
% 1.1583111574981046
% 5.100082961050347
% 0.00195991108466842
% -0.004602146117344219
% -0.26371037604012165
% 0.0010610727981537893
% [1.13312431318163,4.03840691576644,0.0130660508694370,-0.00685105257407714,-0.265995177360047,0.00336900176351898]
% [1.10984738111557,2.96413787713957,0.0225860084312159,-0.00546978134357840,-0.269429524078906,0.00206575648361142]

% x0 = 1.10984738111557;
% y0 = 2.96413787713957;
% z0 = 0.0225860084312159;
% 
% r0 = [x0;y0;z0];
% v0 = [-0.00546978134357840; -0.269429524078906; 0.00206575648361142];



x0 = 1.107806576804364;
y0 = 3.142221225079182;
z0 = 8.069775985859451e-04;

r0 = [x0;y0;z0];
v0 = [-0.001158743037575; -0.132814098445661; 0.001165274918660];


% xf=1.10904	;
% yf=1.14889	;
% zf=0.01926;	
% 0.0010	-0.13294	0.00126
% x10dot = -0.0012;
% x20dot = -0.1328;
% x30dot = 0.0012;
% xf = 1.118180642596087;
% yf =-0.047648256169940;
% zf = 0.030534007796103;

% 
% xf=1.08075284143593;
% yf=-0.0460421447605426;
% zf =0.0295095509738497;

xmin = -2;
xmax = 2;
ymin = -2;
ymax = 4;
zmin = -1;
zmax = 1;
vmin = -2;
vmax = -vmin;

param_min = [];
param_max = [];
path_min = 0;
path_max = 1;
event_min = [0,0,0]';
event_max = [0.00,0.00,0.00]';
duration_min = [];
duration_max = [];

iphase = 1;
limits(iphase).time.min = [t0 tf];
limits(iphase).time.max = [t0 tf];
limits(iphase).state.min(1,:) = [r0(1) xmin xmin];
limits(iphase).state.max(1,:) = [r0(1) xmax xmax];
limits(iphase).state.min(2,:) = [r0(2) ymin ymin];
limits(iphase).state.max(2,:) = [r0(2) ymax ymax];
limits(iphase).state.min(3,:) = [r0(3) zmin zmin];
limits(iphase).state.max(3,:) = [r0(3) zmax zmax];
limits(iphase).state.min(4,:) = [v0(1) vmin vmin];
limits(iphase).state.max(4,:) = [v0(1) vmax vmax];
limits(iphase).state.min(5,:) = [v0(2) vmin vmin];
limits(iphase).state.max(5,:) = [v0(2) vmax vmax];
limits(iphase).state.min(6,:) = [v0(3) vmin vmin];
limits(iphase).state.max(6,:) = [v0(3) vmax vmax];
limits(iphase).control.min(1,:)    = 0;
limits(iphase).control.max(1,:)    = 0.25;
limits(iphase).control.min(2,:)    = -0.03;
limits(iphase).control.max(2,:)    = 0;
limits(iphase).control.min(3,:)    = 0;
limits(iphase).control.max(3,:)    = 0.0025;
limits(iphase).parameter.min  = param_min;
limits(iphase).parameter.max  = param_max;
limits(iphase).path.min       = path_min;
limits(iphase).path.max       = path_max;
limits(iphase).event.min      = event_min;
limits(iphase).event.max      = event_max;
limits(iphase).duration.min    = [];
limits(iphase).duration.max    = [];

guess(iphase).time             = [0; tf];
guess(iphase).state(:,1)      = [x0; x0];
guess(iphase).state(:,2)      = [y0; y0];
guess(iphase).state(:,3)      = [z0; z0];
guess(iphase).state(:,4)      = [v0(1); v0(1)];
guess(iphase).state(:,5)      = [v0(2); v0(2)];
guess(iphase).state(:,6)      = [v0(3); v0(3)];
guess(iphase).control(:,1)    = [0.3;0.03];
guess(iphase).control(:,2)    = [0.02;0];
guess(iphase).control(:,3)    = [0.01;0];
guess(iphase).event.min      = event_min;
guess(iphase).event.max      = event_max;
guess(iphase).parameter       = [];

setup.name  = 'Spacecraft-Problem';
setup.funcs.cost = 'spacecraftCost';
setup.funcs.dae = 'spacecraftDae';
setup.funcs.event = 'spacecraftEvent';
setup.limits = limits;
setup.guess = guess;
setup.linkages = [];
setup.derivatives = 'finite-difference';
setup.checkDerivatives = 0;
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-2;
setup.mesh.iteration = 2;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
[output,gpopsHistory] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;

ResultPlot

