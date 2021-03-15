clear setup limits guess CONSTANTS

global CONSTANTS TEST

TEST = 1;

CONSTANTS.g = 1.6;

h0 = 10;
hf = 0;
v0 = -2;
vf = 0;

hmin = -20;
hmax =  20;
vmin = -20;
vmax =  20;
umin = -10;
umax =  10;
t0min = 0;
t0max = 0;
tfmin = 0;
tfmax = 1000;

% Phase 1 Information
iphase = 1;
limits(iphase).time.min        = [t0min tfmin];
limits(iphase).time.max        = [t0max tfmax];
limits(iphase).state.min(1,:) = [h0 hmin hf];
limits(iphase).state.max(1,:) = [h0 hmax hf];
limits(iphase).state.min(2,:) = [v0 vmin vf];
limits(iphase).state.max(2,:) = [v0 vmax vf];
limits(iphase).control.min    = 0;
limits(iphase).control.max    = 3;
limits(iphase).parameter.min  = [];
limits(iphase).parameter.max  = [];
limits(iphase).path.min       = [];
limits(iphase).path.max       = [];
limits(iphase).event.min      = [];
limits(iphase).event.max      = [];
limits(iphase).duration.min    = [];
limits(iphase).duration.max    = [];
guess(iphase).time             = [t0min; tfmax];
guess(iphase).state(:,1)      = [h0; h0];
guess(iphase).state(:,2)      = [v0; v0];
guess(iphase).control         = [0; 0];
guess(iphase).parameter       = [];

linkages = [];
setup.name  = 'Moon-Lander-Problem';
setup.funcs.cost = 'moonlanderCost';
setup.funcs.dae = 'moonlanderDae';
setup.limits = limits;
setup.guess = guess;
setup.linkages = linkages;
setup.derivatives = 'finite-difference';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-7;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 10;

[output,gpopsHistory] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;

%--------------------------------%
% END: function moonlanderMain.m %
%--------------------------------%
