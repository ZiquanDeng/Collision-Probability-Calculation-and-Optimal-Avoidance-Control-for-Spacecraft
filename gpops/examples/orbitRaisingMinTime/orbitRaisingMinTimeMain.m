global CONSTANTS

r0 = 1;
rf = 1.5;
u0 = 0;
uf = 0;
v0 = 1;
vf = 0.816496581;
m0 = 1;
mf = 0;
CONSTANTS.m0 = m0;
CONSTANTS.mdot = 0.0749;
CONSTANTS.mu  = 1; 
CONSTANTS.T   = 0.1405;
CONSTANTS.Ve  = 1.8758;

rmin = 1;
rmax = 2;
umin = -0.1;
umax = 0.5;
vmin = 0.5;
vmax = 1.1;
mmax = 1.05;
mmin = 0.5;

uumin =  -4*pi;
uumax =   4*pi;
t0min = 0;
t0max = 0;
tfmin = 0.5;
tfmax = 10;

% Phase 1 Information
iphase = 1;
limits(iphase).intervals = 10;
limits(iphase).nodesperint = 3;
limits(iphase).time.min = [t0min tfmin];
limits(iphase).time.max = [t0max tfmax];
limits(iphase).state.min(1,:) = [r0 rmin rf];
limits(iphase).state.max(1,:) = [r0 rmax rf];
limits(iphase).state.min(2,:) = [u0 umin uf];
limits(iphase).state.max(2,:) = [u0 umax uf];
limits(iphase).state.min(3,:) = [v0 vmin vf];
limits(iphase).state.max(3,:) = [v0 vmax vf];
limits(iphase).state.min(4,:) = [m0 mmin mmin];
limits(iphase).state.max(4,:) = [m0 mmax mmax];
limits(iphase).control.min(1,:) = -1;
limits(iphase).control.max(1,:) =  1;
limits(iphase).control.min(2,:) = -1;
limits(iphase).control.max(2,:) =  1;
limits(iphase).parameter.min = [];
limits(iphase).parameter.max = [];
limits(iphase).path.min = 1;
limits(iphase).path.max = 1;
limits(iphase).event.min = [];
limits(iphase).event.max = [];
guess(iphase).time = [t0min; 5];
guess(iphase).state(:,1) = [r0; rf];
guess(iphase).state(:,2) = [u0; uf];
guess(iphase).state(:,3) = [v0; vf];
guess(iphase).state(:,4) = [m0; m0];
guess(iphase).control(:,1) = [1; 1];
guess(iphase).control(:,2) = [0; 0];
guess(iphase).parameter = [];

linkages = [];
setup.name  = 'Orbit-Raising-Min-Time-Problem';
setup.funcs.cost = 'orbitRaisingMinTimeCost';
setup.funcs.dae = 'orbitRaisingMinTimeDae';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-7;
setup.derivatives = 'finite-difference';
setup.limits = limits;
setup.guess = guess;
setup.linkages = linkages;
[output,gpopsHistory] = gpops(setup);
solution = output.solution;
