
function events = spacecraftEvent(sol);
global CONSTANTS events;
t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
% rf = xf(:,1:3);
% vf = xf(:,4:6);
% rfdot(:,1) = vf(:,1);
% rfdot(:,2) = vf(:,2);
% rfdot(:,3) = vf(:,3);
load('Phi_max');
load('Rr.mat');

events = Phi_max*[xf(1) xf(2) xf(3) xf(4) xf(5) xf(6)].'-Rr.'


