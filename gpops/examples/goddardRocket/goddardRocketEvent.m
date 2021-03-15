function event = goddardRocketEvent(sol);

global CONSTANTS
t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
p  = sol.parameter;
iphase = sol.phase;

if iphase==2,
    h = xf(1);
    v = xf(2);
    m = xf(3);
    D = CONSTANTS.dragk.*(v.^2).*exp(-h/CONSTANTS.H);
    e0 = m*CONSTANTS.g0-(1+v/CONSTANTS.c).*D;
    event = e0;
else
    event = [];
end;

%---------------------------------%
% End File:  goddardRocketEvent.m %
%---------------------------------%