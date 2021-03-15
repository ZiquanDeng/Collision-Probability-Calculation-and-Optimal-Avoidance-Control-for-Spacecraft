function event = hohmannTransferEvent(sol);

global CONSTANTS
t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
p  = sol.parameter;
iphase = sol.phase;
v0_plus = x0(4:6);
dv1 = p(1:3);
dv2 = p(4:6);
vf_plus = xf(4:6) + dv2;
radf = sqrt(xf(1:3).'*xf(1:3));
speedf = sqrt(vf_plus.'*vf_plus);
singf = xf(1:3).'*vf_plus/(radf*speedf);

oe = hohmannTransferRv2oe(xf(1:3),vf_plus,CONSTANTS.mu);
af = oe(1);
ef = oe(2);
nuf = oe(6);
event0 = CONSTANTS.v0 + dv1 -v0_plus;
eventf = [radf; speedf; singf];
event = [event0; eventf];
