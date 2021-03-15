function events = orbitRaisingMaxRadiusEvent(solevent);

global CONSTANTS

t0      = solevent.initial.time;
rf      = solevent.terminal.state(1);
vthetaf = solevent.terminal.state(4);
p       = solevent.parameter;

events(1,:) = sqrt(CONSTANTS.mu/rf)-vthetaf;
events(2,:) = p-t0;

