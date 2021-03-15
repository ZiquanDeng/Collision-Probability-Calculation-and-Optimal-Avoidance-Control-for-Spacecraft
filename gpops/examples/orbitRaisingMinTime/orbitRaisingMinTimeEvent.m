function event = orbitRaisingMinTimeEvent(solevent);

t0 = solevent.initial.time;
p  = solevent.parameter;

event(1,:) = t0-p;
