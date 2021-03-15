function dae = gpopsDaeWrapper(t,xup,extras);
%------------------------------------------------------------------%
% Wrapper function for differential-algebraic equations            %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

setup = extras.setup;
iphase = extras.phase;
nstates = extras.nstates;
ncontrols = extras.ncontrols;
nparameters = extras.nparameters;
sol.time = t;
sol.state = xup(1:nstates).';
sol.control = xup(nstates+1:nstates+ncontrols).';
sol.parameter = xup(nstates+ncontrols+1:nstates+ncontrols+nparameters);
sol.phase = iphase;
dae = feval(setup.funcs.dae,sol).';

