%------------------------------------------------------------------%
% BEGIN: script hyperSensitiveMain.m                               %
%------------------------------------------------------------------%
% This example is taken from the following reference:              %
% Rao, A. V. and Mease, K. D., "Eigenvector Approximate Dichotomic %
% Basis Method for Solving Hypersensitive Optimal Control          %
% Problems," Optimal Control Applications and Methods, Vol. 21,    %
% No. 1, 2000, pp. 1-19.                                           %
% The optimal control problem is described as follows:             %
% Minimize J = 0.5*(x^2+u^2)                                       %
% subject to the dynamic constraint                                % 
%   xdot = -x^3 + u                                                %
% and the boundary conditions                                      %
%   x(0) = 1.5                                                     %
%   x(tf) = 1                                                      %
%------------------------------------------------------------------%
clear all
clc

t0 = 0; tf = 5000; x0 = 1.5; xf = 1;
xmin = -10; xmax =  10; umin = -10; umax =  10;

iphase = 1;
limits(iphase).time.min = [t0 tf];
limits(iphase).time.max = [t0 tf];
limits(iphase).state.min = [x0 xmin xf];
limits(iphase).state.max = [x0 xmax xf];
limits(iphase).control.min = umin;
limits(iphase).control.max = umax;
limits(iphase).parameter.min = [];
limits(iphase).parameter.max = [];
limits(iphase).path.min = [];
limits(iphase).path.max = [];
limits(iphase).event.min = [];
limits(iphase).event.max = [];
guess(iphase).time = [t0; tf];
guess(iphase).state(:,1) = [x0; xf];
guess(iphase).control = [-1; 1];
guess(iphase).parameter = [];

setup.name  = 'HyperSensitive-Problem';
setup.funcs.cost = 'hyperSensitiveCost';
setup.funcs.dae = 'hyperSensitiveDae';
setup.linkages = [];
setup.limits = limits;
setup.guess = guess;
setup.derivatives = 'finite-difference';
setup.checkDerivatives = 0;
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 10;

[output,gpopsHistory] = gpops(setup);
%----------------------------------%
% END: script hyperSensitiveMain.m %
%----------------------------------%
