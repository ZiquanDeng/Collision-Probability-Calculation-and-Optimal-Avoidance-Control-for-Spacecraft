function relativeError = gpopsSolutionError(setup,iphase)
%------------------------------------------------------------------%
% This function determines an approximation to the solution error  %
% on the current mesh (i.e., the error on the last mesh on which   %
% the problem was solved)                                          %
%------------------------------------------------------------------%
% GPOPS Copyright (c) 2008-2012 Anil V. Rao and David A. Benson    %
%------------------------------------------------------------------%

%------------------------------------------------------------------%
% Interpolate the solution on a mesh that consists of one more LGR %
% point in each mesh interval than was used to solve the NLP       %
%------------------------------------------------------------------%
lagrange = gpopsSolutionInterpolation(setup,iphase);
nstates = setup.sizes(iphase,1);
ncontrols = setup.sizes(iphase,2);
tf = setup.solution(iphase).time(end);
t0 = setup.solution(iphase).time(1);
%-------------------------------------------------------------------%
% Transform the time domain [-1,+1] to [t0,tf] in order to evaluate %
% the differential-algebraic functions on the correct time domain   %
%-------------------------------------------------------------------%
sol.time = (tf-t0)/2*lagrange.T(1:end-1,:)+(tf+t0)/2;
sol.state = lagrange.S(1:end-1,:);
sol.phase = iphase;
sol.control = lagrange.C(1:end,:);
sol.parameter = setup.solution(iphase).parameter;
daeout = feval(setup.funcs.dae,sol);
daeout(:,1:nstates) = (tf-t0)*daeout(:,1:nstates)/2;
meshPoints = setup.limits(iphase).meshPoints;
%----------------------------------------------------------------%
% Set the number of LGR points in each interval to one more than %
% the number of LGR points used to solve the NLP, but keep the   %
% number of mesh intervals the same as was used to solve the NLP % 
%----------------------------------------------------------------%
nodesPerInterval = setup.limits(iphase).nodesPerInterval+1;
numMeshIntervals = size(setup.limits(iphase).nodesPerInterval,2);
%------------------------------------------------------------------%
% Compute the composite integration and unity matrices in order to %
% integrate the interpolated solution at the new set of LGR points %
%------------------------------------------------------------------%
RPM = gpopsRPM(numMeshIntervals,meshPoints,nodesPerInterval);
%-----------------------------------------------------%
% Integrate the dynamics on the interpolated solution %
%-----------------------------------------------------%
oderhs = daeout(:,1:setup.sizes(iphase,1));
integratedRHS = [lagrange.S(1,:); RPM.unityMatrix*lagrange.S+RPM.integrationMatrix*oderhs];
%----------------------------------------------------------------%
% Compute the unscaled difference between the interpolated value %
% of the state and the integrated value of the state.            %
%----------------------------------------------------------------%
unscaledError = abs(integratedRHS-lagrange.S);
stateScalesMatrix = repmat(setup.stateScales{iphase},size(unscaledError,1),1);
%------------------------------------------------------------%
% Scale the unscaled difference between the interpolated and %
% integrated state by unity plus the state scale factors     %
%------------------------------------------------------------%
relativeError = unscaledError./(1+stateScalesMatrix);
