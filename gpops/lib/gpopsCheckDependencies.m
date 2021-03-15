function gpopsCheckDependencies(setup,iphase)
%------------------------------------------------------------------%
% Check sparsity of user defined dependencies                      %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

nstates = setup.sizes(iphase,1);
ncontrols = setup.sizes(iphase,2);
dependenciesDetect = zeros(size(setup.dependencies{iphase}));
%--------------%
% Get user DAE %
%--------------%
sol.time = setup.guess(iphase).time;
sol.state = setup.guess(iphase).state;
if isfield(setup.guess(iphase),'control')
    sol.control = setup.guess(iphase).control;
else
    sol.control = [];
end
if isfield(setup.guess(iphase),'parameter')
    sol.parameter = setup.guess(iphase).parameter;
else
    sol.parameter = [];
end
sol.phase = iphase;
[dae] = feval(setup.funcs.dae,sol);
%------------------------------------------------------%
% Check the Derivative of DAEs with Respect to State   %
%------------------------------------------------------%
for istate=1:nstates
    solPert = sol;
    limitRange = (setup.limits(iphase).state.max(istate,2) - setup.limits(iphase).state.min(istate,2));
    limitRange = max(limitRange,1E-6);
    statePert = 0.01*ones(size(solPert.state(:,istate)))*limitRange;
    limitMean = (setup.limits(iphase).state.max(istate,2) + setup.limits(iphase).state.min(istate,2))/2;
    signPert = sign(solPert.state(:,istate) - limitMean);
    signPert(logical(signPert == 0)) = 1;
    solPert.state(:,istate) = solPert.state(:,istate) + signPert.*statePert;
    [daePert] = feval(setup.funcs.dae,solPert);
    dependenciesDetect(:,istate) = any(daePert-dae,1)';
end;
%------------------------------------------------------%
% Check the Derivative of DAEs with Respect to Control %
%------------------------------------------------------%
for icontrol=1:ncontrols
    solPert = sol;    
    limitRange = (setup.limits(iphase).control.max(icontrol) - setup.limits(iphase).control.min(icontrol));
    limitRange = max(limitRange,1E-6);
    controlPert = 0.01*ones(size(solPert.control(:,icontrol)))*limitRange;
    limitMean = (setup.limits(iphase).control.max(icontrol) + setup.limits(iphase).control.min(icontrol))/2;
    signPert = sign(solPert.control(:,icontrol) - limitMean);
    signPert(logical(signPert == 0)) = 1;
    solPert.control(:,icontrol) = solPert.control(:,icontrol) + signPert.*controlPert;
    [daePert] = feval(setup.funcs.dae,solPert);
    dependenciesDetect(:,icontrol+nstates) = any(daePert-dae,1)';
end;
%------------------------------------------------------%
% find if dependency detected and user set to zero     %
%------------------------------------------------------%
[I,J] = find(dependenciesDetect > setup.dependencies{iphase});
if ~isempty(I)
    errstr = '';
    for i = 1:length(I)
        if I(i) <= nstates
          str = sprintf('State %i derivative depends on',I(i));
        else
          str = sprintf(' Path %i derivative depends on',I(i)-nstates);
        end
        if J(i) <= nstates
          str = sprintf('%s State %i',str,J(i));
        else
          str = sprintf('%s Control %i',str,J(i)-nstates);
        end
        errstr = sprintf('%s\n\t%s, setup.limits(%i).dependencies(%i,%i) = 1',errstr,str,iphase,I(i),J(i));

    end
    error(['DAE dependency detected where user value is zero in phase %i',...
            errstr],iphase)
end
