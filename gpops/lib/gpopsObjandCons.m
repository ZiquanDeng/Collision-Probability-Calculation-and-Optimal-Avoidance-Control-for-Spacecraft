function [C] = gpopsObjandCons(x)
%------------------------------------------------------------------%
% Compute the nonlinear constraints and objective function         %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

global mysetup
sumnodes = mysetup.nodes;
Cons = cell(mysetup.numphases,1);
Cost = 0;
for i=1:mysetup.numphases
  nstates = mysetup.sizes(i,1);
  ncontrols = mysetup.sizes(i,2);
  nparameters = mysetup.sizes(i,3);
  npaths = mysetup.sizes(i,4);
  nevents = mysetup.sizes(i,5);
  state_vector = x(mysetup.indices(i).state);
  control_vector = x(mysetup.indices(i).control);
  t0 = x(mysetup.indices(i).time(1));
  tf = x(mysetup.indices(i).time(2));
  tspan = tf-t0;
  t_radau = tspan*(mysetup.ps(i).Points+1)/2+t0;
  state_matrix = reshape(state_vector,sumnodes(i)+1,nstates);    
  state_radau = state_matrix(1:end-1,:);
  x0 = state_matrix(1,:).';
  xf = state_matrix(end,:).';
  control_radau = reshape(control_vector,sumnodes(i),ncontrols);
  if nparameters >0
    parameters = x(mysetup.indices(i).parameter);
  else
    parameters = [];
  end
  %--------------%
  % Get user DAE %
  %--------------%
  sol.time = t_radau;
  sol.state = state_radau;
  sol.control = control_radau;
  sol.parameter = parameters;
  sol.phase = i;
  [dae_out] = feval(mysetup.funcs.dae,sol);
  %-----------------%
  % Get constraints %
  %-----------------%
  odelhs  = mysetup.ps(i).Ddiag*state_matrix;
  % odelhs  = mysetup.ps(i).D*state_matrix;
  oderhs  = tspan*dae_out(:,1:nstates)/2;
  if npaths > 0
    paths = dae_out(:,nstates+1:end);
  else
    paths = [];
  end
  defects = odelhs-oderhs;
  if nevents>0,
    solevents.initial.time = t0;
    solevents.initial.state = x0;
    solevents.terminal.time = tf;
    solevents.terminal.state = xf;
    solevents.parameter = parameters;
    solevents.phase = i;
    [events] = feval(mysetup.funcs.event,solevents);
  else
    events = [];
  end;
  if npaths >0 && nevents > 0
    Cons{i,1} = [defects(:); paths(:); events];
  elseif npaths>0
    Cons{i,1} = [defects(:); paths(:)];
  elseif nevents>0
    Cons{i,1} = [defects(:); events];
  else
    Cons{i,1} = defects(:);
  end
  %---------------%
  % Get user Cost %
  %---------------%
  solcost.initial.time = t0;
  solcost.initial.state = x0;
  solcost.terminal.time = tf;
  solcost.terminal.state = xf;
  solcost.time = t_radau;
  solcost.state = state_radau;
  solcost.control = control_radau;
  solcost.parameter = parameters;
  solcost.phase = i;
  [Mayer,Lagrange] = feval(mysetup.funcs.cost,solcost);
  if ~isempty(Lagrange) 
    integrand = tspan*mysetup.ps(i).Weights*Lagrange/2;
    Cost = Cost + Mayer + integrand;
  else
    Cost = Cost + Mayer;
  end
  solTotal(i) = solcost;
end;
Constraints = vertcat(Cons{:,1});
linkages = mysetup.linkages;
if ~isempty(mysetup.linkages),
  link_out = cell(mysetup.numlinkpairs,1);
  for ipair=1:mysetup.numlinkpairs;
    left_phase  = linkages(ipair).left.phase;
    right_phase = linkages(ipair).right.phase;
    xf_left = solTotal(left_phase).terminal.state;
    p_left  = solTotal(left_phase).parameter;
    x0_right = solTotal(right_phase).initial.state;
    p_right  = solTotal(right_phase).parameter;
    %-------------------------%
    % Get Linkage constraints %
    %-------------------------%
    sollink.left.state = xf_left;
    sollink.left.parameter = p_left;
    sollink.left.phase = left_phase;
    sollink.right.state = x0_right;
    sollink.right.parameter = p_right;
    sollink.right.phase = right_phase;
    [link_out{ipair,1}] = feval(mysetup.funcs.link,sollink);
  end
  Clink = vertcat(link_out{:,1});
  Constraints = [Constraints; Clink];
end
C = [Cost; Constraints; mysetup.initlincons];
