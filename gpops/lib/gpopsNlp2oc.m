function setup = gpopsNlp2oc(setup)
%------------------------------------------------------------------%
% Convert the NLP solution to trajectory format                    %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David A. Benson              %
%------------------------------------------------------------------%

clear snoptcmex
funcs = setup.funcs;
ps = setup.ps;
result = setup.result;
variable_indices = setup.variable_indices;
constraint_indices = setup.constraint_indices;
sizes = setup.sizes;
nodes = setup.nodes;
numnonlin = setup.numnonlin;
numphases = setup.numphases;
numlinks = setup.numlinks;

x = result.x;
lambda  = result.Fmul(1:numnonlin-numlinks);
cost = 0;
varshift = 0;
conshift = 0;
solution = repmat(struct('time',[],'state',[],'control',[],...
                         'parameter',[],'costate',[],'pathmult',[],'Hamiltonian',[],...
                         'Mayer_cost',[],'Lagrange_cost',[]),1,numphases);
for i=1:numphases;
  %------------------------------------------------------------%
  % First construct the primal solution to the optimal control %
  % problem (i.e., states, controls, and parameters)           %
  %------------------------------------------------------------%
  xcurr = x(variable_indices{i});
  nstates = sizes(i,1);
  ncontrols = sizes(i,2);
  nparameters = sizes(i,3);
  npaths = sizes(i,4);
  state_indices = 1:(nodes(i)+1)*nstates;
  if ncontrols>0,
    control_indices = state_indices(end)+1:state_indices(end)+nodes(i)*ncontrols;
    t0_index = control_indices(end)+1;
  else
    control_indices = [];
    t0_index = state_indices(end)+1;
  end;
  tf_index = t0_index+1;   
  parameter_indices = tf_index+1:tf_index+nparameters;
  t0 = xcurr(t0_index);
  tf = xcurr(tf_index);
  tspan = tf - t0;
  state_vector = xcurr(state_indices);
  control_vector = xcurr(control_indices);
  tau_all = [ps(i).Points; 1];
  t_all = (tf-t0)*(tau_all+1)/2+t0;
  t_gauss = t_all(1:end-1);
  if nstates>0,
    state_matrix = reshape(state_vector,nodes(i)+1,nstates);
  else
    state_matrix = [];
  end;
  if ncontrols>0,
    control_matrix = reshape(control_vector,nodes(i),ncontrols);
    control_t0 = zeros(1,ncontrols);
    control_tf = control_t0;
    try,
      for j=1:ncontrols
        control_tf(j) = interp1(tau_all(1:end-1),control_matrix(:,j),+1,'spline','extrap');
      end;
    catch
      keyboard
    end;
    control_matrixTotal = [control_matrix; control_tf];
    control_matrixTotalExtrap = [control_matrix; control_tf];
  else
    control_matrixTotal = [];
    control_matrixTotalExtrap = [];
  end;
  parameter = xcurr(parameter_indices);
  %------------------------------%
  % Next, construct the costates %
  %------------------------------%
  lamcurr = lambda(conshift+1:conshift+(nodes(i))*nstates);
  lamcurr = reshape(lamcurr,nodes(i),nstates);
  Inverse_Weight_Matrix = diag(1./ps(i).Weights);
  costate_gauss = Inverse_Weight_Matrix*lamcurr(1:end,:);
  DD = ps(i).D(:,end);
  costatef = DD.'*lamcurr(1:end,:);
  costate = [costate_gauss; costatef];
  if npaths>0,
    pathmult = lambda(conshift+(nodes(i))*nstates+1:conshift+(nodes(i))*nstates+npaths*nodes(i));
    pathmult = reshape(pathmult,nodes(i),npaths);
    for j=1:npaths
      pathmult(:,j) = 2*Inverse_Weight_Matrix*pathmult(:,j)/(tf-t0);
    end;
    if ncontrols>0,
      pathmultF = zeros(1,npaths);
      for k=1:npaths
        pathmultF(k) = interp1(tau_all(1:end-1),pathmult(:,k),1,'spline','extrap');
        % pathmultF(k) = interp1(t_all(1:end-1),pathmult(:,k),t_all(end),'spline','extrap');
      end;
      pathmultTotal = [pathmult; pathmultF];
    else
      pathmultF = zeros(1,npaths);
      for k=1:npaths
        pathmultF(k) = interp1(tau_all(1:end-1),pathmult(:,k),1,'spline','extrap');
        % pathmultF(k) = interp1(t_all(1:end-1),pathmult(:,k),t_all(end),'spline','extrap');
      end;
      pathmultTotal = [pathmult; pathmultF];
    end
  else
    pathmultTotal = [];
  end;
  %-------------------------------------------------------------------------%
  % Put all of the information from the primal and dual solution            %
  %  into the appropriate location in the cell array SOLUTION               %
  %  i = phase number                                                       %
  %  solution(i).time: vector containing time                               %
  %  solution(i).state: N+1 by nx matrix containing states                  %
  %  solution(i).control: N+1 by nu matrix containing controls              %
  %  solution(i).parameter: vector of length q containing static parameters %
  %  solution(i).costate: N+1 by np matrix containing costates              %
  %  solution(i).pathmult: N+1 by np matrix containing path multipliers     %
  %  solution(i).Hamiltonian: vector of length N+1 containing Hamiltonian   %
  %  solution(i).Mayer_cost: scalar containing the Mayer cost               %
  %  solution(i).Lagrange_cost: scalar containing the Lagrange cost         %
  %-------------------------------------------------------------------------%
  sol.initial.time = t_all(1);
  sol.initial.state = state_matrix(1,:).';
  sol.terminal.time = t_all(end);
  sol.terminal.state = state_matrix(end,:).';
  sol.time = t_all;
  sol.state = state_matrix;
  sol.control = control_matrixTotalExtrap;
  sol.parameter = parameter;
  sol.phase = i;
  [Mayer,Lagrange] = feval(setup.funcs.cost,sol);
  soldae.time = t_all;
  soldae.state = state_matrix;
  soldae.control = control_matrixTotal;
  soldae.parameter = parameter;
  soldae.phase = i;
  dae = feval(setup.funcs.dae,soldae);
  Hamiltonian = Lagrange+sum(costate.*dae(:,1:nstates),2);
  varshift = varshift+length(variable_indices{1,i});
  conshift = conshift+length(constraint_indices{1,i});
  Lagrange_Cost = (tf-t0)*setup.ps(i).Weights*Lagrange(1:end-1)/2;
  cost = cost + Mayer + Lagrange_Cost;
  solution(i).time = soldae.time;
  solution(i).state = soldae.state;
  solution(i).control = control_matrixTotalExtrap;
  solution(i).parameter = soldae.parameter;
  solution(i).costate = costate;
  solution(i).pathmult = pathmultTotal;
  solution(i).Hamiltonian = Hamiltonian;
  solution(i).Mayer_cost = Mayer;
  solution(i).Lagrange_cost = Lagrange_Cost;
end;
setup.solution = solution;
setup.guess = solution;
setup.cost     = cost;
