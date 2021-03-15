function [C J] = gpopsObjandConsFiniteDifference(x)
%------------------------------------------------------------------%
% Compute the nonlinear constraints and objective function using   %
% the complex-step derivative approximation.                       %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

global mysetup

[ncon, nvar] = size(mysetup.sparsity_all);
J = spalloc(ncon,nvar,length(mysetup.iGfun));
pertFactor = mysetup.tolerances(1);
sumnodes = mysetup.nodes;
Cons = cell(mysetup.numphases,1);
Cost = 0;
for i=1:mysetup.numphases
  nnodes = mysetup.limits(i).nodes;
  nstates = mysetup.sizes(i,1);
  ncontrols = mysetup.sizes(i,2);
  nparameters = mysetup.sizes(i,3);
  npaths = mysetup.sizes(i,4);
  nevents = mysetup.sizes(i,5);
  state_vector = x(mysetup.indices(i).state);
  control_vector = x(mysetup.indices(i).control);
  t0 = x(mysetup.indices(i).time(1));
  tf = x(mysetup.indices(i).time(2));
  if isequal(mysetup.autoscale,'on'),
    stateScales = repmat(mysetup.stateScales{i},nnodes,1);
    controlScales = repmat(mysetup.controlScales{i},nnodes,1);
    parameterScales = mysetup.parameterScales{i};
    x0Scales = mysetup.stateScales{i}.';
    xfScales = mysetup.stateScales{i}.';
    t0Scales = mysetup.t0Scales{i};
    tfScales = mysetup.tfScales{i};
    timeScales = repmat((t0Scales+tfScales)/2,nnodes,1);
  else
    stateScales = ones(nnodes,nstates);
    controlScales = ones(nnodes,ncontrols);
    parameterScales = ones(nparameters,1);
    x0Scales = ones(nstates,1);
    xfScales = ones(nstates,1);
    t0Scales = 1;
    tfScales = 1;
    timeScales = ones(nnodes,1);
  end;
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
  %-------------------------------------------------------%
  % Compute Finite-Difference Derivatives of DAE Function %
  %-------------------------------------------------------%
  pertt0 = pertFactor*((t0Scales+tfScales)/2+1);
  perttf = pertFactor*((t0Scales+tfScales)/2+1);
  pertx0 = pertFactor*(x0Scales+1);
  pertxf = pertFactor*(xfScales+1);
  pertTime = pertFactor*(timeScales+1);
  pertState = pertFactor*(stateScales+1);
  pertControl = pertFactor*(controlScales+1);
  if nparameters>0,
    pertParameter = pertFactor*(parameterScales+1);
  else
    pertParameter = [];
  end;
  tRadauPert = t_radau + pertTime;
  stateRadauPert = state_radau + pertState;
  controlRadauPert = control_radau + pertControl;
  parameterPert = parameters + pertParameter;
  DDae.time  = zeros(nnodes,nstates+npaths);
  DDae.state = cell(1,nstates);
  DDae.control = cell(1,ncontrols);
  DDae.parameter = cell(1,nparameters);
  %------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to Time. %
  %------------------------------------------------------%
  sol.time = tRadauPert;
  denominator = repmat(pertTime,1,nstates+npaths);
  DDae.time = (feval(mysetup.funcs.dae,sol)-dae_out)./denominator;
  sol.time = t_radau;
  %------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to State %
  %------------------------------------------------------%
  for istate=1:nstates
    sol.state(:,istate) = stateRadauPert(:,istate);
    denominator = repmat(pertState(:,istate),1,nstates+npaths);
    DDae.state{istate}=(feval(mysetup.funcs.dae,sol)-dae_out)./denominator;
    sol.state(:,istate) = state_radau(:,istate);
  end;
  %---------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to Control. %
  %---------------------------------------------------------%
  for icontrol=1:ncontrols
    sol.control(:,icontrol) = controlRadauPert(:,icontrol);
    denominator = repmat(pertControl(:,icontrol),1,nstates+npaths);
    DDae.control{icontrol} = (feval(mysetup.funcs.dae,sol)-dae_out)./denominator;
    sol.control(:,icontrol) = control_radau(:,icontrol);
  end;
  %-----------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to Parameter. %
  %-----------------------------------------------------------%
  for iparameter=1:nparameters
    sol.parameter(iparameter) = parameterPert(iparameter);
    denominator = repmat(pertParameter(iparameter),nnodes,nstates+npaths);
    DDae.parameter{iparameter} = (feval(mysetup.funcs.dae,sol)-dae_out)./denominator;
    sol.parameter(iparameter) = parameters(iparameter);
  end;
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
  %-----------------------------------%
  % Get Derivatives of Event Function %
  %-----------------------------------%
  if nevents>0,
    DEvent.t0 = zeros(nevents,1);
    DEvent.x0 = cell(1,nstates);
    DEvent.tf = zeros(nevents,1);
    DEvent.xf = cell(1,nstates);
    DEvent.parameter = cell(1,nparameters);
    solevents.initial.time = t0;
    solevents.initial.state = x0;
    solevents.terminal.time = tf;
    solevents.terminal.state = xf;
    solevents.parameter = parameters;
    solevents.phase = i;
    [events] = feval(mysetup.funcs.event,solevents);
    t0Pert = t0 + pertt0;
    x0Pert = x0 + pertx0;
    tfPert = tf + perttf;
    xfPert = xf + pertxf;
    %-------------------------------------------------------%
    % Get Derivatives of Event with Respect to Initial Time %
    %-------------------------------------------------------%
    solevents.initial.time = t0Pert;
    DEvent.t0 = (feval(mysetup.funcs.event,solevents)-events)./pertt0;
    solevents.initial.time = t0;
    %--------------------------------------------------------%
    % Get Derivatives of Event with Respect to Terminal Time %
    %--------------------------------------------------------%
    solevents.terminal.time = tfPert;
    DEvent.tf = (feval(mysetup.funcs.event,solevents)-events)./perttf;
    solevents.terminal.time = tf;
    %---------------------------------------------------------------------%
    % Get Derivatives of Event with Respect to Initial and Terminal State %
    %---------------------------------------------------------------------%
    for istate=1:nstates;
      solevents.initial.state(istate) = x0Pert(istate);
      DEvent.x0{istate} = (feval(mysetup.funcs.event,solevents)-events)./pertx0(istate);
      solevents.initial.state(istate) = x0(istate);
      solevents.terminal.state(istate) = xfPert(istate);
      DEvent.xf{istate} = (feval(mysetup.funcs.event,solevents)-events)./pertxf(istate);
      solevents.terminal.state(istate) = xf(istate);
    end
    %-----------------------------------------------------%
    % Get Derivatives of Event with Respect to Parameters %
    %-----------------------------------------------------%
    for iparameter=1:nparameters;
      solevents.parameter(iparameter) = parameterPert(iparameter);
      DEvent.parameter{iparameter} = (feval(mysetup.funcs.event,solevents)-events)./pertParameter(iparameter);
      solevents.parameter(iparameter) = parameters(iparameter);
    end;
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
  DMayer.t0 = zeros(1,1);
  DMayer.x0 = zeros(1,nstates);
  DMayer.tf = zeros(1,1);
  DMayer.xf = zeros(1,nstates);
  DMayer.parameter = zeros(1,nparameters);
  DLagrange.time = zeros(nnodes,1);
  DLagrange.state = zeros(nnodes,nstates);
  DLagrange.control = zeros(nnodes,ncontrols);
  DLagrange.parameter = zeros(nnodes,nparameters);
  [Mayer,Lagrange] = feval(mysetup.funcs.cost,solcost);
  t0Pert = t0 + pertt0;
  x0Pert = x0 + pertx0;
  tfPert = tf + perttf;
  xfPert = xf + pertxf;
  if ~isempty(Mayer),
    %------------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Initial Time %
    %------------------------------------------------------------%
    solcost.initial.time = t0Pert;
    [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
    DMayer.t0 = (ComplexMayer-Mayer)./pertt0;
    solcost.initial.time = t0;
    %-------------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Terminal Time %
    %-------------------------------------------------------------%
    solcost.terminal.time = tfPert;
    [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
    DMayer.tf = (ComplexMayer-Mayer)./perttf;
    solcost.terminal.time = tf;
    %--------------------------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Initial and Terminal State %
    %--------------------------------------------------------------------------%
    for istate=1:nstates;
      solcost.initial.state(istate) = x0Pert(istate);
      [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
      DMayer.x0(istate) = (ComplexMayer-Mayer)/pertx0(istate);
      solcost.initial.state(istate) = x0(istate);
      solcost.terminal.state(istate) = xfPert(istate);
      [ComplexMayer,Lagrange] = feval(mysetup.funcs.cost,solcost);
      DMayer.xf(istate) = (ComplexMayer-Mayer)/pertxf(istate);
      solcost.terminal.state(istate) = xf(istate);
    end;
    %----------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Parameters %
    %----------------------------------------------------------%
    for iparameter=1:nparameters;
      solcost.parameter(iparameter) = parameterPert(iparameter);
      [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
      DMayer.parameter(iparameter) = (ComplexMayer-Mayer)./pertParameter(iparameter);
      solcost.parameter(iparameter) = parameters(iparameter);
    end;
    Cost = Cost + Mayer;
  end;
  %----------------------------------%
  % Get Derivatives of Lagrange Cost %
  %----------------------------------%
  if ~isempty(Lagrange)
    solcost.time = tRadauPert;
    [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
    DLagrange.time = (ComplexLagrange-Lagrange)./pertTime;
    solcost.time = t_radau;
    for istate=1:nstates
      solcost.state(:,istate) = stateRadauPert(:,istate);
      [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
      DLagrange.state(:,istate) = (ComplexLagrange-Lagrange)./pertState(:,istate);
      solcost.state(:,istate) = state_radau(:,istate);
    end;
    for icontrol=1:ncontrols
      solcost.control(:,icontrol) = controlRadauPert(:,icontrol);
      [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
      DLagrange.control(:,icontrol) = (ComplexLagrange-Lagrange)./pertControl(:,icontrol);
      solcost.control(:,icontrol) = control_radau(:,icontrol);
    end;
    for iparameter=1:nparameters
      solcost.parameter(iparameter) = parameterPert(iparameter);
      [ComplexMayer,ComplexLagrange] = feval(mysetup.funcs.cost,solcost);
      DLagrange.parameter(:,iparameter) = (ComplexLagrange-Lagrange)./pertParameter(iparameter);
      solcost.parameter(iparameter) = parameters(iparameter);
    end;
    integrand = tspan*mysetup.ps(i).Weights*Lagrange/2;
    Cost = Cost + integrand;     
  end;
  %-----------------------------------%
  % Insert derivatives of constraints %
  %-----------------------------------%
  nnodes = mysetup.nodes(i);
  nrows = (nstates+npaths)*nnodes+nevents;
  ncols = nstates*(nnodes+1)+ncontrols*nnodes+nparameters+2;
  ndaevar = (nstates+ncontrols+nparameters+2);
  nnonzeros = (nstates+npaths)*ndaevar+nevents*(2*nstates+nparameters+2);
  Jcon = spalloc(nrows,ncols,nnonzeros);
  rowshift = 0;
  for ii = 1:nstates
    rows = rowshift+nnodes*(ii-1)+1:rowshift+nnodes*ii;
    colshift = 0;
  %--------------%
  % Insert df/dx %
  %--------------%
    for jj = 1:nstates
      cols = colshift+(jj-1)*nnodes+jj:colshift+jj*nnodes+jj;
      cols = cols(1:end-1);
      if ii == jj
        Jcon(rows,cols) = mysetup.ps(i).Ddiag(:,1:end-1) - (tf-t0)/2*diag(DDae.state{jj}(:,ii));                
        % Jcon(rows,cols) = mysetup.ps(i).D-(tf-t0)/2*[diag(DDae.state{jj}(:,ii)) zeros(length(rows),1)];
      else
        Jcon(rows,cols) = - (tf-t0)/2*diag(DDae.state{jj}(:,ii));
        % Jcon(rows,cols) = -(tf-t0)/2*[diag(DDae.state{jj}(:,ii)) zeros(length(rows),1)];
      end
    end
  %--------------%
  % Insert df/du &
  %--------------%
    colshift = nstates*(nnodes+1);
    for jj = 1:ncontrols
      cols = colshift+(jj-1)*nnodes+1:colshift+jj*nnodes;
      Jcon(rows,cols) = -(tf-t0)/2*diag(DDae.control{jj}(:,ii));            
    end
    colshift = colshift + ncontrols*nnodes;
    %---------------%
    % Insert df/dto %
    %---------------%
    cols = colshift+1;
    Jcon(rows,cols) = 1/2*dae_out(:,ii)-(tf-t0)/2*DDae.time(:,ii).*(-mysetup.ps(i).Points/2+1/2);
    colshift = colshift+1;
    %---------------%
    % Insert df/dtf %
    %---------------%
    cols = colshift+1;
    Jcon(rows,cols) = -1/2*dae_out(:,ii)-(tf-t0)/2*DDae.time(:,ii).*(mysetup.ps(i).Points/2+1/2);
    colshift = colshift+1;
    %------------------%
    % Insert df/dparam %
    %------------------%    
    for jj=1:nparameters
      cols = colshift+jj;
      Jcon(rows,cols) = -(tf-t0)/2*DDae.parameter{jj}(:,ii);            
    end
  end
  rowshift = rowshift + nstates*nnodes;
  for ii = 1:npaths
    rows = rowshift+(ii-1)*nnodes+1:rowshift+ii*nnodes;
    colshift = 0;
    %--------------%
    % Insert dc/dx %
    %--------------%
    for jj = 1:nstates
      cols = colshift+(jj-1)*nnodes+jj:colshift+jj*nnodes+jj-1;
      Jcon(rows,cols) = diag(DDae.state{jj}(:,ii+nstates));
    end
    colshift = colshift + nstates*(nnodes+1);
    %--------------%
    % Insert dc/du %
    %--------------%
    for jj = 1:ncontrols
      cols = colshift+nnodes*(jj-1)+1:colshift+nnodes*jj;
      Jcon(rows,cols) = diag(DDae.control{jj}(:,ii+nstates));
    end
    colshift = colshift + ncontrols*nnodes;
    %---------------%
    % Insert dc/dto %
    %---------------%
    cols = colshift + 1;
    Jcon(rows,cols) = DDae.time(:,ii+nstates).*(-mysetup.ps(i).Points/2+1/2);
    colshift = colshift + 1;
    %---------------%
    % Insert dc/dtf %
    %---------------%
    cols = colshift + 1;
    Jcon(rows,cols) = DDae.time(:,ii+nstates).*(mysetup.ps(i).Points/2+1/2);
    colshift = colshift + 1;
    %------------------%
    % Insert dc/dparam %
    %------------------%
    for jj = 1:nparameters
      cols = colshift + jj;
      Jcon(rows,cols) = DDae.parameter{jj}(:,ii+nstates);
    end
  end
  rowshift = rowshift+npaths*nnodes;
  for ii = 1:nevents
    rows = rowshift+ii;
    %-----------%
    % dEvent/dx %
    %-----------%
    for jj = 1:nstates
      col0 = nnodes*(jj-1)+jj;
      colF = nnodes*jj+jj;
      Jcon(rows,col0) = DEvent.x0{jj}(ii);
      Jcon(rows,colF) = DEvent.xf{jj}(ii);
    end
    %------------%
    % dEvent/dt0 %
    %------------%
    cols = nstates*(nnodes+1)+ncontrols*nnodes+1;
    Jcon(rows,cols) = DEvent.t0(ii,1);
    %------------%
    % dEvent/dtf %
    %------------%
    cols = nstates*(nnodes+1)+ncontrols*nnodes+2;
    Jcon(rows,cols) = DEvent.tf(ii,1);
    %---------------%
    % dEvent/dparam %
    %---------------%
    for jj = 1:nparameters
      cols = nstates*(nnodes+1)+ncontrols*nnodes+2+jj;
      Jcon(rows,cols) = DEvent.parameter{jj}(ii);
    end
  end
  %------------------------------------------------%
  % Map constraint derivatives to sparsity pattern %
  %------------------------------------------------%
  J(mysetup.constraint_indices{i}+1,mysetup.variable_indices{i}) = Jcon;
  %----------------------%
  % Get Cost derivatives %
  %----------------------%
  Jcost = zeros(1,nstates*(nnodes+1)+ncontrols*nnodes+nparameters+2);
  %----------%
  % dCost/dx %
  %----------%
  for jj = 1:nstates
    col0 = nnodes*(jj-1)+jj;
    colF = nnodes*jj+jj;
    cols = col0:colF;
    cols = cols(1:end-1);
    Jcost(1,col0) = DMayer.x0(jj);
    Jcost(1,cols) = (tf-t0)*mysetup.ps(i).Weights.*DLagrange.state(:,jj).'/2;
    Jcost(1,colF) = DMayer.xf(jj);
  end
  colshift = nstates*(nnodes+1);
  %----------%
  % dCost/du %
  %----------%
  for jj = 1:ncontrols
    cols = colshift+(jj-1)*nnodes+1:colshift+jj*nnodes;
    Jcost(1,cols) = (tf-t0)*mysetup.ps(i).Weights.*DLagrange.control(:,jj).'/2;
  end
  colshift = colshift + ncontrols*nnodes;
  %-----------%
  % dCost/dt0 %
  %-----------%
  cols = colshift + 1;
  Jcost(1,cols) = DMayer.t0-mysetup.ps(i).Weights*Lagrange/2+(tf-t0)/2*mysetup.ps(i).Weights*diag(DLagrange.time)*(-mysetup.ps(i).Points/2+1/2);
  %-----------%
  % dCost/dtf %
  %-----------%  
  colshift = colshift+1;
  cols = colshift + 1;
  Jcost(1,cols) = DMayer.tf+mysetup.ps(i).Weights*Lagrange/2+(tf-t0)/2*mysetup.ps(i).Weights*diag(DLagrange.time)*(mysetup.ps(i).Points/2+1/2);
  %--------------%
  % dCost/dparam %
  %--------------%
  colshift = colshift + 1;
  for jj = 1:nparameters
    cols = colshift+jj;
    Jcost(1,cols) = DMayer.parameter(:,jj) + (tf-t0)*mysetup.ps(i).Weights*DLagrange.parameter(:,jj)/2;
  end
  %------------------------------------------%
  % Map cost derivatives to sparsity pattern %
  %------------------------------------------%
  J(1,mysetup.variable_indices{i}) = Jcost;
  %----------------------------------%
  % save solcost for use in linkages %
  %----------------------------------%
  solTotal(i) = solcost;
end;
Constraints = vertcat(Cons{:,1});
linkages = mysetup.linkages;
if ~isempty(mysetup.linkages),
  link_out = cell(mysetup.numlinkpairs,1);
  link_row = length(Constraints) + 1;
  for ipair=1:mysetup.numlinkpairs;
    nlinks = length(linkages(ipair).min);
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
    %-------------------------%
    % Get Linkage derivatives %
    %-------------------------%
    link_row = link_row(end)+ (1:nlinks);
    nstatesLeft = mysetup.sizes(left_phase,1);
    ncontrolsLeft = mysetup.sizes(left_phase,2);
    nparametersLeft = mysetup.sizes(left_phase,3);
    nstatesRight = mysetup.sizes(right_phase,1);
    ncontrolsRight = mysetup.sizes(right_phase,2);
    nparametersRight = mysetup.sizes(right_phase,3);
    nnodesLeft = mysetup.limits(left_phase).nodes;
    nnodesRight = mysetup.limits(right_phase).nodes;
    DLink.xf_left = zeros(length(link_row),nstatesLeft);
    DLink.p_left = zeros(length(link_row),nparametersLeft);
    DLink.x0_right = zeros(length(link_row),nstatesRight);
    DLink.p_right = zeros(length(link_row),nparametersRight);
    pertxfLeft= pertFactor*(1+abs(xf_left));
    pertpLeft= pertFactor*(1+abs(p_left));
    pertx0Right= pertFactor*(1+abs(x0_right));
    pertpRight= pertFactor*(1+abs(p_right));
    xfLeftPert = xf_left+pertxfLeft;
    pLeftPert = p_left+pertpLeft;
    x0RightPert = x0_right+pertx0Right;
    pRightPert = p_right+pertpRight;
    %----------------------------------------------------%
    % Get Derivatives of Linkage with Respect to xf_left %
    %----------------------------------------------------%
    for ii=1:nstatesLeft
      sollink.left.state(ii) = xfLeftPert(ii);
      DLink.xf_left(:,ii) = (feval(mysetup.funcs.link,sollink)-link_out{ipair,1})/pertxfLeft(ii);
      sollink.left.state(ii) = xf_left(ii);
    end;
    for ii=1:nparametersLeft
      sollink.left.parameter(ii) = pLeftPert(ii);
      DLink.p_left(:,ii) = (feval(mysetup.funcs.link,sollink)-link_out{ipair,1})/pertpLeft(ii);
      sollink.left.parameter(ii) = p_left(ii);
    end;
    for ii=1:nstatesRight
      sollink.right.state(ii) = x0RightPert(ii);
      DLink.x0_right(:,ii) = (feval(mysetup.funcs.link,sollink)-link_out{ipair,1})/pertx0Right(ii);
      sollink.right.state(ii) = x0_right(ii);
    end;
    for ii=1:nparametersRight
      sollink.right.parameter(ii) = pRightPert(ii);
      DLink.p_right(:,ii) = (feval(mysetup.funcs.link,sollink)-link_out{ipair,1})/pertpRight(ii);
      sollink.right.parameter(ii) = p_right(ii);
    end;
    for ii=1:nlinks
      for jj = 1:nstatesLeft
        colFLeft = sum(horzcat(mysetup.variable_indices{1:left_phase-1}))+nnodesLeft*jj+jj;
        Jcon(link_row(ii),colFLeft) = DLink.xf_left(ii,jj);
      end;
      for jj = 1:nstatesRight
        col0Right = sum(horzcat(mysetup.variable_indices{1:right_phase-1}))+nnodesRight*(jj-1)+jj;
        Jcon(link_row(ii),col0Right) = DLink.x0_right(ii,jj);
      end;
      for jj = 1:nparametersLeft
        colFLeft = sum(horzcat(mysetup.variable_indices{1:left_phase-1}))+nstatesLeft*(nnodesLeft+1)+ncontrolsLeft*nnodesLeft+jj;
        Jcon(link_row(ii),colFLeft) = DLink.p_left(ii,jj);
      end;
      for jj = 1:nparametersRight
        col0Right = sum(horzcat(mysetup.variable_indices{1:right_phase-1}))+nstatesRight*(nnodesRight+1)+ncontrolsRight*nnodesRight+jj;
        Jcon(link_row(ii),col0Right) = DLink.p_right(ii,jj);
      end;
      [rows,cols] = find(mysetup.sparsity_all(link_row,:));
      Dlink_out = [DLink.xf_left DLink.p_left DLink.x0_right DLink.p_right];
      J(link_row,1:max(cols)) = sparse(rows,cols,Dlink_out(:));
    end
  end
  Clink = vertcat(link_out{:,1});
  Constraints = [Constraints; Clink];
end
% Constraints = Constraints-mysetup.Alinear_augmented(2:mysetup.numnonlin+1,:)*mysetup.invDx*mysetup.column_shifts;
C = [Cost; Constraints; mysetup.initlincons];
