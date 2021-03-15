function [C J] = gpopsObjandConsAnalytic(x)
%------------------------------------------------------------------%
% Compute the nonlinear constraints and objective function using   %
% analytic derivatives.                                            %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

global mysetup

[ncon, nvar] = size(mysetup.sparsity_all);
J = spalloc(ncon,nvar,length(mysetup.iGfun));
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
  [daeOut,DdaeOut] = feval(mysetup.funcs.dae,sol);
  %----------------------------------------------%
  % Compute Analytic Derivatives of DAE Function %
  %----------------------------------------------%
  nnodes = mysetup.limits(i).nodes;
  DDae.time  = zeros(nnodes,nstates+npaths);
  DDae.state = cell(1,nstates);
  DDae.control = cell(1,ncontrols);
  DDae.parameter = cell(1,nparameters);
  %------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to State %
  %------------------------------------------------------%
  for istate=1:nstates
    DDae.state{istate}=reshape(DdaeOut(:,istate),nnodes,nstates+npaths);
  end;
  %---------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to Control. %
  %---------------------------------------------------------%
  for icontrol=1:ncontrols
    DDae.control{icontrol}=reshape(DdaeOut(:,icontrol+nstates),nnodes,nstates+npaths);
  end;
  %------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to Time. %
  %------------------------------------------------------%
  DDae.time = reshape(DdaeOut(:,nstates+ncontrols+1),nnodes,nstates+npaths);
  %-----------------------------------------------------------%
  % Compute the Derivative of DAEs with Respect to Parameter. %
  %-----------------------------------------------------------%
  for iparameter=1:nparameters
    DDae.parameter{iparameter}=reshape(DdaeOut(:,iparameter+nstates+ncontrols+1),nnodes,nstates+npaths);
  end;
  %-----------------%
  % Get constraints %
  %-----------------%
  odelhs  = mysetup.ps(i).Ddiag*state_matrix;
  oderhs  = tspan*daeOut(:,1:nstates)/2;
  if npaths > 0
    paths = daeOut(:,nstates+1:end);
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
    [eventOut,DeventOut] = feval(mysetup.funcs.event,solevents);
    %-------------------------------------------------------%
    % Get Derivatives of Event with Respect to Initial Time %
    %-------------------------------------------------------%
    DEvent.t0 = DeventOut(:,nstates+1);
    %--------------------------------------------------------%
    % Get Derivatives of Event with Respect to Terminal Time %
    %--------------------------------------------------------%
    DEvent.tf = DeventOut(:,2*nstates+2);
    %---------------------------------------------------------------------%
    % Get Derivatives of Event with Respect to Initial and Terminal State %
    %---------------------------------------------------------------------%
    for istate=1:nstates;
      DEvent.x0{istate} = DeventOut(:,istate);
      DEvent.xf{istate} = DeventOut(:,istate+nstates+1);
    end
    %-----------------------------------------------------%
    % Get Derivatives of Event with Respect to Parameters %
    %-----------------------------------------------------%
    for iparameter=1:nparameters;
      DEvent.parameter{iparameter} = DeventOut(:,iparameter+2*(nstates+1));
    end;
  else
    DEvent = [];
  end;
  if npaths >0 && nevents > 0
    Cons{i,1} = [defects(:); paths(:); eventOut];
  elseif npaths>0
    Cons{i,1} = [defects(:); paths(:)];
  elseif nevents>0
    Cons{i,1} = [defects(:); eventOut];
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
  [MayerOut,LagrangeOut,DMayerOut,DLagrangeOut] = feval(mysetup.funcs.cost,solcost);
  if ~isempty(MayerOut),
    %------------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Initial Time %
    %------------------------------------------------------------%
    DMayer.t0 = DMayerOut(:,nstates+1);
    %-------------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Terminal Time %
    %-------------------------------------------------------------%
    DMayer.tf = DMayerOut(:,2*nstates+2);
    %-------------------------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Initial and Terminal Time %
    %-------------------------------------------------------------------------%
    for istate=1:nstates;
      DMayer.x0(istate) = DMayerOut(:,istate);
      DMayer.xf(istate) = DMayerOut(:,istate+nstates+1);
    end
    %----------------------------------------------------------%
    % Get Derivatives of Mayer Cost with Respect to Parameters %
    %----------------------------------------------------------%
    for iparameter=1:nparameters;
      DMayer.parameter(iparameter) = DMayerOut(:,iparameter+2*(nstates+1));
    end;
    Cost = Cost + MayerOut;
  end;
  %----------------------------------%
  % Get Derivatives of Lagrange Cost %
  %----------------------------------%
  if ~isempty(LagrangeOut)
    DLagrange.time = DLagrangeOut(:,end);
    for istate=1:nstates
      DLagrange.state(:,istate) = DLagrangeOut(:,istate);
    end;
    for icontrol=1:ncontrols
      DLagrange.control(:,icontrol) = DLagrangeOut(:,icontrol+nstates);
    end;
    for iparameter=1:nparameters
      DLagrange.parameter(:,iparameter) = DLagrangeOut(:,iparameter+nstates+ncontrols);
    end;
    integrand = tspan*mysetup.ps(i).Weights*LagrangeOut/2;
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
    % Insert df/dx &
    %--------------%
    for jj = 1:nstates
      cols = colshift+(jj-1)*nnodes+jj:colshift+jj*nnodes+jj;
      cols = cols(1:end-1);
      if ii == jj
        Jcon(rows,cols) = mysetup.ps(i).Ddiag(:,1:end-1) - (tf-t0)/2*diag(DDae.state{jj}(:,ii));                
      else
        Jcon(rows,cols) = - (tf-t0)/2*diag(DDae.state{jj}(:,ii));
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
    % Insert df/dto &
    %---------------%
    cols = colshift+1;
    Jcon(rows,cols) = 1/2*daeOut(:,ii)-(tf-t0)/2*DDae.time(:,ii).*(-mysetup.ps(i).Points/2+1/2);
    colshift = colshift+1;
    %---------------%
    % Insert df/dtf &
    %---------------%
    cols = colshift+1;
    Jcon(rows,cols) = -1/2*daeOut(:,ii)-(tf-t0)/2*DDae.time(:,ii).*(mysetup.ps(i).Points/2+1/2);
    colshift = colshift+1;
    %------------------%
    % Insert df/dparam &
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
    % Insert dc/dparam &
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
    % dEvent/dx &
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
  %------------------------------------------------&
  % Map constraint derivatives to sparsity pattern &
  %------------------------------------------------&
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
  cols = colshift+1;
  Jcost(1,cols) = DMayer.t0-mysetup.ps(i).Weights*LagrangeOut/2+(tf-t0)/2*mysetup.ps(i).Weights*diag(DLagrange.time)*(-mysetup.ps(i).Points/2+1/2);
  %-----------%
  % dCost/dtf %
  %-----------%
  colshift = colshift+1;
  cols = colshift+1;
  Jcost(1,cols) = DMayer.tf+mysetup.ps(i).Weights*LagrangeOut/2+(tf-t0)/2*mysetup.ps(i).Weights*diag(DLagrange.time)*(mysetup.ps(i).Points/2+1/2);
  %--------------%
  % dCost/dparam %
  %--------------%
  colshift = colshift+1;
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
  linkOut = cell(mysetup.numlinkpairs,1);
  link_row = length(Constraints) + 1;
  for ipair=1:mysetup.numlinkpairs;
    nlinks = length(linkages(ipair).min);
    left_phase  = linkages(ipair).left.phase;
    right_phase = linkages(ipair).right.phase;
    xf_left = solTotal(left_phase).terminal.state;
    p_left  = solTotal(left_phase).parameter;
    x0_right = solTotal(right_phase).initial.state;
    p_right  = solTotal(right_phase).parameter;
    %----------------------------%
    % Get Connection constraints %
    %----------------------------%
    sollink.left.state = xf_left;
    sollink.left.parameter = p_left;
    sollink.left.phase = left_phase;
    sollink.right.state = x0_right;
    sollink.right.parameter = p_right;
    sollink.right.phase = right_phase;
    %-----------------------------------------%
    % Get Linkage Constraints and Derivatives %
    %-----------------------------------------%
    [linkOut{ipair,1},DLinkOut] = feval(mysetup.funcs.link,sollink);
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
    DLinkLeft = DLinkOut(:,1:nstatesLeft+nparametersLeft);
    indexStart = nstatesLeft+nparametersLeft;
    DLinkRight = DLinkOut(:,indexStart+1:indexStart+nstatesRight+nparametersRight);
    DLink.xf_left = DLinkLeft(:,1:nstatesLeft);
    DLink.p_left = DLinkLeft(:,nstatesLeft+1:nstatesLeft+nparametersLeft);
    DLink.x0_right = DLinkRight(:,1:nstatesRight);
    DLink.p_right = DLinkRight(:,nstatesRight+1:nstatesRight+nparametersRight);
    for ii=1:nlinks
      cols = mysetup.indices(left_phase).state(nnodesLeft+1:nnodesLeft+1:end);
      Jcon(link_row,cols) = DLink.xf_left;
      cols = mysetup.indices(right_phase).state(nnodesRight+1:nnodesRight+1:end);
      Jcon(link_row,cols) = DLink.x0_right;
      cols = mysetup.indices(left_phase).parameter;
      Jcon(link_row,cols) = DLink.p_left;
      cols = mysetup.indices(right_phase).parameter;
      Jcon(link_row,cols) = DLink.p_right;
      [rows,cols] = find(mysetup.sparsity_all(link_row,:));
      DlinkOut = [DLink.xf_left DLink.p_left DLink.x0_right DLink.p_right];
      J(link_row,1:max(cols)) = sparse(rows,cols,DlinkOut(:));
    end
  end
  Clink = vertcat(linkOut{:,1});
  Constraints = [Constraints; Clink];
end
C = [Cost; Constraints; mysetup.initlincons];
