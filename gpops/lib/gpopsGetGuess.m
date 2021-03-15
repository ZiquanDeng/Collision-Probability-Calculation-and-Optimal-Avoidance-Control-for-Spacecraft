function setup = gpopsGetGuess(setup)

%--------------------------------------------------------------%
% Calculuate the guess used by the NLP solver in solving a     %
% multiple-phase optimal control problem using the hp-adaptive %
% Radau pseudospectral method                                  %
%--------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David A. Benson          %
%--------------------------------------------------------------%
global igrid

guess = setup.guess;
numphases = setup.numphases;
if length(guess)<numphases
  error('Number of phases in Guess is less than number of phases in limits')
end
if length(guess)>numphases,
  guess = guess(1:numphases);
end;
sizes = setup.sizes;
nodes = setup.nodes;
nlpGuess = cell(numphases,1);
ps = repmat(struct('D',[],'Points',[],'Weights',[],'Ddiag',[]),1,numphases);
for iphase=1:numphases;  
  nstates     = sizes(iphase,1);
  ncontrols   = sizes(iphase,2);
  nparameters = sizes(iphase,3);
  npaths      = sizes(iphase,4);
  nevents     = sizes(iphase,5);
  tGuess      = guess(iphase).time;
  xGuess      = guess(iphase).state;
  %--------------------------------------------%
  % Get the guess in each phase of the problem %
  %--------------------------------------------%
  switch isfield(guess(iphase),'control')
   case true
    uGuess    = guess(iphase).control;
   case false
    uGuess = [];
  end
  switch isfield(guess(iphase),'parameter'),
   case true
    pGuess    = guess(iphase).parameter;
   case false
    pGuess = [];
  end;
  %-------------------------------%
  % Check guess for proper format %
  %-------------------------------%
  if igrid == 1
    [length_guess col_t] = size(tGuess);
    if col_t ~= 1
      error('Guess for time in phase %i must be a column vector',iphase)
    elseif length_guess < 2
      error('Guess in phase %i must have a least two points',iphase)
    elseif length_guess ~= length(unique(tGuess))
      error('Guess for time in phase %i does not contain unique values',iphase);
    end
    [row_state col_state] = size(xGuess);
    if row_state ~= length_guess
      error('Size of state does not match size of time in guess for phase %i',iphase) 
    elseif col_state ~= nstates
      error('Number of states in guess does not match limits in phase %i',iphase) 
    end
    [row_control col_control] = size(uGuess);
    if row_control ~= 0 && row_control ~= length_guess
      error('Size of control does not match size of time in guess for phase %i',iphase) 
    elseif col_control ~= ncontrols
      error('Number of controls in guess does not match limits in phase %i',iphase) 
    end
    [row_param col_param] = size(pGuess);
    if col_param ~= 0 && col_param ~= 1
      error('Guess for parameters in phase %i must be a column vector',iphase)
    elseif row_param ~= nparameters
      error('Number of parameters in guess does not match limits in phase %i',iphase) 
    end
    %---------------------------------------%
    % Check Cost function for proper format %
    %---------------------------------------%
    if ~isempty(setup.funcs.cost),
      if ~(isa(setup.funcs.cost,'char') || isa(setup.funcs.cost,'function_handle'))
        error('Invalid Cost function in setup.funcs.cost')
      end
      clear sol
      sol.initial.time = tGuess(1);
      sol.initial.state = xGuess(1,:).';
      sol.terminal.time = tGuess(end);
      sol.terminal.state = xGuess(end,:).';
      sol.time = tGuess;
      sol.state = xGuess;
      sol.control = uGuess;
      sol.parameter = pGuess;
      sol.phase = iphase;
      %-----------------------------------%
      % Evaluate the user's cost function %
      %-----------------------------------%
      [Mayer, Lagrange] = feval(setup.funcs.cost,sol);
      %-------------------------------------------------%
      % Check to see if the outputs are the proper size %
      %-------------------------------------------------%
      if ~isscalar(Mayer)
        fname = setup.funcs.cost; if isa(fname,'function_handle'); fname = func2str(fname); end
        error('Cost function "%s" did not return scalar for Mayer cost using guess in phase %i',fname,iphase)
      end
      [row_L col_L] = size(Lagrange);
      if row_L ~= length_guess || col_L ~= 1
        fname = setup.funcs.cost; if isa(fname,'function_handle'); fname = func2str(fname); end
        error('Cost function "%s" returned invalid size of column vector for Lagrange cost using guess in phase %i',fname,iphase)
      end
    else
      error('Must Specify a Cost Function for Problem');
    end;
    %---------------------------------------------------------%
    % Check Differential-Algebraic function for proper format %
    %---------------------------------------------------------%
    switch isfield(setup.funcs','dae')
      case true
       if ~isempty(setup.funcs.dae),
         if ~(isa(setup.funcs.dae,'char') || isa(setup.funcs.dae,'function_handle'))
           error('Invalid Dae function in setup.funcs.dae')
         end;
       end
       clear sol
       sol.time      = tGuess;
       sol.state     = xGuess;
       sol.control   = uGuess;
       sol.parameter = pGuess; 
       sol.phase     = iphase;
       %-----------------------------------------------------%
       % Evaluate User DAE Function and Test Size of Outputs %
       %-----------------------------------------------------%
       dae = feval(setup.funcs.dae,sol);
       [row_dae col_dae] = size(dae);
       if col_dae ~= (nstates+npaths)
         fname = setup.funcs.dae; if isa(fname,'function_handle'); fname = func2str(fname); end
         error('Dae function "%s" returned invalid number of columns using guess in phase %i',fname,iphase) 
       elseif row_dae ~= length_guess
         fname = setup.funcs.dae; if isa(fname,'function_handle'); fname = func2str(fname); end
         error('Dae function "%s" returned invalid number of rows using guess in phase %i',fname,iphase)
       end
     case false
      if (nstates>0),
        error('Must Specify a Differential-Algebraic Function When # of States is Nonzero');
      end
    end
    % ---------------------------------------%
    % Check Event function for proper format %
    % ---------------------------------------%
    if (nevents > 0) || isfield(setup.funcs,'event'),
      switch isempty(setup.funcs.event),
       case false
        if ~(isa(setup.funcs.event,'char') || isa(setup.funcs.event,'function_handle'))
          error('Invalid Event function in setup.funcs.event')
        end
        clear sol
        sol.initial.time = tGuess(1);
        sol.initial.state = xGuess(1,:).';
        sol.terminal.time = tGuess(2);
        sol.terminal.state = xGuess(end,:).';
        sol.parameter = pGuess;
        sol.phase = iphase;
        %-------------------------------------------------------%
        % Evaluate User Event Function and Test Size of Outputs %
        %-------------------------------------------------------%
        event = feval(setup.funcs.event,sol);
        [row_event col_event] = size(event);
        if nevents > 0 && col_event ~= 1
          fname = setup.funcs.event; if isa(fname,'function_handle'); fname = func2str(fname); end
          error('Event function "%s" must return column vector in phase %i',fname,iphase)
        elseif row_event ~= nevents
          fname = setup.funcs.event; if isa(fname,'function_handle'); fname = func2str(fname); end
          error('Event function "%s" returned invalid number of events using guess in phase %i',fname,iphase)
        end
      case true
       if (nevents>0),
         error('Must Specify Event Function When # of Events is Nonzero');
       end;
      end;
    end % if (nevents > 0) || isfield(setup.funcs,'event'),
  end % if igrid == 1
  % ------------------%
  % Interpolate Guess %
  % ------------------%
  nodesPerInterval = setup.limits(iphase).nodesPerInterval;
  meshPoints = setup.limits(iphase).meshPoints;
  numseg = length(setup.limits(iphase).nodesPerInterval);
  t0Guess            = tGuess(1);
  tfGuess            = tGuess(end);
  RPM                = gpopsRPM(numseg,meshPoints,nodesPerInterval);
  ps(iphase).D       = RPM.differentiationMatrix;
  ps(iphase).Points  = RPM.Points;
  ps(iphase).Weights = RPM.Weights.';
  ps(iphase).Ddiag   = RPM.differentiationMatrixDiag;
  numRadau           = length(RPM.Points);
  tau_plus_ends = [RPM.Points;1];
  if isequal(t0Guess,tfGuess),
    %-----------------------------------------------------------%
    % If the initial and final times in the phase are the same, %
    % then arbitrarily interpolate on an evenly spaced grid of  %
    % of length tGuess between -1 and +1.                       %
    %-----------------------------------------------------------%
    tauGuess = linspace(-1,+1,length(tGuess));
  else
    %-----------------------------------------------------------%
    % If the initial and terminal times are not the same, then  %
    % transform the time to the interval between -1 and +1.     %
    %-----------------------------------------------------------%
    tauGuess = 2*(tGuess-t0Guess)/(tfGuess-t0Guess)-1;
  end;  
  if nstates>0,
    xinterp = interp1(tauGuess,xGuess,tau_plus_ends,'spline','extrap');
  else
    xinterp = [];
  end;
  if ncontrols>0,
    uinterp = interp1(tauGuess,uGuess,tau_plus_ends(1:end-1),'spline','extrap');
  else
    uinterp = [];
  end; 
  nlpGuess{iphase,1} = [xinterp(:); uinterp(:); t0Guess; tfGuess; pGuess];
  % ---------------------------------------------%
  % Check DAE function for proper use of t0 & tf %
  % ---------------------------------------------%
  if igrid == 1 % check only on first iteration
    clear sol
    sol.time      = (tfGuess-t0Guess)*(RPM.Points+1)/2+t0Guess;
    sol.state     = xinterp(1:end-1,:);
    sol.control   = uinterp;
    sol.parameter = pGuess; 
    sol.phase     = iphase;
    %--------------------------------%
    % Evaluate the User DAE Function %
    %--------------------------------%
    dae = feval(setup.funcs.dae,sol);
    %----------------------------%
    % Find the State Perurbation %
    %----------------------------%
    limitRange = (setup.limits(iphase).state.max(:,2) - setup.limits(iphase).state.min(:,2));
    limitRange = max(limitRange,1E-6);
    statePert = 0.01*limitRange;
    limitMean = (setup.limits(iphase).state.max(:,2) + setup.limits(iphase).state.min(:,2))/2;
    signPert = sign(sol.state - repmat(limitMean.',length(sol.time),1));
    signPert(logical(signPert == 0)) = 1;
    statePert = signPert .* repmat(statePert',length(sol.time),1);
    %-------------------------------%
    % Find the Control Perturbation %
    %-------------------------------%
    if ncontrols > 0
      limitRange = (setup.limits(iphase).control.max - setup.limits(iphase).control.min);
      limitRange = max(limitRange,1E-6);
      controlPert = 0.01*limitRange;
      limitMean = (setup.limits(iphase).control.max + setup.limits(iphase).control.min)/2;
      signPert = sign(sol.control - repmat(limitMean.',length(sol.time),1));
      signPert(logical(signPert == 0)) = 1;
      controlPert = signPert .* repmat(controlPert',length(sol.time),1);
    end
    %----------------------------%
    % Find the Time Perturbation %
    %----------------------------%
    timePert = (tfGuess-t0Guess)*(-RPM.Points/10)/2;
    timePert(timePert == 0) = min(abs(timePert(abs(timePert)>0)))/2;
    %-------------------------%
    % Perturb Each Time Point %
    %-------------------------%
    for j = 1:length(sol.time)
      solPert = sol;
      solPert.time(j) = solPert.time(j) + timePert(j);
      solPert.state(j,:) = solPert.state(j,:) + statePert(j,:);
      if ncontrols > 0
        solPert.control(j,:) = solPert.control(j,:) + controlPert(j,:);
      end
      daePert = feval(setup.funcs.dae,solPert);
      % check for changes in non-j row
      daeErr = daePert - dae;
      daeErr(j,:) = 0;
      if any(any(daeErr))
        fname = setup.funcs.dae; if isa(fname,'function_handle'); fname = func2str(fname); end
        error(['Invalid use of state, control, or time in dae function (%s) in phase: %i\n',...
               'dae function at time j may only depend on state, control, and time at time j.\n',...
               'If delta state or time from a boundary (t0,tf) is needed in dae equations,\n',...
               'add parameters to the problem and an event to constrain the parameters \n',...
               'to be equal to the desired boundary state or time, then use the parameters \n',...
               'within the dae function.'], fname, iphase)
      end
    end
  end
end;
% init_vector = vertcat(init{:,1});
% nlpGuess = vertcat(nlpGuess{:,1});
% setup.init_vector = init_vector;
setup.nlpGuess = vertcat(nlpGuess{:,1});
setup.ps = ps;
%------------------------------------------%
% Check Linkage function for proper format %
%------------------------------------------%
if igrid == 1
  numlinks = setup.numlinks;
  linkages = setup.linkages;
  numlinkpairs = setup.numlinkpairs;
  if numlinks > 0 
    if ~(ischar(setup.funcs.link) || isa(setup.funcs.link,'function_handle')),
      error('Invalid Linkage function in setup.funcs.link')
    end
    for ipair = 1:numlinkpairs;
      leftPhase = linkages(ipair).left.phase;
      rightPhase = linkages(ipair).right.phase;
      nlink = length(setup.linkages(ipair).min);
      clear sol
      sol.left.state       = guess(leftPhase).state(end,:).';
      if isfield(guess(leftPhase),'parameter'),
        sol.left.parameter   = guess(leftPhase).parameter;
      else
        sol.left.parameter   = [];
      end
      sol.left.phase        = leftPhase;
      sol.right.state      = guess(rightPhase).state(1,:).';
      if isfield(guess(rightPhase),'parameter'),
        sol.right.parameter   = guess(rightPhase).parameter;
      else
        sol.right.parameter   = [];
      end
      sol.right.phase       = rightPhase;
      %---------------------------%
      % Try user linkage function %
      %---------------------------%
      link = feval(setup.funcs.link,sol);
      %-----------------------%
      % Check size of outputs %
      %-----------------------%
      [row_con col_con] = size(link);
      if nlink > 0 && col_con ~= 1
        fname = setup.funcs.link; if isa(fname,'function_handle'); fname = func2str(fname); end
        error('Linkage function "%s" must return column vector for linkage %i',fname,ipair)
      elseif row_con ~= nlink
        fname = setup.funcs.link; if isa(fname,'function_handle'); fname = func2str(fname); end
        error('Linkage function "%s" returned invalid number of linkage constraints for linkage %i',fname,ipair)
      end
    end
  else
    if isfield(setup.funcs,'link'),
      if ~isempty(setup.funcs.link),
        if ischar(setup.funcs.link) || isa(setup.funcs.link,'function_handle')
          error('Linkage function defined in setup.funcs.link with no linkages')
        end
      end;
    end;
  end
end % if igrid == 1
