function setup = gpopsScaleNlp(setup)
%------------------------------------------------------------------%
% Determine the row and column scales for a non-sequential         %
% multiple-phase optimal control problem                           %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

nodes = setup.nodes;
sizes = setup.sizes;
constraint_indices = setup.constraint_indices;
limits = setup.limits;
xmin = setup.varbounds_min;
xmax = setup.varbounds_max;

numvars = setup.numvars;
numnonlin = setup.numnonlin;
numphases = setup.numphases;
threshTol = 1e-3;

colscales     = ones(numvars,1);
colshifts      = zeros(numvars,1);
rowscales     = ones(numnonlin,1);
dependencies  = cell(1,numphases);
stateScales   = cell(1,numphases);
controlScales = cell(1,numphases);
odeScales     = cell(1,numphases);
pathScales    = cell(1,numphases);
t0Scales      = cell(1,numphases);
tfScales      = cell(1,numphases);
parameterScales      = cell(1,numphases);
if isfield(setup,'autoscale') && isequal(setup.autoscale,'on')
  %------------------------------------------------------------%
  % Compute the Column Scales Using Maximum and Minimum Values % 
  % of Time, State, Control, and Parameters Provided by User   % 
  %------------------------------------------------------------%
  numvarsCum = 0;
  for iphase=1:numphases
    numLGR = length(setup.ps(iphase).Points);
    nstates = sizes(iphase,1);
    ncontrols = sizes(iphase,2);
    nparameters = sizes(iphase,3);
    numvarsCurr = nstates*(numLGR+1)+ncontrols*numLGR+nparameters+2;
    indicesCurr = 1:numvarsCurr;
    indicesCurr = indicesCurr+numvarsCum;
    tmincurr = setup.limits(iphase).time.min;
    tmaxcurr = setup.limits(iphase).time.max;
    xmincurr = setup.limits(iphase).state.min;
    xmaxcurr = setup.limits(iphase).state.max;
    umincurr = setup.limits(iphase).control.min;
    umaxcurr = setup.limits(iphase).control.max;
    pmincurr = setup.limits(iphase).parameter.min;
    pmaxcurr = setup.limits(iphase).parameter.max;
    xscale = ones(nstates,1);
    xshift = zeros(nstates,1);
    uscale = ones(ncontrols,1);
    ushift = zeros(ncontrols,1);
    tscale = 1;
    stateScalesCurr = ones(numLGR+1,nstates);
    stateShiftsCurr = zeros(numLGR+1,nstates);
    for istate=1:nstates
      if ~isequal(xmincurr(istate,2),xmaxcurr(istate,2)),
        xscalecurr = xmaxcurr(istate,2)-xmincurr(istate,2);
        xshiftcurr = 0.5-xmaxcurr(istate,2)/(xmaxcurr(istate,2)-xmincurr(istate,2));
        stateScalesCurr(:,istate) = xscalecurr*stateScalesCurr(:,istate);
        stateShiftsCurr(:,istate) = xshiftcurr*ones(numLGR+1,1);
      else
        xscalecurr = 1;
        xshiftcurr = 0;
      end;
    end;
    controlScalesCurr = ones(numLGR,ncontrols);
    controlShiftsCurr = zeros(numLGR,ncontrols);
    for icontrol=1:ncontrols
      if ~isequal(umincurr(icontrol),umaxcurr(icontrol)),
        uscalecurr = umaxcurr(icontrol)-umincurr(icontrol);
        ushiftcurr = 0.5-umaxcurr(icontrol)/(umaxcurr(icontrol)-umincurr(icontrol));
        controlScalesCurr(:,icontrol) = uscalecurr*controlScalesCurr(:,icontrol);
        controlShiftsCurr(:,icontrol) = ushiftcurr*ones(numLGR,1);
      else
        uscalecurr = 1;
        ushiftcurr = 0;
      end;
    end;
    tScaleCurr = 1;
    tshiftCurr = 0;
    tscaleCurr = tmaxcurr(2)-tmincurr(1);
    tshiftCurr = 0.5-tmaxcurr(2)/(tmaxcurr(2)-tmincurr(1));
    paramScalesCurr = ones(nparameters,1);
    paramShiftsCurr = zeros(nparameters,1);
    for iparam=1:nparameters
      if ~isequal(pmincurr(iparam),pmaxcurr(iparam)),
        paramScaleCurr(iparam) = pmaxcurr(iparam)-pmincurr(iparam);
        paramShiftsCurr(iparam) = 0.5-pmaxcurr(iparam)/(pmaxcurr(iparam)-pmincurr(iparam));
      else
        paramScalesCurr(iparam) = 1;
        paramShiftsCurr(iparam) = 0;
      end;
    end;
    colscales(indicesCurr) = [stateScalesCurr(:);controlScalesCurr(:);tscaleCurr;tscaleCurr;paramScalesCurr];
    colshifts(indicesCurr) = [stateShiftsCurr(:);controlShiftsCurr(:);tshiftCurr;tshiftCurr;paramShiftsCurr];
    numvarsCum = numvarsCum + numvarsCurr;
  end;
  %------------------------------------------------------------------------%
  % Now Compute the Row Scales Using the Method on Page 166 of Betts' Book %
  %------------------------------------------------------------------------%
  BIG_NUM = 1E20;
  % rand('state',0);
  ntrials = 93;
  extras.setup = setup;
  xupall = cell(1,numphases);
  phasecons = zeros(1,numphases);
  for iphase=1:numphases;
    tlowlimits = limits(iphase).time.min;
    tupplimits = limits(iphase).time.max;
    xlowlimits = limits(iphase).state.min;
    xupplimits = limits(iphase).state.max;
    ulowlimits = limits(iphase).control.min;
    uupplimits = limits(iphase).control.max;
    plowlimits = limits(iphase).parameter.min;
    pupplimits = limits(iphase).parameter.max;
    nstates = sizes(iphase,1);
    ncontrols = sizes(iphase,2);
    nparameters = sizes(iphase,3);
    npaths = sizes(iphase,4);
    nevents = sizes(iphase,5);
    colscales_curr = colscales(setup.variable_indices{1,iphase});
    state_scales = colscales_curr(1:(nodes(iphase)+1)*nstates);
    control_scales = colscales_curr((nodes(iphase)+1)*nstates+1:(nodes(iphase)+1)*nstates+nodes(iphase)*ncontrols);
    state_scales_matrix = reshape(state_scales,nodes(iphase)+1,nstates);
    control_scales_matrix = reshape(control_scales,nodes(iphase),ncontrols);
    state_scales_use = state_scales_matrix(2,:);
    control_scales_use = control_scales_matrix(1,:);
    stateScales{iphase} = state_scales_use;
    controlScales{iphase} = control_scales_use;
    odeScales{iphase} = state_scales_use;
    parameterIndexStart = nstates*(nodes(iphase)+1)+ncontrols*nodes(iphase)+1;
    parameterIndexFinish = nstates*(nodes(iphase)+1)+ncontrols*nodes(iphase)+nparameters;
    t0Index = nstates*(nodes(iphase)+1)+ncontrols*nodes(iphase)+nparameters+1;
    tfIndex = nstates*(nodes(iphase)+1)+ncontrols*nodes(iphase)+nparameters+2;
    t0Scale = colscales_curr(t0Index);
    tfScale = colscales_curr(tfIndex);
    parameterScale = colscales_curr(parameterIndexStart:parameterIndexFinish);
    parameterScales{iphase} = parameterScale;
    t0Scales{iphase} = t0Scale; 
    tfScales{iphase} = tfScale; 
    state_scales_matrix = repmat(state_scales_use,nodes(iphase),1);
    %---------------------------------------%
    % Set bounds and remove infinite values %
    %---------------------------------------%
    tlow = tlowlimits(1);
    tlow(isinf(tlow)) = BIG_NUM * sign(tlow(isinf(tlow)));
    tupp = tupplimits(2);
    tupp(isinf(tupp)) = BIG_NUM * sign(tupp(isinf(tupp)));
    xlow = xlowlimits(:,2).';
    xlow(isinf(xlow)) = BIG_NUM * sign(xlow(isinf(xlow)));
    xupp = xupplimits(:,2).';
    xupp(isinf(xupp)) = BIG_NUM * sign(xupp(isinf(xupp)));
    ulow = ulowlimits.';
    ulow(isinf(ulow)) = BIG_NUM * sign(ulow(isinf(ulow)));
    uupp = uupplimits.';
    uupp(isinf(uupp)) = BIG_NUM * sign(uupp(isinf(uupp)));
    plow = plowlimits.';
    plow(isinf(plow)) = BIG_NUM * sign(plow(isinf(plow)));
    pupp = pupplimits.';
    pupp(isinf(pupp)) = BIG_NUM * sign(pupp(isinf(pupp)));
    %--------------------------------------------------%
    % Generate random sample of variables using bounds %
    %--------------------------------------------------%
    randtime = rand(ntrials,1);
    trand = randtime*tupp+(1-randtime)*tlow;
    if 0,
      if setup.igrid > 1,
        ntrials = length(setup.solution(iphase).time);
        trand = setup.solution(iphase).time;
      end;
    end;
    if nstates>0,
      randstate = rand(nstates,ntrials).';
      xrand = randstate.*repmat(xupp,ntrials,1)+(1-randstate).*repmat(xlow,ntrials,1);
      if 0,
        if setup.igrid > 1,
          xrand = setup.solution(iphase).state;
        end;
      end;
    else
      xrand = [];
    end;
    if ncontrols>0,
      randcontrol = rand(ncontrols,ntrials).';
      urand = randcontrol.*repmat(uupp,ntrials,1)+(1-randcontrol).*repmat(ulow,ntrials,1);
      if 0,
        if setup.igrid > 1,
          urand = setup.solution(iphase).control;
        end;
      end;
    else
      urand = [];
    end;
    if nparameters>0,
      randparameters = rand(nparameters,ntrials).';
      prand = randparameters.*repmat(pupp,ntrials,1)+(1-randparameters).*repmat(plow,ntrials,1);
      if 0,
        if setup.igrid > 1,
          prand = setup.solution(iphase).parameter;
        end;
      end;
    else
      prand = [];
    end;
    %----------------------%
    % Initialize variables %
    %----------------------%
    xuptot = zeros(nstates+ncontrols+nparameters,ntrials);
    dae_norm = zeros(nstates+npaths,ntrials);
    event_norm = zeros(nevents,ntrials);
    extras.phase = iphase;
    extras.nstates = nstates;
    extras.ncontrols = ncontrols;
    extras.nparameters = nparameters;
    for j=1:ntrials;
      %-----------------------------------%
      % Evaluate the Jacobian of the DAEs %
      %-----------------------------------%
      t = trand(j);
      xup = zeros(nstates+ncontrols+nparameters,1);
      if nstates>0,
        xup(1:nstates) = xrand(j,:).';
      end;
      if ncontrols>0,
        xup(nstates+1:nstates+ncontrols) = urand(j,:).';
      end;
      if nparameters>0,
        xup(nstates+ncontrols+1:nstates+ncontrols+nparameters) = prand(j,:).';
      end;
      xuptot(:,j) = xup;
      %----------------------------------------%
      % Evaluate Jacobian using random numbers %
      %----------------------------------------%
      fty = feval('gpopsDaeWrapper',t,xup,extras);
      thresh{1} = threshTol*ones(size(xup));
      thresh{2} = xup;
      [daejac] = numjac('gpopsDaeWrapper',t,xup,fty,thresh,[],0,[],[],extras);
      daejac(isnan(daejac)) = 1;
      dae_norm(:,j) = sqrt(dot(daejac,daejac,2));
      %-------------------------------------%
      % Evaluate the Jacobian of the Events %
      %-------------------------------------%
      if nevents>0,
        xupevent = zeros(2+2*nstates+nparameters,1);
        if nstates>0,
          initevent = [t; xrand(j,:).'];
          termevent = [t; xrand(j,:).'];
          xupevent(1:2+2*nstates) = [initevent; termevent];
        end;
        if nparameters>0,
          xupevent(2+2*nstates+1:end) = prand(j,:).';
        end;
        fty = feval('gpopsEventWrapper',t,xupevent,extras);
        thresh{1} = threshTol*ones(size(xupevent));
        thresh{2} = xupevent;
        [eventjac] = numjac('gpopsEventWrapper',t,xupevent,fty,thresh,[],0,[],[],extras);
        event_norm(:,j) = sqrt(dot(eventjac,eventjac,2));
      end;
    end;
    xupall{iphase} = xuptot;
    if (nstates+npaths)>0,
      dae_norm_average = mean(dae_norm,2);
    end;
    if nstates>0,
      tmaxcurr = setup.limits(iphase).time.max;
      tmincurr = setup.limits(iphase).time.min;
      tscaleCurr = tmaxcurr(2)-tmincurr(1);
      ode_norm_matrix = state_scales_matrix;
    else
      ode_norm_matrix = [];
    end;
    if npaths>0,
      path_norm_average = dae_norm_average(nstates+1:end);
      path_norm_average(logical(path_norm_average<eps)) = 1;
      path_norm_matrix = repmat(path_norm_average,1,nodes(iphase)).';
      pathScales{iphase} = path_norm_matrix(1,:);
    else
      path_norm_matrix = [];
    end;
    if isempty(path_norm_matrix),
      dae_scales{iphase} = state_scales_use;
    else
      dae_scales{iphase} = [state_scales_use path_norm_matrix(1,:)];
    end;
    
    dae_norm_vector = [ode_norm_matrix(:); path_norm_matrix(:)];
    if nevents>0,
      event_norm_vector = mean(event_norm,2);
      event_norm_vector(logical(event_norm_vector<eps)) = 1;
    else
      event_norm_vector = [];
    end;
    rowscales(constraint_indices{1,iphase}) = [dae_norm_vector; event_norm_vector];
    phasecons(iphase) = (nodes(iphase))*nstates+nodes(iphase)*npaths+nevents;
  end;
  %------------------------------%
  % Determine the linkage scales %
  %------------------------------%
  numlinkpairs = setup.numlinkpairs;
  linkages = setup.linkages;
  link_scales = [];
  last_phase_index = sum(phasecons);
  linkage_indices = (last_phase_index+1:numnonlin).';
  t = 0;
  for ipair=1:numlinkpairs;
    left_phase               = linkages(ipair).left.phase;
    right_phase               = linkages(ipair).right.phase;
    nstates_left             = sizes(left_phase,1);
    nstates_right            = sizes(right_phase,1);
    ncontrols_left           = sizes(left_phase,2);
    ncontrols_right          = sizes(right_phase,2);
    nparameters_left         = sizes(left_phase,3);
    nparameters_right        = sizes(right_phase,3);
    nlinks                   = length(linkages(ipair).min);
    extras.left.phase        = left_phase;
    extras.left.nstates      = nstates_left;
    extras.left.nparameters  = nparameters_left;
    extras.right.phase       = right_phase;
    extras.right.nstates     = nstates_right;
    extras.right.nparameters = nparameters_right;
    link_norm                = zeros(nlinks,ntrials);
    for j=1:ntrials;
      xupleft  = xupall{left_phase}(:,ipair);
      xupright = xupall{right_phase}(:,ipair);
      xleft = xupleft(1:nstates_left);
      pleft = xupleft(nstates_left+ncontrols_left+1:end);
      xright = xupright(1:nstates_right);
      pright = xupright(nstates_right+ncontrols_right+1:end);
      xplink = [xleft; pleft; xright; pright];
      fty = feval('gpopsLinkWrapper',t,xplink,extras);
      thresh{1} = threshTol*ones(size(xplink));
      thresh{2} = xplink;
      [linkjac] = numjac('gpopsLinkWrapper',t,xplink,fty,thresh,[],0,[],[],extras);
      link_norm(:,j) = sqrt(dot(linkjac,linkjac,2));
    end;
    link_norm_average = mean(link_norm,2);
    link_norm_average(logical(link_norm_average<eps)) = 1;
    link_scales = [link_scales; link_norm_average];
  end;
  rowscales(linkage_indices) = link_scales;
else  % END:  isfield(setup,'autoscale') && isequal(setup.autoscale,'on')
  for iphase=1:numphases;
    nstates = sizes(iphase,1);
    ncontrols = sizes(iphase,2);
    nparameters = sizes(iphase,3);
    npaths = sizes(iphase,4);
    nevents = sizes(iphase,5);
    stateScales{iphase} = ones(1,nstates);
    controlScales{iphase} = ones(1,ncontrols);
    parameterScales{iphase} = ones(1,nparameters);
    odeScales{iphase} = ones(1,nstates);
    pathScales{iphase} = ones(1,npaths);
  end
end

setup.column_scales = 1./colscales;
setup.column_shifts = colshifts;
setup.row_scales = 1./rowscales;

if exist('dae_norm','var')
  setup.autoscaleinfo.dae_norm = dae_norm;
  setup.autoscaleinfo.dae_scales = dae_scales;
  setup.autoscaleinfo.xuptot = xuptot;
end
if exist('event_norm','var')
  setup.autoscaleinfo.event_norm = event_norm;
end
if exist('link_norm','var')
  setup.autoscaleinfo.link_norm = link_norm;
end
%----------------------------------%
% Find worst-case sparsity pattern %
%----------------------------------%
for iphase = 1:numphases
  % dependencies = (nstates + npaths) x (nstates + ncontrols)
  dependencies{iphase} = ones(sizes(iphase,1)+sizes(iphase,4),...
                              sizes(iphase,1)+sizes(iphase,2));
end
setup.dependencies = dependencies;
%---------------------------------%
% Check for user sparsity pattern %
%---------------------------------%
useUserSparsity = 0;
for iphase=1:numphases;
  if isfield(setup.limits(iphase),'dependencies') && ~isempty(setup.limits(iphase).dependencies)
    useUserSparsity = 1;
    %---------------------------%
    % Override Default Sparsity %
    %---------------------------%
    setup.dependencies{iphase} = setup.limits(iphase).dependencies;
    %------------%
    % Check Size %
    %------------%
    %-------------------------%
    % drow = nstates + npaths %
    %-------------------------%
    drow = sizes(iphase,1) + sizes(iphase,4); 
    %----------------------------%
    % dcol = nstates + ncontrols %
    %----------------------------%
    dcol = sizes(iphase,1) + sizes(iphase,2); 
    [drow_in,dcol_in] = size(setup.dependencies{iphase});
    if ~isequal(drow,drow_in) || ~isequal(dcol,dcol_in)
      error(['User defined sparsity in phase %i has incorrect size:',...
             '\n\tUser: (%i X %i) Needs: (%i X %i)\n'],...
            iphase,drow_in,dcol_in,drow,dcol)
    elseif ~all(setup.limits(iphase).dependencies(:) == 0 | setup.limits(iphase).dependencies(:) == 1)
      error(['User defined sparsity in phase %i has incorrect values',...
             '\n\tAll values must be 0 or 1\n'],iphase)
    end
    %--------------------------------%
    % Check for Missing Dependencies %
    %--------------------------------%
    gpopsCheckDependencies(setup,iphase);
  end
end
if setup.printoff == 0
  fprintf('\n');
  fprintf(' _____________________________________________________________________________\n');
  fprintf('|                                                                             |\n');
  if useUserSparsity == 1
    fprintf('| Using User Defined Sparsity                                                 |\n');
  else
    fprintf('| Using Default Sparsity                                                      |\n');
  end
end
setup.stateScales  = stateScales;     
setup.controlScales  = controlScales;     
setup.odeScales  = odeScales;     
setup.pathScales = pathScales;
setup.t0Scales = t0Scales;
setup.tfScales = tfScales;
setup.parameterScales = parameterScales;

%---------------------------------------------------------------------------------------------%
% BEGIN:  LEGACY CODE FOR COMPUTING COLUMN SCALES                                             %
% minmaxnotzero = find((xmin~=0) & (xmax~=0) & ~isinf(xmin) & ~isinf(xmax));                  %
% colscales(minmaxnotzero)          = max(abs(xmin(minmaxnotzero)),abs(xmax(minmaxnotzero))); %
% colshifts(minmaxnotzero) = -xmax(minmaxnotzero)./colscales(minmaxnotzero);                  %
% minzero                           = find((xmin==0) & (xmax~=0) & ~isinf(xmax));             %
% maxzero                           = find((xmin~=0) & ~isinf(xmin) & (xmax==0));             %
% colscales(minzero)            = abs(xmax(minzero));                                         %
% colscales(maxzero)            = abs(xmin(maxzero));                                         %
% mininfmaxnotzero                = find(isinf(xmin) & ~isinf(xmax) & (xmax~=0));             %
% colscales(mininfmaxnotzero) = abs(xmax(mininfmaxnotzero));                                  %
% maxinfminnotzero                = find(isinf(xmax) & ~isinf(xmin) & (xmin~=0));             %
% colscales(maxinfminnotzero) = abs(xmin(maxinfminnotzero));                                  %
% END:  LEGACY CODE FOR COMPUTING COLUMN SCALES                                               %
%---------------------------------------------------------------------------------------------%
