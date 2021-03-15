function setup = gpopsGetBounds(setup)
%------------------------------------------------------------------%
% Get all bounds in a multiple-phase optimal control problem       %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

limits = setup.limits;
linkages = setup.linkages;
sizes = setup.sizes;
numphases = setup.numphases;
nodes = horzcat(limits(:).nodes);
variable_offset = 0;
constraint_offset = 0;
nlplimits = cell(numphases,4);
variables = zeros(1,numphases);
constraints = zeros(1,numphases);
variable_indices = cell(1,numphases);
constraint_indices = cell(1,numphases);
indices = repmat(struct('state',[],'control',[],'time',[],'parameter',[]),numphases,1);
for iphase=1:numphases;
    % ---------------------------------------------------------
    % Get the lower and upper limits on the initial and terminal
    % time in the current phase
    % ---------------------------------------------------------
    tMin = limits(iphase).time.min;
    tMax = limits(iphase).time.max;
    t0_min = tMin(1);
    t0_max = tMax(1);
    tf_min = tMin(2);
    tf_max = tMax(2);
    if size(tMin(:),1) ~= 2 || size(tMax(:),1) ~= 2
        error('Time Bound Vectors Must Have Length 2, in Phase: %i',iphase)
    end
    if any(isnan(tMin)) || any(isnan(tMax))
        error('Bounds on Time are NaN in Phase: %i',iphase)
    end
    % -------------------------------------------------------------------
    % Get the lower and upper limits on the states in the current phase
    % -------------------------------------------------------------------
    state_matrix_min = zeros(nodes(iphase)+1,sizes(iphase,1));
    state_matrix_max = zeros(nodes(iphase)+1,sizes(iphase,1));
    if ~isequal(sizes(iphase,1),0),
        state0Min = limits(iphase).state.min(:,1);
        stateMin = limits(iphase).state.min(:,2);
        statefMin = limits(iphase).state.min(:,3);
        state0Max = limits(iphase).state.max(:,1);
        stateMax = limits(iphase).state.max(:,2);
        statefMax = limits(iphase).state.max(:,3);
        state_matrix_min(1,:) = state0Min;
        state_matrix_min(2:nodes(iphase),:) = repmat(stateMin,1,nodes(iphase)-1).';
        state_matrix_min(nodes(iphase)+1,:) = statefMin;
        state_matrix_max(1,:) = state0Max;
        state_matrix_max(2:nodes(iphase),:) = repmat(stateMax,1,nodes(iphase)-1).';
        state_matrix_max(nodes(iphase)+1,:) = statefMax;
        % -------------------------------------------------------------------
        % Check the lower and upper limits on the states in the current phase
        % -------------------------------------------------------------------
        if (any(state0Max - stateMin < 0) || any(statefMax - stateMin < 0) || any(stateMax - state0Min < 0) || any(stateMax - statefMin < 0) || any(stateMax(:,:) - stateMin(:,:) < 0))
            error('Bounds on State are Inconsistent (i.e. max < min) in Phase: %i',iphase)
        end
        if any(any(isnan(stateMin))) || any(any(isnan(stateMax)))
            error('Bounds on State are NaN in Phase: %i',iphase)
        end
    else
        state_matrix_min = [];
        state_matrix_max = [];
    end;
    % -------------------------------------------------------------------
    % Get the lower and upper limits on the controls in the current phase
    % -------------------------------------------------------------------
    if ~isequal(sizes(iphase,2),0),
        controlMin = limits(iphase).control.min;
        controlMax = limits(iphase).control.max;
        control_matrix_min = repmat(controlMin,1,nodes(iphase)).';
        control_matrix_max = repmat(controlMax,1,nodes(iphase)).';
        % -------------------------------------------------------------------
        % Check the lower and upper limits on the controls in the current phase
        % -------------------------------------------------------------------
        if any(controlMax(:,:) - controlMin(:,:) < 0)
            error('Bounds on Control are Inconsistent (i.e. max < min) in Phase: %i',iphase)
        end
        if any(isnan(controlMin)) || any(isnan(controlMax))
            error('Bounds on Control are NaN in Phase: %i',iphase)
        end
    else
        control_matrix_min = [];
        control_matrix_max = [];
    end;
    % ----------------------------------------------------------------------------
    % Get the lower and upper limits on the static parameters in the current phase
    % ----------------------------------------------------------------------------
    if ~isequal(sizes(iphase,3),0),
        parameter_min = limits(iphase).parameter.min;
        parameter_max = limits(iphase).parameter.max;
        % -------------------------------------------------------------------
        % Check the lower and upper limits on the static parameters in the current phase
        % -------------------------------------------------------------------
        if any(parameter_max(:,:) - parameter_min(:,:) < 0)
            error('Bounds on Parameters are Inconsistent (i.e. max < min) in Phase: %i',iphase)
        end
        if any(isnan(parameter_min)) || any(isnan(parameter_max))
            error('Bounds on Parameters are NaN in Phase: %i',iphase)
        end
    else
        parameter_min = [];
        parameter_max = [];
    end;
    % ---------------------------------------------------------------------------
    % Get the lower and upper limits on the path constraints in the current phase
    % ---------------------------------------------------------------------------
    if ~isequal(sizes(iphase,4),0),
        pathMin = limits(iphase).path.min;
        pathMax = limits(iphase).path.max;
        path_matrix_min = repmat(pathMin,1,nodes(iphase)).';
        path_matrix_max = repmat(pathMax,1,nodes(iphase)).';
        % -------------------------------------------------------------------
        % Check the lower and upper limits on the path constraints in the current phase
        % -------------------------------------------------------------------
        if any(pathMax(:,:) - pathMin(:,:) < 0)
            error('Bounds on Path Constraints are Inconsistent (i.e. max < min) in Phase: %i',iphase)
        end
        if any(isnan(pathMin)) || any(isnan(pathMax))
            error('Bounds on Path Constraints are NaN in Phase: %i',iphase)
        end
    else
        path_matrix_min = [];
        path_matrix_max = [];
    end;
    % ----------------------------------------------------------------------------
    % Get the lower and upper limits on the event constraints in the current phase
    % ----------------------------------------------------------------------------
    if ~isequal(sizes(iphase,5),0),
        event_vector_min = limits(iphase).event.min;
        event_vector_max = limits(iphase).event.max;
        % -------------------------------------------------------------------
        % Check the lower and upper limits on the event constraints in the current phase
        % -------------------------------------------------------------------
        if any(event_vector_max(:,:) - event_vector_min(:,:) < 0)
            error('Bounds on Event Constraints are Iconsistent (i.e. max < min) in Phase: %i',iphase)
        end
        if any(isnan(event_vector_min)) || any(isnan(event_vector_max))
            error('Bounds on Event Constraints are NaN in Phase: %i',iphase)
        end
    else
        event_vector_min = [];
        event_vector_max = [];
    end;
    state_vector_min = state_matrix_min(:);
    state_vector_max = state_matrix_max(:);
    control_vector_min = control_matrix_min(:);
    control_vector_max = control_matrix_max(:);
    ode_vector_min = zeros((nodes(iphase))*sizes(iphase,1),1);
    ode_vector_max = zeros((nodes(iphase))*sizes(iphase,1),1);
    path_vector_min = path_matrix_min(:);
    path_vector_max = path_matrix_max(:);
    % ------------------------------------------------------------------------------------------------------
    % The cell array NLPLIMITS contains the lower and upper limits on the NLP variables in the current phase
    % ------------------------------------------------------------------------------------------------------
    nlplimits{iphase,1} = [state_vector_min; control_vector_min; t0_min; tf_min; parameter_min];
    nlplimits{iphase,2} = [state_vector_max; control_vector_max; t0_max; tf_max; parameter_max];
    nlplimits{iphase,3} = [ode_vector_min; path_vector_min; event_vector_min];
    nlplimits{iphase,4} = [ode_vector_max; path_vector_max; event_vector_max];
    variables(iphase) = length(nlplimits{iphase,1});
    constraints(iphase) = length(nlplimits{iphase,3});
    variable_indices{iphase} = variable_offset+1:variable_offset+variables(iphase);
    constraint_indices{iphase} = constraint_offset+1:constraint_offset+constraints(iphase);
    nstates = sizes(iphase,1);
    ncontrols = sizes(iphase,2);
    nparameters = sizes(iphase,3);
    state_indices = variable_offset+1:variable_offset+(nodes(iphase)+1)*nstates;
    control_indices = state_indices(end)+1:state_indices(end)+nodes(iphase)*ncontrols;
    if ~isempty(control_indices),
        t0_index = control_indices(end)+1;
    else
        t0_index = state_indices(end)+1;
    end;
    tf_index = t0_index+1;
    parameter_indices = tf_index+1:tf_index+nparameters;
    indices(iphase).state = state_indices;
    indices(iphase).control = control_indices;
    indices(iphase).time = [t0_index; tf_index];
    indices(iphase).parameter = parameter_indices;
    variable_offset = variable_offset+variables(iphase);
    constraint_offset = constraint_offset+constraints(iphase);
end;
numlinkpairs = setup.numlinkpairs;
linkMinTot = zeros(setup.numlinks,1);
linkMaxTot = zeros(setup.numlinks,1);
linkRowShift = 1;
for ipair=1:numlinkpairs;
    % ------------------------------------------------------------%
    % Check the lower and upper limits in the linkage constraints %
    % ------------------------------------------------------------%
    linkMin = linkages(ipair).min;
    linkMax = linkages(ipair).max;
    if any(linkMin - linkMax > 0)
        error('Bounds on Linkage Constraints are Inconsistent (i.e. max < min) in Linkage Pair %i',ipair)
    end
    if any(isnan(linkMin)) || any(isnan(linkMax))
        error('Bounds on Linkage Constraints are NaN in Linkage Pair %i',ipair)
    end
    linkRows = linkRowShift:linkRowShift-1+length(linkMin);
    linkMinTot(linkRows) = linkMin;
    linkMaxTot(linkRows) = linkMax;
    linkRowShift = linkRowShift+length(linkMin);
end;
varbounds_min = vertcat(nlplimits{:,1});
varbounds_max = vertcat(nlplimits{:,2});
conbounds_min = vertcat(nlplimits{:,3});
conbounds_max = vertcat(nlplimits{:,4});
conbounds_min = [conbounds_min; linkMinTot];
conbounds_max = [conbounds_max; linkMaxTot];
setup.varbounds_min = varbounds_min;
setup.varbounds_max = varbounds_max;
setup.conbounds_min = conbounds_min;
setup.conbounds_max = conbounds_max;
setup.variables = variables;
setup.constraints = constraints;
setup.variable_indices = variable_indices;
setup.constraint_indices = constraint_indices;
setup.numphases = numphases;
setup.numnonlin = length(conbounds_min);
setup.indices = indices;

% --------------------------------
% Get Bounds on Linear Constraints
% --------------------------------

numvars = length(varbounds_min);
Alinear_I = zeros(2*(numphases+numlinkpairs),1);
Alinear_J = zeros(2*(numphases+numlinkpairs),1);
Alinear_V = zeros(2*(numphases+numlinkpairs),1);
AlinRowshift = 0;
Alinmin = zeros(numphases+numlinkpairs,1);
Alinmax = zeros(numphases+numlinkpairs,1);
% ---------------------------------------------
% Part 1:  Monotonicity of Independent Variable
% ---------------------------------------------
for iphase=1:numphases;
    nstates = sizes(iphase,1);
    ncontrols = sizes(iphase,2);
    if ~isequal(iphase,1),
        ishift = variable_indices{iphase-1}(end);
    else
        ishift =0;
    end;
    t0_index = ishift+(nodes(iphase)+1)*nstates+nodes(iphase)*ncontrols+1;
    tf_index = t0_index+1;
    AlinRowshift = AlinRowshift + 1;
    Alinear_I(AlinRowshift) = iphase;
    Alinear_J(AlinRowshift) = t0_index;
    Alinear_V(AlinRowshift) = -1;
    AlinRowshift = AlinRowshift + 1;
    Alinear_I(AlinRowshift) = iphase;
    Alinear_J(AlinRowshift) = tf_index;
    Alinear_V(AlinRowshift) = 1;
    if isfield(limits(iphase),'duration'),
        % Check if only ONE of the lower and upper bounds on the phase duration are specified
        if (isfield(limits(iphase).duration,'min') && ~ isfield(limits(iphase).duration,'max')) || (~isfield(limits(iphase).duration,'min') && isfield(limits(iphase).duration,'max')),
            error('Must Specify Both Lower and Upper Bounds on Phase Duration');
        else
            if isempty(limits(iphase).duration.min) && isempty(limits(iphase).duration.max),
                % Neither Minimum nor Maximum Duration Specified
                Alinmin(iphase) = 0;
                Alinmax(iphase) = Inf;
            elseif ~isempty(limits(iphase).duration.min) && isempty(limits(iphase).duration.max),
                % Minimum Duration Specified But Maximum Duration Unspecified
                if isequal(size(limits(iphase).duration.min),[1 1]),
                    Alinmin(iphase) = limits(iphase).duration.min;
                    Alinmax(iphase) = Inf;
                end;
            elseif isempty(limits(iphase).duration.min) && ~isempty(limits(iphase).duration.max),
                % Minimum Duration Unspecified But Maximum Duration Specified
                if isequal(size(limits(iphase).duration.max),[1 1]),
                    Alinmin(iphase) = 0;
                    Alinmin(iphase) = limits(iphase).duration.max;
                end;
            else
                % Both Minimum and Maximum Duration Specified
                if ~(isequal(size(limits(iphase).duration.min),[1 1]) && isequal(size(limits(iphase).duration.max),[1 1])),
                    errStr1 = 'Lower and Upper Bounds on Phase Duration in Phase ';
                    errStr2 = ' Must be Scalars';
                    error('%s, in Linkage Pair %i %s',errStr1,ipair,errStr2)
                else
                    Alinmin(iphase) = limits(iphase).duration.min;
                    Alinmax(iphase) = limits(iphase).duration.max;
                end;
            end;
        end;
    else
        Alinmin(iphase) = 0;
        Alinmax(iphase) = Inf;
    end;
end;
istart = numphases;
% --------------------------------------
% Part 2:  Linkage of Time Across Phases
% --------------------------------------
for ipair=1:numlinkpairs;
    left_phase = linkages(ipair).left.phase;
    right_phase = linkages(ipair).right.phase;
    nparameters_left = sizes(left_phase,3);
    nparameters_right = sizes(right_phase,3);
    tf_index_left = variable_indices{left_phase}(end-nparameters_left);
    t0_index_right = variable_indices{right_phase}(end-nparameters_right-1);
    AlinRowshift = AlinRowshift + 1;
    Alinear_I(AlinRowshift) = istart+ipair;
    Alinear_J(AlinRowshift) = tf_index_left;
    Alinear_V(AlinRowshift) = -1;
    AlinRowshift = AlinRowshift + 1;
    Alinear_I(AlinRowshift) = istart+ipair;
    Alinear_J(AlinRowshift) = t0_index_right;
    Alinear_V(AlinRowshift) = 1;
    Alinmin(istart+ipair) = 0;
    Alinmax(istart+ipair) = 0;
end;
setup.numvars   = numvars;
setup.numlin    = numphases+numlinkpairs;
setup.Alinear   = sparse(Alinear_I,Alinear_J,Alinear_V,numphases+numlinkpairs,numvars);
setup.Alinmin   = Alinmin;
setup.Alinmax   = Alinmax;
