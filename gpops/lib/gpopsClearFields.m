function setup = gpopsClearFields(setup)
%------------------------------------------------------------------%
% Clear fields from the structure setup at end of a NLP solver run %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

fields_to_remove = {
                    'sizes'
                    'varbounds_min'
                    'varbounds_max'
                    'conbounds_min'
                    'conbounds_max'
                    'nlplimits'
                    'variables'
                    'constraints'
                    'variable_indices'
                    'constraint_indices'
                    'numphases'
                    'numnonlin'
                    'numlinkpairs'
                    'numlinks'
                    'indices'
                    'numconnections'
                    'numvars'
                    'numlin'
                    'Alinear'
                    'Alinmin'
                    'Alinmax'
                    'nlpGuess'
                    'ps'
                    'column_scales'
                    'row_scales'
                    'dependencies'
                    'sparsity_nonlinear'
                    'sparsity_jac'
                    'result'
                    'sparsity_constant'
                    'sparsity_all'
                    'Dx'
                    'invDx'
                    'DF'
                    'invDF'
                    'iGfun'
                    'jGvar'
                    'initlincons'
                    'seed'
                    'color_groups'
                    'hpert'
                    'deltaxmat'
                    'Jaczeros'
                    'autoscaleinfo'
                    'stateScales'
                    'controlScales'
                    'odeScales'
                    'pathScales'
                    't0Scales'
                    'tfScales'
                    'parameterScales'
                    };
for i = 1:length(fields_to_remove) 
    if isfield(setup,fields_to_remove{i})
        setup = rmfield(setup,fields_to_remove{i});
    end
end

