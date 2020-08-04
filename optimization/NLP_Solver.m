%
%   Copyright 2020 Janis Maczijewski
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
%

classdef NLP_Solver
    
    properties
        variables_pack
        variables_unpack
        
        parameters_pack
        parameters_unpack
        
        slack_init_lower_bound_fn
        slacks_unpack
        g_hard_fn
        g_hard_unpack_fn
        
        lbg
        ubg
        
        objective_parts_fn
        objective_terms
        
        nlpsolver
        nlp_struct
        nlp_fn
        nlp_options
    end
    
    methods
        function obj = NLP_Solver(nlp, sweep_config)
            variables_flat = struct2cell(map_struct(nlp.variables, @(e)e(:)));
            variables_flat = vertcat(variables_flat{:});
            
            obj.variables_pack = casadi.Function('variables_pack', ...
                struct2cell(nlp.variables), {variables_flat},...
                fieldnames(nlp.variables), {'variables_flat'});
            
            obj.variables_unpack = casadi.Function('variables_unpack', ...
                {variables_flat}, struct2cell(nlp.variables),...
                {'variables_flat'}, fieldnames(nlp.variables));
            
            parameters_flat = struct2cell(map_struct(nlp.parameters, @(e)e(:)));
            parameters_flat = vertcat(parameters_flat{:});
            parameter_init_variables_flat = casadi.SX.sym('parameter_init_variables_flat',size(variables_flat,1), size(variables_flat,2));
            parameter_init_variables_flat_weight = casadi.SX.sym('parameter_init_variables_flat_weight',1,1);
            
            obj.parameters_pack = casadi.Function('parameters_pack', ...
                struct2cell(nlp.parameters), {parameters_flat},...
                fieldnames(nlp.parameters), {'parameters_flat'});
            
            obj.parameters_unpack = casadi.Function('parameters_unpack', ...
                {parameters_flat}, struct2cell(nlp.parameters),...
                {'parameters_flat'}, fieldnames(nlp.parameters));
            
            
            
            equations = casadi.SX([]); % each entry == 0
            g_hard = casadi.SX([]); % each entry <= 0
            g_soft = casadi.SX([]); % each entry <= slack_variable
            g_hard_placeholders_flat = casadi.SX([]);
            g_hard_placeholders = struct;
            slack_variables_flat = casadi.SX([]);
            slack_variables = struct;
            slack_objective = casadi.SX(0);
            slack_init_lower_bounds = [];
            
            
            
            constraint_names = fieldnames(nlp.constraints);
            for i = 1:numel(constraint_names)
                constraint = nlp.constraints.(constraint_names{i});
                if constraint.is_hard_constraint
                    if constraint.is_equation
                        equations = [equations; constraint.g(:)];
                    else
                        g_hard = [g_hard; (constraint.g(:) - constraint.max_value(:))];
                        g_hard_ub_placeholder = casadi.SX.sym([constraint_names{i} '_ub'], size(constraint.g, 1), size(constraint.g, 2)); % These are used to unpack the flat g_hard vector, and inspect which constraints are active.
                        g_hard_placeholders_flat = [g_hard_placeholders_flat; g_hard_ub_placeholder(:)];
                        g_hard_placeholders.([constraint_names{i} '_ub']) = g_hard_ub_placeholder;
                        
                        if constraint.has_lower_bound
                            g_hard = [g_hard; (constraint.min_value(:) - constraint.g(:))];
                            g_hard_lb_placeholder = casadi.SX.sym([constraint_names{i} '_lb'], size(constraint.g, 1), size(constraint.g, 2));
                            g_hard_placeholders_flat = [g_hard_placeholders_flat; g_hard_lb_placeholder(:)];
                            g_hard_placeholders.([constraint_names{i} '_lb']) = g_hard_lb_placeholder;
                        end
                    end
                else
                    % Handle soft constraints (exact penalty method, with slack variables)
                    slack_var = casadi.SX.sym(constraint_names{i}, size(constraint.g, 1), size(constraint.g, 2));
                    slack_variables.(constraint_names{i}) = slack_var;
                    slack_var_flat = slack_var(:);
                    weighted_slack = constraint.weight .* slack_var;
                    slack_objective = slack_objective + sum(weighted_slack(:));
                    slack_variables_flat = [slack_variables_flat; slack_var_flat];
                    
                    
                    g_soft_entry_ub = constraint.g - constraint.max_value;
                    g_soft_entry_ub_flat = g_soft_entry_ub(:);
                    g_soft = [g_soft; (g_soft_entry_ub_flat-slack_var_flat)];
                    
                    slack_init_lower_bound_entry = casadi.SX.zeros(size(g_soft_entry_ub_flat, 1), 3);
                    slack_init_lower_bound_entry(:,1) = g_soft_entry_ub_flat;
                    
                    if constraint.has_lower_bound
                        g_soft_entry_lb = constraint.min_value - constraint.g;
                        g_soft_entry_lb_flat = g_soft_entry_lb(:);
                        g_soft = [g_soft; (g_soft_entry_lb_flat-slack_var_flat)];
                        
                        slack_init_lower_bound_entry(:,2) = g_soft_entry_lb_flat;
                    end
                    
                    slack_init_lower_bounds = [slack_init_lower_bounds; slack_init_lower_bound_entry];
                end
            end
            
            obj.slack_init_lower_bound_fn = casadi.Function('slack_init_lower_bound_fn', {variables_flat, parameters_flat}, {slack_init_lower_bounds});
            obj.g_hard_fn = casadi.Function('g_hard_fn', {variables_flat, parameters_flat}, {g_hard});
            
            obj.slacks_unpack = casadi.Function('slacks_unpack', ...
                {slack_variables_flat}, struct2cell(slack_variables), ...
                {'slacks_flat'}, fieldnames(slack_variables));
            
            obj.g_hard_unpack_fn = casadi.Function('g_hard_unpack_fn', ...
                {g_hard_placeholders_flat}, struct2cell(g_hard_placeholders), ...
                {'g_hard'}, fieldnames(g_hard_placeholders));
            
            g = [equations; g_hard; g_soft];
            obj.lbg = [zeros(size(equations));   -inf(size(g_hard));   -inf(size(g_soft))];
            obj.ubg = [zeros(size(equations));  zeros(size(g_hard));  zeros(size(g_soft))];
            
            objective_soft_trust_region = parameter_init_variables_flat_weight * sum((variables_flat - parameter_init_variables_flat).^2);
            objective = nlp.objective + slack_objective + objective_soft_trust_region;
            
            obj.objective_terms = nlp.objective_terms;
            obj.objective_terms.solver_nlp = nlp.objective;
            obj.objective_terms.solver_slack = slack_objective;
            obj.objective_terms.solver_soft_trust_region = objective_soft_trust_region;
            obj.objective_terms.solver_total = objective;
            
            objective_term_names = fieldnames(obj.objective_terms);
            for i = 1:length(objective_term_names)
                objective_terms_SX{i} = obj.objective_terms.(objective_term_names{i});
            end
            
            obj.objective_parts_fn = casadi.Function('objective_parts_fn', ...
                {[variables_flat; slack_variables_flat], parameters_flat, parameter_init_variables_flat, parameter_init_variables_flat_weight}, ...
                objective_terms_SX, ...
                {'x', 'parameters_flat', 'parameter_init_variables_flat', 'parameter_init_variables_flat_weight'}, ...
                objective_term_names);
            
            obj.nlp_struct = struct(...
                'x', [variables_flat; slack_variables_flat], ...
                'f', objective, ...
                'g', g, ...
                'p', [parameters_flat;parameter_init_variables_flat;parameter_init_variables_flat_weight] ...
            );
        
            obj.nlp_fn = casadi.Function('nlp_fn', ...
                {obj.nlp_struct.x, obj.nlp_struct.p},...
                {obj.nlp_struct.f, obj.nlp_struct.g},...
                {'x','p'},{'f','g'});
        
            % https://www.coin-or.org/Bonmin/option_pages/options_list_ipopt.html#sec:Convergence
            options = struct;
            options.ipopt.max_iter = 400;
            options.ipopt.tol = sweep_config.ipopt_tol;
            options.ipopt.max_cpu_time = 1800;
            options.ipopt.mu_strategy = 'adaptive';
            
            options.ipopt.nlp_scaling_method = sweep_config.ipopt_nlp_scaling_method;
            
            
            options.ipopt.bound_frac = sweep_config.ipopt_bound_push;
            options.ipopt.bound_push = sweep_config.ipopt_bound_push;
            options.ipopt.slack_bound_push = sweep_config.ipopt_bound_push;
            
            options.ipopt.warm_start_init_point = 'yes';
            options.ipopt.warm_start_mult_bound_push = 1e-6;
            options.ipopt.warm_start_bound_frac = 1e-6;
            options.ipopt.warm_start_bound_push = 1e-6;
            options.ipopt.warm_start_slack_bound_push = 1e-6;
            options.ipopt.warm_start_slack_bound_frac = 1e-6;
            
            
            options.ipopt.mu_max = sweep_config.ipopt_mu_max;
                
            switch sweep_config.solver
                case 'ipopt_ma57'
                    options.ipopt.linear_solver = 'ma57';
                    obj.nlpsolver = casadi.nlpsol('nlpsolver', 'ipopt', obj.nlp_struct, options);
                    
                case 'ipopt_mumps'
                    options.ipopt.linear_solver = 'mumps';
                    options.ipopt.mumps_permuting_scaling = 7;
                    options.ipopt.mumps_scaling = 8;
                    obj.nlpsolver = casadi.nlpsol('nlpsolver', 'ipopt', obj.nlp_struct, options);
                    
                case 'ipopt_ma27'
                    options.ipopt.linear_solver = 'ma27';
                    obj.nlpsolver = casadi.nlpsol('nlpsolver', 'ipopt', obj.nlp_struct, options);
                    
                case 'worhp'
                    obj.nlpsolver = casadi.nlpsol('nlpsolver', 'worhp', obj.nlp_struct);
                    
                otherwise
                    error(['Invalid solver name "' solver_name '"']);
            end
            
            obj.nlp_options = options;
        end
        
        function parameters_struct = get_parameters_struct(obj, init_value)
            parameters_struct = obj.parameters_unpack.call(struct('parameters_flat', init_value));
            parameters_struct = map_struct(parameters_struct, @(e)full(e));
        end
        
        function variables_struct = get_variables_struct(obj, init_value)
            variables_struct = obj.variables_unpack.call(struct('variables_flat', init_value));
            variables_struct = map_struct(variables_struct, @(e)full(e));
        end
        
        function results = solve(obj, initialization, lower_bounds, upper_bounds, parameters, initialization_weight, previous_result)
            
            initialization_flat = obj.variables_pack.call(initialization);
            initialization_flat = initialization_flat.variables_flat;
            
            lower_bounds_flat = obj.variables_pack.call(lower_bounds);
            lower_bounds_flat = lower_bounds_flat.variables_flat;
            
            upper_bounds_flat = obj.variables_pack.call(upper_bounds);
            upper_bounds_flat = upper_bounds_flat.variables_flat;
            
            parameters_flat_values = obj.parameters_pack.call(parameters);
            parameters_flat_values = parameters_flat_values.parameters_flat;
            
            assert(~any(isnan(full(initialization_flat))));
            assert(~any(isnan(full(lower_bounds_flat))));
            assert(~any(isnan(full(upper_bounds_flat))));
            assert(~any(isnan(full(parameters_flat_values))));
            assert(all(full(initialization_flat) < 1e9));
            assert(all(full(initialization_flat) > -1e9));
            assert(all(full(initialization_flat) < full(upper_bounds_flat) + 1e-4));
            assert(all(full(initialization_flat) > full(lower_bounds_flat) - 1e-4));
            
            slack_init_lower_bound = obj.slack_init_lower_bound_fn.call({initialization_flat, parameters_flat_values});
            slack_init_lower_bound = slack_init_lower_bound{1};
            initialization_slacks = max(full(slack_init_lower_bound),[],2) + 1e-5;
            
            nlp_fn_values = obj.nlp_fn.call({[initialization_flat; initialization_slacks], [parameters_flat_values;initialization_flat;initialization_weight]});
            g_value = nlp_fn_values{2};
            
            max_infeasibility_initialization = max(max(full(g_value) - obj.ubg), max(obj.lbg - full(g_value)));
            fprintf('Initialization infeasibility %e\n', max_infeasibility_initialization);
            
            nlp_inputs = struct;
            nlp_inputs.x0 = [initialization_flat; initialization_slacks];
            nlp_inputs.lbg = obj.lbg;
            nlp_inputs.ubg = obj.ubg;
            nlp_inputs.lbx = [lower_bounds_flat; zeros(size(initialization_slacks))];
            nlp_inputs.ubx = [upper_bounds_flat; inf * ones(size(initialization_slacks))];
            nlp_inputs.p = [parameters_flat_values;initialization_flat;initialization_weight];
            
            
            %% Run solver, warm-start if possible
            if     isstruct(previous_result) ...
                && isfield(previous_result, 'raw_solution') ...
                && numel(previous_result.raw_solution.lam_g) == numel(nlp_inputs.lbg) ...
                && numel(previous_result.raw_solution.lam_x) == numel(nlp_inputs.x0) ...
                
                nlp_result = obj.nlpsolver(...
                    'x0', nlp_inputs.x0, ...
                    'lbg', nlp_inputs.lbg, ...
                    'ubg', nlp_inputs.ubg, ...
                    'lbx', nlp_inputs.lbx, ...
                    'ubx', nlp_inputs.ubx, ...
                    'p', nlp_inputs.p, ...
                    'lam_x0', previous_result.raw_solution.lam_x, ...
                    'lam_g0', previous_result.raw_solution.lam_g ...
                );
            else
                nlp_result = obj.nlpsolver(...
                    'x0', nlp_inputs.x0, ...
                    'lbg', nlp_inputs.lbg, ...
                    'ubg', nlp_inputs.ubg, ...
                    'lbx', nlp_inputs.lbx, ...
                    'ubx', nlp_inputs.ubx, ...
                    'p', nlp_inputs.p ...
                );
            end

            objective_term_values = obj.objective_parts_fn.call(...
                struct(...
                    'x', nlp_result.x, ...
                    'parameters_flat', parameters_flat_values, ...
                    'parameter_init_variables_flat', initialization_flat, ...
                    'parameter_init_variables_flat_weight', initialization_weight ...
                ) ...
            );
            objective_term_values = map_struct(objective_term_values, @(e)full(e));
            
            % Check constraint solution
            nlp_fn_values = obj.nlp_fn.call({nlp_result.x, [parameters_flat_values;initialization_flat;initialization_weight]});
            g_value = nlp_fn_values{2};
            max_infeasibility_solution = max(max(full(g_value) - obj.ubg), max(obj.lbg - full(g_value)));
            assert(max_infeasibility_solution < 1e-4, 'Infeasible result from CasADi/nlpsol()');
            
            % Unpack results
            primal_soln = full(nlp_result.x);       
            solution_flat = primal_soln(1:length(initialization_flat));
            slacks_flat = full(primal_soln((1+length(initialization_flat)):end));
            
            assert(all(full(solution_flat) < 1e9));
            assert(all(full(solution_flat) > -1e9));
            assert(all(full(solution_flat) < full(upper_bounds_flat) + 1e-7));
            assert(all(full(solution_flat) > full(lower_bounds_flat) - 1e-7));
            
            solution_values = obj.variables_unpack.call(struct('variables_flat', solution_flat));
            solution_values = map_struct(solution_values, @(e)full(e));
            
            delta_values = obj.variables_unpack.call(struct('variables_flat', solution_flat - initialization_flat));
            delta_values = map_struct(delta_values, @(e)full(e));
            
            slack_values = obj.slacks_unpack.call(struct('slacks_flat', slacks_flat));
            slack_values = map_struct(slack_values, @(e)full(e));
            
            % Unpack hard constraints, to see which are active
            g_hard_value = obj.g_hard_fn.call({solution_flat, parameters_flat_values});
            g_hard_unpacked = obj.g_hard_unpack_fn.call(struct('g_hard',g_hard_value{1}));
            g_hard_unpacked = map_struct(g_hard_unpacked, @(e)full(e));
            
            results = struct;
            results.solution = solution_values;
            results.g_hard_unpacked = g_hard_unpacked;
            results.delta_values = delta_values;
            results.delta_flat = full(solution_flat - initialization_flat);
            results.max_delta = max(abs(results.delta_flat));
            results.solution_flat = full(solution_flat);
            results.slacks = slack_values;
            results.max_infeasibility_solution = max_infeasibility_solution;
            results.raw_solution = map_struct(nlp_result, @(e)full(e));
            results.max_slack = max(slacks_flat);
            results.stats = obj.nlpsolver.stats();
            results.objective_terms = objective_term_values;
        end
    end
end

