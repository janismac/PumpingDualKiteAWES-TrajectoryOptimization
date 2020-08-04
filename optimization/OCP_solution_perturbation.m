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

function OCP_solution_perturbation(parameter_name)

    perturbation_config = struct;
    perturbation_config.parameter_name = parameter_name;


    if strcmp(parameter_name, 'period')
        perturbation_config.increment = 0.5;
        perturbation_config.min_value = 40.0;
        perturbation_config.max_value = 80.0;
    end
    
    if strcmp(parameter_name, 'wind_speed')
        perturbation_config.increment = 0.2;
        perturbation_config.min_value = 4.0;
        perturbation_config.max_value = 25.0;
    end
    
    if strcmp(parameter_name, 'secondary_tether_length')
        perturbation_config.increment = 2;
        perturbation_config.min_value = 50;
        perturbation_config.max_value = 200;
    end
    
    if strcmp(parameter_name, 'minimum_aircraft_distance')
        perturbation_config.increment = 2;
        perturbation_config.min_value = 10;
        perturbation_config.max_value = 250;
    end
    
    if strcmp(parameter_name, 'max_power')
        perturbation_config.increment = 200;
        perturbation_config.min_value = 10000;
        perturbation_config.max_value = 35600;
    end
    
    
    assert(exist('output','dir') ~= 0)
    assert(exist('output/sweep_perturbation_base','dir') ~= 0)
    base_solution = load('output/sweep_perturbation_base/OCP_results2.mat');
    
    solution_checkpoint = base_solution.checkpoint2(end);
    ocp2_config = base_solution.ocp2_config;
    sweep_config = base_solution.sweep_config;
    
    ocp2 = OCP_create(ocp2_config);
    ocp2 = combine_structs(ocp2, base_solution.ocp2_);
    ocp2.solver = NLP_Solver(ocp2.nlp, sweep_config);

    %% Perturb the optimal trajectory, force a parameter up and down
    if strcmp(perturbation_config.parameter_name, 'period')
        perturbation_start_value = solution_checkpoint.variables.t;
    elseif strcmp(perturbation_config.parameter_name, 'wind_speed')
        perturbation_start_value = ocp2.parameters.p(ocp2.model_info.ip_wind_ref_speed);
    elseif strcmp(perturbation_config.parameter_name, 'secondary_tether_length')
        perturbation_start_value = ocp2.parameters.p(ocp2.model_info.ip_length1);
    elseif strcmp(perturbation_config.parameter_name, 'minimum_aircraft_distance')
        perturbation_start_value = ocp2.parameters.minimum_aircraft_distance;
    elseif strcmp(perturbation_config.parameter_name, 'max_power')
        perturbation_start_value = ocp2.parameters.maximum_power;
    else
        error('not implemented');
    end
    
    perturbation_down_values = (perturbation_start_value-perturbation_config.increment):(-perturbation_config.increment):perturbation_config.min_value;
    perturbation_up_values = (perturbation_start_value+perturbation_config.increment):perturbation_config.increment:perturbation_config.max_value;
    
    assert(numel(perturbation_down_values) > 1)
    assert(numel(perturbation_up_values) > 1)
    
    checkpoints_perturbation_up = solution_checkpoint;
    checkpoints_perturbation_down = solution_checkpoint;
    
    for i = 1:length(perturbation_up_values)
        c = checkpoints_perturbation_up(end);
        assert(c.solver_output.stats.success);
        [ocp2_i, c] = apply_parameter_perturbation(ocp2, c, perturbation_up_values(i), perturbation_config);
        checkpoints_perturbation_up(end+1) = OCP_step(ocp2_i, c, 0);
        ocp_up(i) = ocp_delete_casadi_objects(ocp2_i);
    end
    
    for i = 1:length(perturbation_down_values)
        c = checkpoints_perturbation_down(end);
        assert(c.solver_output.stats.success);
        [ocp2_i, c] = apply_parameter_perturbation(ocp2, c, perturbation_down_values(i), perturbation_config);
        checkpoints_perturbation_down(end+1) = OCP_step(ocp2_i, c, 0);
        ocp_down(i) = ocp_delete_casadi_objects(ocp2_i);
    end
    
    perturbation_values = [fliplr(perturbation_down_values) perturbation_start_value perturbation_up_values];
    perturbation_checkpoints = [fliplr(checkpoints_perturbation_down(2:end)) solution_checkpoint checkpoints_perturbation_up(2:end)];
    perturbation_ocp = [fliplr(ocp_down) ocp_delete_casadi_objects(ocp2) ocp_up];
    
    
    disp('Saving OCP solutions');
    save(['output/sweep_perturbation_base/OCP_' perturbation_config.parameter_name '_perturbation__no_eval.mat'],...
        'perturbation_ocp','ocp2_config','sweep_config','perturbation_values','perturbation_checkpoints','perturbation_config',...
        '-v7.3');
    
    
    disp('Running trajectory interpolation/evaluation');
    for i = 1:length(perturbation_checkpoints)
        perturbation_evaluation(i) = OCP_evaluate(perturbation_ocp(i), perturbation_checkpoints(i), true);
    end
    
    
    disp('Saving trajectory interpolation/evaluation');
    save(['output/sweep_perturbation_base/OCP_' perturbation_config.parameter_name '_perturbation.mat'],...
        'perturbation_ocp','ocp2_config','sweep_config','perturbation_values','perturbation_checkpoints','perturbation_evaluation','perturbation_config',...
        '-v7.3');
    
    disp('Running trajectory validation');
    for i = 1:length(perturbation_checkpoints)
        openloop_3sec_position_error2 = OCP_validate_openloop(perturbation_ocp(i), perturbation_evaluation(i));
        [~, tracking2] = tracking_validation_lagrangian_inelastic(...
            perturbation_evaluation(i), perturbation_ocp(i).parameters.p, 2 * perturbation_ocp(i).upper_bounds.u(:,1), 10);

        fprintf('openloop_3sec_position_error2:   %f\n', openloop_3sec_position_error2);
        fprintf('tracking2.max_dist_deviation:    %f\n', tracking2.max_distance_deviation);

        % Result validation
        assert(openloop_3sec_position_error2 < 3, ...
            'Large open loop validation error.');

        assert(tracking2.max_distance_deviation < 3, ...
            'Large tracking validation error.');
    end
    
end

function ocp = ocp_delete_casadi_objects(ocp)
    ocp = rmfield(ocp, 'solver');
    ocp = rmfield(ocp, 'nlp');
end

function [ocp2, checkpoint] = apply_parameter_perturbation(ocp2, checkpoint, value, perturbation_config)
    
    if strcmp(perturbation_config.parameter_name, 'period')
        checkpoint.variables.t(:) = value;
        ocp2.lower_bounds.t(:) = value;
        ocp2.upper_bounds.t(:) = value;
    elseif strcmp(perturbation_config.parameter_name, 'wind_speed')
        ocp2.parameters.p(ocp2.model_info.ip_wind_ref_speed) = value;
    elseif strcmp(perturbation_config.parameter_name, 'secondary_tether_length')
        ocp2.parameters.p(ocp2.model_info.ip_length) = value;
    elseif strcmp(perturbation_config.parameter_name, 'minimum_aircraft_distance')
        ocp2.parameters.minimum_aircraft_distance(:) = value;
    elseif strcmp(perturbation_config.parameter_name, 'max_power')
        ocp2.parameters.maximum_power(:) = value;
    else
        error('not implemented');
    end
    
end