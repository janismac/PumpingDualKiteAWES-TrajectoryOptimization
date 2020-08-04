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

function OCP_runall(sweep_config)

    if ~exist(sweep_config.output_directory, 'dir')
        mkdir(sweep_config.output_directory);
    end

    ocp_config = struct;
    ocp_config.is_ODE_constraint_soft = true;
    ocp_config.is_path_constraint_soft = true;
    ocp_config.has_power_maximization = false;
    ocp_config.has_reference_path = true;
    ocp_config.has_winch_speed_bias_penalty = false;
    ocp_config.n_nodes = sweep_config.collocation_n_nodes;
    ocp_config.n_intervals = max(round(sweep_config.collocation_n_timesteps/(sweep_config.collocation_n_nodes-1)),2);

    ocp = OCP_create(ocp_config);
    ocp.solver = NLP_Solver(ocp.nlp, sweep_config);
    ocp.units = units;
    ocp = OCP_set_parameters(ocp);
    ocp.parameters.max_curvature(:) = sweep_config.max_curvature;
    ocp.parameters.p(ocp.model_info.ip_wind_ref_speed) = sweep_config.wind_speed;
    ocp.parameters.p(ocp.model_info.ip_length) = sweep_config.secondary_tether_length;
    ocp.parameters.p(ocp.model_info.ip_tether_drag_diameter) = sweep_config.tether_diameter;
    ocp.upper_bounds.x(ocp.model_info.ix_theta0, :) = sweep_config.theta0_upper_bound;
    ocp.parameters.minimum_aircraft_distance(:) = sweep_config.minimum_aircraft_distance;
    
    ocp.parameters.path_constraint_weight(:)          = sweep_config.path_constraint_weight;
    ocp.parameters.ode_weight(:)                      = sweep_config.ode_weight;
    
    ocp.initialization = OCP_initialize(ocp, sweep_config.initialization_period_shift);
    ocp.parameters.x_ref(:,:) = ocp.initialization.x;
    ocp.parameter_units.x_ref(:,:) = ocp.variable_units.x;
    
    % Fix the phase split during the initialization problem
    ocp.lower_bounds.reelout_phase_fraction(:) = ocp.initialization.reelout_phase_fraction;
    ocp.upper_bounds.reelout_phase_fraction(:) = ocp.initialization.reelout_phase_fraction;

    checkpoint = struct;
    checkpoint.variables = ocp.initialization;
    checkpoint.solver_output = [];

    % Soft trust region loop
    for a = [10.^(5:-0.5:-5) 0]
        disp(['Relaxed dynamics OCP with weight: ' num2str(a)]);
        checkpoint(end+1) = OCP_step(ocp, checkpoint(end), a);
        save([sweep_config.output_directory '/OCP_checkpoint.mat'],'checkpoint');
    end
    

    % Calculate secondary results
    evaluation = OCP_evaluate(ocp, checkpoint(end));
    openloop_3sec_position_error = OCP_validate_openloop(ocp, evaluation);
    [~, tracking] = tracking_validation_lagrangian_inelastic(...
        evaluation, ocp.parameters.p, 1.1 * ocp.upper_bounds.u(:,1), 10);
    
    
    % Print some stats
    fprintf('max_slack:                       %e\n', checkpoint(end).solver_output.max_slack);
    fprintf('openloop_3sec_position_error:    %f\n', openloop_3sec_position_error);
    fprintf('tracking.max_dist_deviation:     %f\n', tracking.max_distance_deviation);
    fprintf('tracking.median_dist_deviation:  %f\n', tracking.median_distance_deviation);
    fprintf('mean_power:                      %f W\n', mean(evaluation.power));

    
    % Save results
    ocp_ = rmfield(ocp, 'solver');
    ocp_ = rmfield(ocp_, 'nlp');
    save([sweep_config.output_directory '/OCP_results.mat'],...
        'ocp_','evaluation','checkpoint','ocp_config','sweep_config');

    
    % Result validation
    assert(checkpoint(end).solver_output.max_slack < 1e-4, ...
        'Large slacks, problem not solved.');
    
    assert(openloop_3sec_position_error < 3, ...
        'Large open loop validation error.');
    
    assert(tracking.median_distance_deviation < 2, ...
        'Large tracking validation error.');
    
    % OCP_render(ocp_, evaluation)
    
    %% Second OCP, power cycle with multiple revolutions %%
    n_revolutions = sweep_config.n_revolutions;
    ocp2_config = struct;
    ocp2_config.is_ODE_constraint_soft = false;
    ocp2_config.is_path_constraint_soft = false;
    ocp2_config.has_power_maximization = true;
    ocp2_config.has_reference_path = false;
    ocp2_config.has_winch_speed_bias_penalty = true;
    ocp2_config.n_nodes = sweep_config.collocation_n_nodes;
    ocp2_config.n_intervals = n_revolutions * max(round(sweep_config.collocation_n_timesteps/(sweep_config.collocation_n_nodes-1)),2);

    ocp2 = OCP_create(ocp2_config);
    ocp2.solver = NLP_Solver(ocp2.nlp, sweep_config);
    ocp2.units = units;
    ocp2 = OCP_set_parameters(ocp2);
    ocp2.parameters.p = ocp_.parameters.p;
    ocp2.parameters.max_curvature(:) = sweep_config.max_curvature;
    ocp2.parameters.minimum_aircraft_distance(:) = sweep_config.minimum_aircraft_distance;
    ocp2.upper_bounds.x(ocp2.model_info.ix_theta0, :) = sweep_config.theta0_upper_bound;
    ocp2.initialization = OCP_initialize_multiple_revolutions(ocp2, evaluation, n_revolutions);
    
    checkpoint2 = struct;
    checkpoint2.variables = ocp2.initialization;
    checkpoint2.solver_output = [];
    
    
    % Soft trust region loop
    for a = [10.^(5:-0.5:-5) 0]
        disp(['Power cycle OCP with weight: ' num2str(a)]);
        checkpoint2(end+1) = OCP_step(ocp2, checkpoint2(end), a);
        save([sweep_config.output_directory '/OCP_checkpoint2.mat'],'checkpoint2');
    end
    
    % Calculate secondary results
    evaluation2 = OCP_evaluate(ocp2, checkpoint2(end));
    openloop_3sec_position_error2 = OCP_validate_openloop(ocp2, evaluation2);
    [~, tracking2] = tracking_validation_lagrangian_inelastic(...
        evaluation2, ocp2.parameters.p, 1.1 * ocp2.upper_bounds.u(:,1), 10);
    
    
    % Print some stats
    fprintf('openloop_3sec_position_error2:   %f\n', openloop_3sec_position_error2);
    fprintf('tracking2.max_dist_deviation:    %f\n', tracking2.max_distance_deviation);
    fprintf('tracking2.median_dist_deviation: %f\n', tracking2.median_distance_deviation);
    fprintf('mean_power:                      %f W\n', mean(evaluation2.power));

    % Save results
    ocp2_ = rmfield(ocp2, 'solver');
    ocp2_ = rmfield(ocp2_, 'nlp');
    save([sweep_config.output_directory '/OCP_results2.mat'],...
        'ocp2_','evaluation2','checkpoint2','ocp2_config','sweep_config');
    
    
    % Result validation
    assert(openloop_3sec_position_error2 < 3, ...
        'Large open loop validation error.');
    
    assert(tracking2.max_distance_deviation < 3, ...
        'Large tracking validation error.');
    
    assert(checkpoint2(end).solver_output.stats.success, ...
        'NLP solver reported failure.');
    
    % OCP_render(ocp2_, evaluation2);

end
