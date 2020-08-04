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

function ocp = OCP_create(ocp_config)
    
    import casadi.*
    
    model_info = load('models/model_lagrangian_inelastic_info');
    
    n_nodes = ocp_config.n_nodes;
    n_intervals = ocp_config.n_intervals;
    n_timesteps = n_intervals * (n_nodes-1);
    
    n_intervals_reelout = floor(n_intervals * 0.8);
    n_intervals_reelin = n_intervals - n_intervals_reelout;
    
    assert(n_intervals_reelout >= 1);
    assert(n_intervals_reelin >= 1);
    
    collocation = collocation_constants(n_nodes);
    
    interval_slices = (1:collocation.n_nodes) + (collocation.n_nodes-1)*(0:n_intervals-1)';
    interval_slices(end) = 1;
    
    nlp = NLP_Builder();
    x = nlp.add_variable('x', model_info.x_names, n_timesteps);
    z = nlp.add_variable('z', model_info.z_names, n_timesteps);
    u = nlp.add_variable('u', model_info.u_names, n_timesteps);
    u_dot = nlp.add_variable('u_dot', sym_names(total_time_derivative(sym(model_info.u_names))), n_timesteps);
    t = nlp.add_variable('t');
    reelout_phase_fraction = nlp.add_variable('reelout_phase_fraction');
    p = nlp.add_parameter('p', model_info.p_names);
    
    interval_duration = casadi.SX.zeros(n_intervals, 1);
    interval_duration(1:n_intervals_reelout) = t * reelout_phase_fraction / n_intervals_reelout;
    interval_duration(n_intervals_reelout + (1:n_intervals_reelin)) = t * (1-reelout_phase_fraction) / n_intervals_reelin;
    
    interval_fractions = casadi.SX.zeros(n_intervals, 1);
    interval_fractions(1:n_intervals_reelout) = reelout_phase_fraction / n_intervals_reelout;
    interval_fractions(n_intervals_reelout + (1:n_intervals_reelin)) = (1-reelout_phase_fraction) / n_intervals_reelin;
    
    

    timestep_integration_weights = casadi.SX.zeros(1,n_timesteps);
    for i = 1:n_intervals
        for j = 1:n_nodes
            timestep_integration_weights(interval_slices(i,j)) = timestep_integration_weights(interval_slices(i,j)) + ...
                interval_duration(i) * collocation.integration_weights(j);
        end
    end
    timestep_integration_weights = simplify(timestep_integration_weights);
    
    
    min_airspeed = nlp.add_parameter('min_airspeed');
    max_airspeed = nlp.add_parameter('max_airspeed');
    
    min_tension_main_tether = nlp.add_parameter('min_tension_main_tether');
    max_tension_main_tether = nlp.add_parameter('max_tension_main_tether');
    
    min_tension_secondary_tether = nlp.add_parameter('min_tension_secondary_tether');
    max_tension_secondary_tether = nlp.add_parameter('max_tension_secondary_tether');
    
    min_proper_acceleration_body = nlp.add_parameter('min_proper_acceleration_body', 3);
    max_proper_acceleration_body = nlp.add_parameter('max_proper_acceleration_body', 3);
    
    input_rate_weight = nlp.add_parameter('input_rate_weight');
    minimum_aircraft_distance = nlp.add_parameter('minimum_aircraft_distance');
    minimum_altitude = nlp.add_parameter('minimum_altitude');
    
    max_curvature = nlp.add_parameter('max_curvature');
    
    if ocp_config.has_power_maximization
        average_power_weight = nlp.add_parameter('average_power_weight');
        maximum_power = nlp.add_parameter('maximum_power');
    end
    
    if ocp_config.has_reference_path
        reference_trajectory_weight = nlp.add_parameter('reference_trajectory_weight');
        x_ref = nlp.add_parameter('x_ref', size(x,1), size(x,2));
    end
    
    if ocp_config.is_ODE_constraint_soft
        ode_weight = nlp.add_parameter('ode_weight');
        ode_weights = ode_weight * repmat(interval_duration.',n_nodes-1,1) .* repmat(collocation.integration_weights(2:end),1,n_intervals);
    else
        ode_weight = inf;
        ode_weights = inf(n_nodes-1, n_intervals);
    end
    
    
    if ocp_config.is_path_constraint_soft
        path_constraint_weights = nlp.add_parameter('path_constraint_weight') * timestep_integration_weights;
    else
        path_constraint_weights = inf(1,n_timesteps);
    end
    
    if ocp_config.has_winch_speed_bias_penalty
        phase_switch_speed = nlp.add_variable('phase_switch_speed');
        phase_switch_weight = nlp.add_parameter('phase_switch_weight');
    end
    
    
    % Evaluate the model for every timestep
    casadi_model_wrapper_fn = casadi_model_wrapper();
    for k = 1:n_timesteps
        model(k) = casadi_model_wrapper_fn.call(...
            struct('x',x(:,k),'u',u(:,k),'p',p,'z',z(:,k)));
    end
    
    % Apply the ODE collocation for each interval
    for i = 1:n_intervals
        dxdt_interval = [model(interval_slices(i,:)).dxdt].';
        x_interval = x(:,interval_slices(i,:))';
        delta_x = x_interval(2:end,:) - repmat(x_interval(1,:), n_nodes-1,1);
        dxdt_integrals = interval_duration(i) * (collocation.integration_matrix * dxdt_interval);
        nlp.add_constraint(['defect_ode_' num2str(i)], delta_x - dxdt_integrals, 0, 0, ode_weights(:,i));
        
        u_dot_interval = u_dot(:,interval_slices(i,:))';
        u_interval = u(:,interval_slices(i,:))';
        delta_u = u_interval(2:end,:) - repmat(u_interval(1,:), n_nodes-1,1);
        u_dot_integrals = interval_duration(i) * (collocation.integration_matrix * u_dot_interval);
        nlp.add_constraint(['defect_input_ode_' num2str(i)], delta_u - u_dot_integrals, 0, 0, ode_weights(:,i));
    end
    
    % Lagrangian mechanics EOM for the tether system
    nlp.add_constraint('defect_EOM', [model.defect_EOM], 0, 0, inf);
    
    % Tether tension limits
    nlp.add_constraint('tension_main_tether', [model.tether_tension0], min_tension_main_tether, max_tension_main_tether, path_constraint_weights);
    nlp.add_constraint('tension_secondary_tether', [model.tether_tension], min_tension_secondary_tether, max_tension_secondary_tether, path_constraint_weights);
    
    % Airspeed limits
    nlp.add_constraint('airspeed', vertcat(model.aero_speed).', min_airspeed, max_airspeed, path_constraint_weights);
    
    % Reel-in/out phase bias
    if ocp_config.has_winch_speed_bias_penalty
        for i = 1:n_intervals
            L0dot = x(model_info.ix_length0_dot,interval_slices(i,:));
            penalty = L0dot - phase_switch_speed;
            if i <= n_intervals_reelout
                penalty = -penalty;
            end
            
            nlp.add_constraint(['phase_switch_penalty' num2str(i)], penalty, ...
                -inf, 0, phase_switch_weight * interval_duration(i) * collocation.integration_weights');
        end
    end
    
    % Penalize input rate (du/dt)^2
    for i = 1:n_intervals
        U_dot_interval = u_dot(:,interval_slices(i,:))';
        integral_dudt_squared = interval_fractions(i) * (collocation.integration_weights' * U_dot_interval.^2);
        nlp.add_objective(['input_rate' num2str(i)], input_rate_weight * sum(integral_dudt_squared));
    end
    
    % Minimum altitude
    aircraft_position_z = casadi.SX.zeros(n_aircraft, n_timesteps);
    for k = 1:n_timesteps
        for i = 1:n_aircraft
            aircraft_position_z(i,k) = model(k).position(3,i);
        end
    end
    nlp.add_constraint('minimum_altitude_constraint', minimum_altitude + aircraft_position_z, -inf, 0, inf);
    
    
    % Penalize reference path deviation
    if ocp_config.has_reference_path
        slice = model_info.ix_theta0:model_info.ix_length0_dot;
        penalty = sum1((x(slice,:)-x_ref(slice,:)).^2);
        nlp.add_objective('path_deviation', reference_trajectory_weight * sum(timestep_integration_weights .* penalty) );
    end
    
    % Constrain the pairwise distances between aircraft
    for i = 1:n_aircraft
        for j = (i+1):n_aircraft
            distances = casadi.SX.zeros(1, n_timesteps);
            for k = 1:n_timesteps
                position_delta = model(k).position(:,i) - model(k).position(:,j);
                distances(k) = norm(position_delta);
            end
            
            nlp.add_constraint(['distance_constraint_' num2str(i) '_' num2str(j)], ...
                minimum_aircraft_distance - distances, -inf, 0, inf);
        end
    end
    
    % Average power maximization, maximum power constraint
    if ocp_config.has_power_maximization
        main_tether_tensions = [model.tether_tension0];
        instantaneous_power = main_tether_tensions .* x(model_info.ix_length0_dot,:);
        mean_power = cell(n_intervals,1);
        for i = 1:n_intervals
            instantaneous_power_interval = instantaneous_power(:,interval_slices(i,:))';
            mean_power{i} = collocation.integration_weights' * instantaneous_power_interval;
        end
        mean_power = [mean_power{:}];
        mean_power = mean_power * interval_fractions;
        nlp.add_objective('mean_power', -average_power_weight * mean_power);
        
        
        nlp.add_constraint('maximum_power_constraint', ...
            instantaneous_power, -inf, maximum_power, inf);
    end
    
    % Tether/aircraft collision avoidance
    for i = 1:n_aircraft
        d_tether_dot_d_lift = casadi.SX.zeros(1, n_timesteps);
        for k = 1:n_timesteps
            d_tether_dot_d_lift(k) = model(k).direction(:,i)' * model(k).aero_direction_lift(:,i);
        end

        nlp.add_constraint(['tether_collision_constraint' num2str(i)], -d_tether_dot_d_lift, -inf, -cos(55/180*pi), inf);
    end
    
    % Path curvature constraint
    curvatures = casadi.SX.zeros(n_aircraft, n_timesteps);
    for k = 1:n_timesteps
        curvatures(:,k) = model(k).curvature.';
    end
    nlp.add_constraint('curvature_constraint', curvatures, -inf, max_curvature, inf);
    
    % Load limit constraints
    for i = 1:n_aircraft
        for d = 1:3
            proper_acceleration_body_component = casadi.SX.zeros(1, n_timesteps);
            for k = 1:n_timesteps
                R = reshape(model(k).orientation_matrix_flat(:,i), 3, 3);
                proper_acceleration_body = (R.') * model(k).proper_acceleration(:,i);
                proper_acceleration_body_component(k) = proper_acceleration_body(d);
            end
            
            nlp.add_constraint(['proper_acceleration_constraint_'  num2str(i) '_'  num2str(d)], ...
                proper_acceleration_body_component, ...
                min_proper_acceleration_body(d), max_proper_acceleration_body(d), inf);
        end
    end
    
    ocp = struct;
    ocp.nlp = nlp;
    ocp.model_info      = model_info;
    ocp.collocation     = collocation;
    ocp.n_intervals     = n_intervals;
    ocp.n_timesteps     = n_timesteps;
    ocp.interval_slices = interval_slices;
    
    ocp.interval_duration_fn = casadi.Function(...
        'interval_duration_fn', ...
        {t, reelout_phase_fraction}, ...
        {interval_duration}, ...
        {'t', 'reelout_phase_fraction'}, ...
        {'interval_duration'});
end


function fn = casadi_model_wrapper()
    model_info = load('models/model_lagrangian_inelastic_info');

    x = casadi.SX.sym('x', length(model_info.x_names));
    z = casadi.SX.sym('z', length(model_info.z_names));
    u = casadi.SX.sym('u', length(model_info.u_names));
    p = casadi.SX.sym('p', length(model_info.p_names));
    w = 0*casadi.SX.sym('w', length(model_info.w_names));
    
    
    aero_inputs = model_lagrangian_inelastic_aero_inputs_wrapper_generated(x,u,p,w);
    aero_outputs = model_aero_forces(aero_inputs);
    [M,f,C] = model_lagrangian_inelastic_generated(x,u,p,w,aero_outputs.applied_force0,aero_outputs.applied_force);
    
    M = M.sparsify();
    f = f.sparsify();
    C = sparse(C);
    
    defect_EOM = M*z - f;
    defect_EOM = defect_EOM.simplify();
    dxdt = C * [x; z; u];
    
    outputs = model_lagrangian_inelastic_outputs_wrapper_generated(x,u,p,w,z,aero_outputs.applied_force0,aero_outputs.applied_force);
        
    fn = casadi.Function('model_wrapper', ...
        {x,u,p,z}, ...
        vertcat({defect_EOM; dxdt}, struct2cell(outputs), struct2cell(aero_outputs)),...
        {'x','u','p','z'},...
        vertcat({'defect_EOM'; 'dxdt'}, fieldnames(outputs), fieldnames(aero_outputs)));
end
