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

function evaluation = OCP_evaluate(ocp, checkpoint, has_state_space_matrices)

    if nargin < 3
        has_state_space_matrices = true; 
    end

    info = ocp.model_info;
    evaluation = struct;
    
    % Interpolate results
    [evaluation.x, evaluation.t] = interpolate_collocation(ocp, checkpoint.variables.x, checkpoint.variables);
    evaluation.u = interpolate_collocation(ocp, checkpoint.variables.u, checkpoint.variables);
    evaluation.u_dot = interpolate_collocation(ocp, checkpoint.variables.u_dot, checkpoint.variables);
    
    
    w = zeros(length(info.w_names), 1);
    dt = median(diff(evaluation.t));
    assert(std(diff(evaluation.t)) < 1e-12);
    
    
    for i = 1:length(evaluation.t)
        aero_inputs = model_lagrangian_inelastic_aero_inputs_wrapper_generated(...
            evaluation.x(:,i), evaluation.u(:,i), ocp.parameters.p, w);
        
        aero_outputs(i) = model_aero_forces(aero_inputs);
        
        [M,f,C] = model_lagrangian_inelastic_generated(...
            evaluation.x(:,i), evaluation.u(:,i), ocp.parameters.p, w,aero_outputs(i).applied_force0,aero_outputs(i).applied_force);
    
        outputs(i) = model_lagrangian_inelastic_outputs_wrapper_generated(...
            evaluation.x(:,i), evaluation.u(:,i), ocp.parameters.p, w, M\f, aero_outputs(i).applied_force0, aero_outputs(i).applied_force);
        
        evaluation.z(:,i) = M\f;
    end
    
    
    
    % Unpack state space vectors
    info_names = fieldnames(info);
    for i = 1:length(info_names)
        name = info_names{i};
        if length(name) > 2 && name(1) == 'i' && name(3) == '_' && name(2) ~= 'p' && name(2) ~= 'w'
            V = evaluation.(name(2));
            I = info.(name);
            evaluation.(name(4:end)) = nan(length(I), length(evaluation.t));
            for k = 1:length(I)
                evaluation.(name(4:end))(k,:) = V(I(k),:);
            end
        end
    end
    
    
    
    if has_state_space_matrices
        state_space_fn = linearize_model_lagrangian_inelastic;
        for i = 1:length(evaluation.t)
            % Calculate continuous state space matrices
            [A,B,G] = state_space_fn(evaluation.x(:,i), evaluation.u(:,i), evaluation.z(:,i), ocp.parameters.p);
            evaluation.A(:,:,i) = A;
            evaluation.B(:,:,i) = B;
            evaluation.G(:,:,i) = G;
            evaluation.max_real_pole(i) = max(real(eig(A)));
            evaluation.unstable_pole_sum(i) = sum(max(real(eig(A)),0));
        end
        
        % Tracking controller gains based on finize horizon tracking LQR
        Q = ones(size(evaluation.x,1),1);
        Q([info.ix_theta0 info.ix_phi0 info.ix_theta info.ix_phi]) = 10;
        Q = diag(Q);
        R = 0.1*eye(size(evaluation.u,1));
        evaluation.K = tune_tracking_controller(evaluation.A,evaluation.B,dt,Q,R);
    end
    
    for i = 1:length(evaluation.t)
        evaluation.orientation_matrix(:,:,:,i) = reshape(aero_outputs(i).orientation_matrix_flat, 3, 3, n_aircraft);
    end
    
    % Estimate angular velocities using finite differences
    n_shift = 2;
    N = length(evaluation.t);
    for I = [circshift(2:N,-n_shift);circshift(2:N,0);circshift(2:N,n_shift)]
        
        i_plus = I(1);
        i = I(2);
        i_minus = I(3);
        
        for k = 1:n_aircraft
            dR = evaluation.orientation_matrix(:,:,k,i_plus) / evaluation.orientation_matrix(:,:,k,i_minus);
            dR = (dR - dR')/2;
            evaluation.angular_velocity(:,k,i) = [dR(3,2);dR(1,3);dR(2,1)] / (2*n_shift*dt);
            evaluation.angular_velocity_body(:,k,i) = evaluation.orientation_matrix(:,:,k,i)' * evaluation.angular_velocity(:,k,i);
        end
    end
    
    evaluation.power = [outputs.tether_tension0] .* evaluation.length0_dot;
    
    % Combine outputs into array
    output_fields = fieldnames(outputs);
    combined_outputs = struct;
    for i = 1:numel(output_fields)
        name = output_fields{i};
        I = ndims(outputs(1).(name));
        combined_outputs.(name) = cat(I+1,outputs.(name));
    end
    
    
    % Combine aero outputs into array
    aero_output_fields = fieldnames(aero_outputs);
    combined_aero_outputs = struct;
    for i = 1:numel(aero_output_fields)
        name = aero_output_fields{i};
        I = ndims(aero_outputs(1).(name));
        combined_aero_outputs.(name) = cat(I+1,aero_outputs.(name));
    end
    
    evaluation = combine_structs(combined_outputs, evaluation);
    evaluation = combine_structs(combined_aero_outputs, evaluation);
    
    % Calculate proper acceleration in body coordinates
    for i = 1:n_aircraft
        for k = 1:length(evaluation.t)
            R = evaluation.orientation_matrix(:,:,i,k);
            evaluation.proper_acceleration_body(:,i,k) = (R.') * evaluation.proper_acceleration(:,i,k);
        end
    end
    
    
    % Power of individual forces
    evaluation.power_aero_force = squeeze(sum(evaluation.position_dot .* evaluation.aero_force_vector,1));
    evaluation.power_lumped_drag_force = squeeze(sum(evaluation.position_dot .* evaluation.lumped_drag_forces,1));
    evaluation.power_lumped_drag_force0 = squeeze(sum(evaluation.position0_dot .* evaluation.lumped_drag_forces_0,1))';
end


function [X_interp2, t_interp2] = interpolate_collocation(ocp, node_values, variables)

    interval_durations = full(ocp.interval_duration_fn(variables.t, variables.reelout_phase_fraction));
    interval_duration_sum = cumsum([0; interval_durations]);
    t_grid = interval_duration_sum(1:end-1)' + interval_durations' .* ocp.collocation.nodes(1:end-1);
    t_grid = t_grid(:);

    assert(length(t_grid) == size(node_values, 2));
    
    X_interp = cell(ocp.n_intervals, 1);
    t_interp = cell(ocp.n_intervals, 1);
    
    for i = 1:ocp.n_intervals
        X = node_values(:, ocp.interval_slices(i,:));
        
        X_interp_interval = X * ocp.collocation.interpolation_matrix';
        t_interp_interval = interval_duration_sum(i) + ocp.collocation.interpolation_nodes * interval_durations(i);
        
        if i < ocp.n_intervals % dont repeat knots
            X_interp_interval = X_interp_interval(:, 1:end-1);
            t_interp_interval = t_interp_interval(:, 1:end-1);
        end
        
        X_interp{i} = X_interp_interval;
        t_interp{i} = t_interp_interval;
    end
    
    X_interp = horzcat(X_interp{:});
    t_interp = horzcat(t_interp{:});
    
    % Re-interpoalte to equidistant 10 ms grid
    assert(max(diff(t_interp)) < 0.1);
    t_interp2 = min(t_interp):0.01:max(t_interp);
    X_interp2 = interp1(t_interp,X_interp',t_interp2)';
end


