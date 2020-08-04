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

function initialization = OCP_initialize(ocp, initialization_period_shift)
    
    info = ocp.model_info;
    initialization = ocp.solver.get_variables_struct(nan);

    %% Create initialization
    % The two aircraft are placed on the reference circle, 180 degrees apart.
    
    init_speed = 0.6 * ocp.parameters.min_airspeed(1) + 0.4 * ocp.parameters.max_airspeed(1);
    init_radius = ocp.parameters.p(info.ip_length1) * 0.65;
    init_period = 2*pi* init_radius / init_speed;
    
    
    initialization.t(:) = init_period;
    initialization.reelout_phase_fraction(:) = 0.8;
    
    interval_durations = full(ocp.interval_duration_fn(initialization.t, initialization.reelout_phase_fraction));
    interval_duration_sum = cumsum([0; interval_durations]);
    init_t_grid = interval_duration_sum(1:end-1)' + interval_durations' .* ocp.collocation.nodes(1:end-1);
    init_t_grid = init_t_grid(:);
    
    init_ref_circle_center = [450; 1; -300];

    init_trajectory = initialization_generalized_coordinates(...
        init_ref_circle_center, ...
        init_radius, ...
        init_t_grid, ...
        init_period, ...
        full(ocp.parameters.p(ocp.model_info.ip_length1)),...
        initialization_period_shift);

    initialization.t(:) = init_period;

    initialization.x(info.ix_theta0,      :) = init_trajectory.theta0;
    initialization.x(info.ix_phi0,        :) = init_trajectory.phi0;
    initialization.x(info.ix_theta,       :) = init_trajectory.theta;
    initialization.x(info.ix_phi,         :) = init_trajectory.phi;
    initialization.x(info.ix_theta0_dot,  :) = init_trajectory.theta0_dot;
    initialization.x(info.ix_phi0_dot,    :) = init_trajectory.phi0_dot;
    initialization.x(info.ix_theta_dot,   :) = init_trajectory.theta_dot;
    initialization.x(info.ix_phi_dot,     :) = init_trajectory.phi_dot;
    initialization.x(info.ix_length0,     :) = init_trajectory.length0;
    initialization.x(info.ix_length0_dot, :) = init_trajectory.length0_dot;
    
    

    initialization.z(info.iz_theta0_dot_dot, :) = init_trajectory.theta0_dot_dot;
    initialization.z(info.iz_phi0_dot_dot,   :) = init_trajectory.phi0_dot_dot;
    initialization.z(info.iz_theta_dot_dot,  :) = init_trajectory.theta_dot_dot;
    initialization.z(info.iz_phi_dot_dot,    :) = init_trajectory.phi_dot_dot;

    initialization.u(info.iu_length0_dot_dot,       :) = init_trajectory.length0_dot_dot;

    % Zero-initialize the uninitialized variables
    initialization = map_struct(initialization, @(ini) nan_to_zero(ini));
    
end

function x = nan_to_zero(x)
    x(isnan(x)) = 0;
end

function result = initialization_generalized_coordinates(ref_circle_center, ref_circle_radius, t_grid, T_period, length1, initialization_period_shift)

    
    Lc = sqrt(length1.^2 - ref_circle_radius.^2);
    Lr = norm(ref_circle_center);
    length0 = Lr - Lc;
    position0 = ref_circle_center / Lr * length0;
    position0_dot = [0;0;0];
    position0_dot_dot = [0;0;0];
    
    [length0, theta0, phi0, ~, ~, ~, ~, ~, ~] = convert_cartesian_to_spherical(position0,position0_dot,position0_dot_dot);
    

    result = struct;
    result.theta0           = theta0 * ones(size(t_grid));
    result.phi0             = phi0 * ones(size(t_grid));
    result.theta0_dot       = zeros(size(t_grid));
    result.phi0_dot         = zeros(size(t_grid));
    result.length0          = length0 * ones(size(t_grid));
    result.length0_dot      = zeros(size(t_grid));
    result.theta0_dot_dot   = zeros(size(t_grid));
    result.phi0_dot_dot     = zeros(size(t_grid));
    result.length0_dot_dot  = zeros(size(t_grid));
    result.t_grid           = t_grid;
    
    e1 = ref_circle_center ./ norm(ref_circle_center);
    e2 = [0;1;0];
    e2 = e2-(e2'*e1)*e1;
    e2 = e2 ./ norm(e2);
    e3 = cross(e1,e2);
    R = [e1 e2 e3];
    
    omega = 2*pi/T_period;
        
    for k = 1:length(t_grid)
        for i = 1:n_aircraft
            t = t_grid(k) + (initialization_period_shift + (i-1)/n_aircraft) * T_period;

            c = cos(omega * t);
            s = sin(omega * t);

            c_dot = -sin(omega * t)*(omega);
            s_dot = cos(omega * t)*(omega);

            c_dot_dot = -cos(omega * t)*(omega)^2;
            s_dot_dot = -sin(omega * t)*(omega)^2;

            position = R * [Lc; ref_circle_radius*c; ref_circle_radius*s];
            position_dot = R * [0; ref_circle_radius*c_dot; ref_circle_radius*s_dot];
            position_dot_dot = R * [0; ref_circle_radius*c_dot_dot; ref_circle_radius*s_dot_dot];

            [~, result.theta(i,k), result.phi(i,k),...
                ~, result.theta_dot(i,k), result.phi_dot(i,k),...
                ~, result.theta_dot_dot(i,k), result.phi_dot_dot(i,k)] = ...
                convert_cartesian_to_spherical(position,position_dot,position_dot_dot);
        end
    end
    
end

