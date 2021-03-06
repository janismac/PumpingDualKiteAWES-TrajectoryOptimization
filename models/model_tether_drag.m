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

function r = model_tether_drag(...
    position0, ...
    position, ...
    position0_dot, ...
    position_dot, ...
    direction0, ...
    direction, ...
    tether_length_0, ...
    secondary_length, ...
    tether_drag_diameter, ...
    air_density_msl, ...
    scale_height, ...
    wind_profile_exponent, ...
    wind_ref_height, ...
    wind_ref_speed, ...
    wind_disturbance0, ...
    wind_disturbance ...
)

    
    r = struct(...
        'lumped_drag_forces_0', nan(3,1), ...
        'lumped_drag_forces', nan(3,n_aircraft) ...
    );
    
    
    % Main tether
    [q00, ~] = lumped_dynamic_pressure_integral(...
        position0, [0;0;0], ...
        position0_dot, [0;0;0], ...
        wind_disturbance0, [0;0;0], ...
        direction0, air_density_msl, scale_height, wind_profile_exponent, wind_ref_height, wind_ref_speed);
    
    % Secondary tethers
    q0 = cell(n_aircraft, 1);
    q = cell(n_aircraft, 1);
    for i = 1:n_aircraft
        [q0{i}, q{i}] = lumped_dynamic_pressure_integral(...
            position0, position(:,i), ...
            position0_dot, position_dot(:,i), ...
            wind_disturbance0, wind_disturbance(:,i), ...
            direction(:,i), air_density_msl, scale_height, wind_profile_exponent, wind_ref_height, wind_ref_speed);
    end
    
    % Drag force from dynamic pressure: F_D = q * L * d
    r.lumped_drag_forces_0 = tether_drag_diameter * q00 * tether_length_0;

    lumped_drag_forces = cell(n_aircraft, 1);
    for i = 1:n_aircraft
        r.lumped_drag_forces_0 = r.lumped_drag_forces_0 ...
            + q0{i} * tether_drag_diameter * secondary_length(i);
        
        lumped_drag_forces{i} = ...
              q{i}  * tether_drag_diameter * secondary_length(i);
    end
    
    r.lumped_drag_forces = horzcat(lumped_drag_forces{:});
end

function [qA, qB] = lumped_dynamic_pressure_integral(...
    pA, ...
    pB, ...
    vA, ...
    vB, ...
    wind_disturbance_A, ...
    wind_disturbance_B, ...
    tether_direction, ...
    air_density_msl, ...
    scale_height, ...
    wind_profile_exponent, ...
    wind_ref_height, ...
    wind_ref_speed ...
)

    qA = [0;0;0];
    qB = [0;0;0];
    
    quadrature_nodes = ([0.11270166 0.5 0.88729833]);
    quadrature_weights = ([5 8 5]/18);
    
    for i = 1:length(quadrature_nodes)
        tau = quadrature_nodes(i);
        w = quadrature_weights(i);
        v_local = vA + tau*(vB - vA);
        wind_disturbance_local = wind_disturbance_A + tau*(wind_disturbance_B - wind_disturbance_A);
        height_local = -(pA(3) + tau*(pB(3) - pA(3)));
        
        atm_local = model_atmosphere( ...
            air_density_msl, ...
            height_local, ...
            scale_height, ...
            v_local, ...
            wind_profile_exponent, ...
            wind_ref_height, ...
            wind_ref_speed, ...
            wind_disturbance_local ...
        );
        
        aero_velocity_normal = atm_local.aero_velocity - (atm_local.aero_velocity.' * tether_direction) * tether_direction;
        
        aero_speed_normal = sqrt(sum(aero_velocity_normal.^2) + 1e-4);
        
        dynamic_pressure_vector = -0.5 * atm_local.air_density * ...
                  aero_velocity_normal * aero_speed_normal;

        qA = qA + w * (1-tau) * dynamic_pressure_vector;
        qB = qB + w * tau     * dynamic_pressure_vector;
    end
    
end

