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

function ocp = OCP_set_parameters(ocp)
    
    ocp.lower_bounds = ocp.solver.get_variables_struct(-inf);
    ocp.upper_bounds = ocp.solver.get_variables_struct(inf);
    ocp.parameters = ocp.solver.get_parameters_struct(nan);
    
    ocp.parameter_units = ocp.solver.get_parameters_struct(1.0);
    ocp.variable_units = ocp.solver.get_variables_struct(1.0);
    
    %% All problem parameters are in SI units. 
    % Unit scaling factors are stored separately for later conversion.
    % Omitted units correspond to dimensionless values (e.g. radian).
    
    u = ocp.units;
    info = ocp.model_info;
    
    ocp.parameters.input_rate_weight(1) = 0.01;
    
    if isfield(ocp.parameters, 'reference_trajectory_weight')
        ocp.parameters.reference_trajectory_weight(1) = 50;
    end
    
    if isfield(ocp.parameters, 'average_power_weight')
        ocp.parameters.average_power_weight(1) = 10;

        ocp.parameters     .maximum_power(:) = 45000;
        ocp.parameter_units.maximum_power(:) = u.W;
    end

    if isfield(ocp.parameters, 'max_curvature')
        ocp.parameters     .max_curvature(1) = nan; % Set by sweep config
        ocp.parameter_units.max_curvature(1) = 1 / u.m;
    end
    
    if isfield(ocp.parameters, 'phase_switch_weight')
        ocp.parameters     .phase_switch_weight(1) = 100;
        ocp.parameter_units.phase_switch_weight(1) = 1;

        ocp.lower_bounds  .phase_switch_speed(:) =  -10;
        ocp.upper_bounds  .phase_switch_speed(:) =  0;
        ocp.variable_units.phase_switch_speed(:) =  u.mps;
    end
    
    ocp.parameters     .max_airspeed(1)              = 45;
    ocp.parameter_units.max_airspeed(1)              = u.mps;

    ocp.parameters     .min_airspeed(1)              = 16;
    ocp.parameter_units.min_airspeed(1)              = u.mps;

    ocp.parameters     .max_tension_main_tether(1)   = 2000 * n_aircraft;
    ocp.parameter_units.max_tension_main_tether(1)   = u.N;

    ocp.parameters     .min_tension_main_tether(1)   = 200;
    ocp.parameter_units.min_tension_main_tether(1)   = u.N;
    
    ocp.parameters     .max_tension_secondary_tether(1) = 2000;
    ocp.parameter_units.max_tension_secondary_tether(1) = u.N;
    
    ocp.parameters     .min_tension_secondary_tether(1) = 200;
    ocp.parameter_units.min_tension_secondary_tether(1) = u.N;
    
    
    ocp.parameters     .max_proper_acceleration_body(:) = [1.5; 3; 3] * 19.80665;
    ocp.parameter_units.max_proper_acceleration_body(:) = u.mps2;
    
    ocp.parameters     .min_proper_acceleration_body(:) = [-1.5; -3; -3] * 19.80665;
    ocp.parameter_units.min_proper_acceleration_body(:) = u.mps2;
    
    
    ocp.parameters     .p(info.ip_g0)               =  9.80665;
    ocp.parameter_units.p(info.ip_g0)               =  u.mps2;

    ocp.parameters     .p(info.ip_length)          =  nan; % Set by sweep config
    ocp.parameter_units.p(info.ip_length)          =  u.m;

    ocp.parameters     .p(info.ip_mass0)            =   1.0;
    ocp.parameter_units.p(info.ip_mass0)            =  u.kg;

    ocp.parameters     .p(info.ip_mass_aircraft)    =  36.8;
    ocp.parameter_units.p(info.ip_mass_aircraft)    =  u.kg;

    ocp.parameters     .p(info.ip_tether_density)   =  0.005;
    ocp.parameter_units.p(info.ip_tether_density)   =  u.kg / u.m;

    ocp.parameters     .p(info.ip_wind_ref_speed)   =  nan; % Set by sweep config
    ocp.parameter_units.p(info.ip_wind_ref_speed)   =  u.mps;

    ocp.parameters     .p(info.ip_wind_ref_height)  =  200.0;
    ocp.parameter_units.p(info.ip_wind_ref_height)  =  u.m;

    ocp.parameters     .p(info.ip_wind_profile_exponent) =  0.14;

    ocp.parameters     .p(info.ip_air_density_msl)  =  1.2243;
    ocp.parameter_units.p(info.ip_air_density_msl)  =  u.kg_p_m3;

    ocp.parameters     .p(info.ip_scale_height)     =  10780;
    ocp.parameter_units.p(info.ip_scale_height)     =  u.m;

    ocp.parameters     .p(info.ip_tether_drag_diameter) = nan; % Set by sweep config
    ocp.parameter_units.p(info.ip_tether_drag_diameter) = u.m;

    ocp.parameters     .p(info.ip_reference_area)   =  3.0;
    ocp.parameter_units.p(info.ip_reference_area)   =  u.m^2;

    ocp.parameters     .p(info.ip_CL0)              =  0.524972;
    ocp.parameters     .p(info.ip_CLalpha)          =  4.589055;
    ocp.parameters     .p(info.ip_CD0)              =  0.036616;
    ocp.parameters     .p(info.ip_CDalpha)          =  0.150993;
    ocp.parameters     .p(info.ip_CDalpha2)         =  0.738638;
    
    ocp.parameters     .minimum_aircraft_distance(1) = nan; % Set by sweep config
    ocp.parameter_units.minimum_aircraft_distance(1) = u.m;
    
    ocp.parameters     .minimum_altitude(:) = 100;
    ocp.parameter_units.minimum_altitude(:) = u.m;
    
    ocp.upper_bounds  .x(info.ix_theta0,            :) =  nan; % Set by sweep config

    ocp.lower_bounds  .x(info.ix_phi0,              :) = -1.4;
    ocp.upper_bounds  .x(info.ix_phi0,              :) =  1.4;

    ocp.lower_bounds  .x(info.ix_phi,               :) = -1.4;
    ocp.upper_bounds  .x(info.ix_phi,               :) =  1.4;

    ocp.variable_units.x(info.ix_theta0_dot,        :) =  1/u.s;
    ocp.variable_units.x(info.ix_phi0_dot,          :) =  1/u.s;
    ocp.variable_units.x(info.ix_theta_dot,         :) =  1/u.s;
    ocp.variable_units.x(info.ix_phi_dot,           :) =  1/u.s;

    ocp.lower_bounds  .x(info.ix_length0,           :) = 50.0;
    ocp.upper_bounds  .x(info.ix_length0,           :) = 500.0;
    ocp.variable_units.x(info.ix_length0,           :) = u.m;

    ocp.lower_bounds  .x(info.ix_length0_dot,       :) = -30.0;
    ocp.upper_bounds  .x(info.ix_length0_dot,       :) =  30.0;
    ocp.variable_units.x(info.ix_length0_dot,       :) = u.mps;
    
    
    ocp.lower_bounds  .x(info.ix_angle_of_attack, :) = -5   /180*pi;
    ocp.upper_bounds  .x(info.ix_angle_of_attack, :) = 10   /180*pi;

    ocp.lower_bounds  .u(info.iu_length0_dot_dot,   :) = -20.0;
    ocp.upper_bounds  .u(info.iu_length0_dot_dot,   :) =  20.0;
    ocp.variable_units.u(info.iu_length0_dot_dot,   :) = u.mps2;

    ocp.lower_bounds  .u(info.iu_angle_of_attack_dot, :) = -0.07;
    ocp.upper_bounds  .u(info.iu_angle_of_attack_dot, :) =  0.07;
    ocp.variable_units.u(info.iu_angle_of_attack_dot, :) = 1/u.s;
    
    ocp.lower_bounds  .u(info.iu_roll_angle_dot,      :) = -0.4;
    ocp.upper_bounds  .u(info.iu_roll_angle_dot,      :) =  0.4;
    ocp.variable_units.u(info.iu_roll_angle_dot,      :) = 1/u.s;
    
    ocp.lower_bounds  .u_dot(info.iu_length0_dot_dot,   :) = -10.0;
    ocp.upper_bounds  .u_dot(info.iu_length0_dot_dot,   :) =  10.0;
    ocp.variable_units.u_dot(info.iu_length0_dot_dot,   :) = u.mps2/u.s;
    
    
    ocp.lower_bounds  .u_dot(info.iu_angle_of_attack_dot,   :) = -0.1;
    ocp.upper_bounds  .u_dot(info.iu_angle_of_attack_dot,   :) =  0.1;
    ocp.variable_units.u_dot(info.iu_angle_of_attack_dot,   :) = 1/(u.s^2);
    
    ocp.lower_bounds  .u_dot(info.iu_roll_angle_dot,   :) = -0.5;
    ocp.upper_bounds  .u_dot(info.iu_roll_angle_dot,   :) =  0.5;
    ocp.variable_units.u_dot(info.iu_roll_angle_dot,   :) = 1/(u.s^2);
    
    ocp.variable_units.z(info.iz_theta0_dot_dot,    :) =  1 / u.s^2;
    ocp.variable_units.z(info.iz_phi0_dot_dot,      :) =  1 / u.s^2;
    ocp.variable_units.z(info.iz_theta_dot_dot,    :)  =  1 / u.s^2;
    ocp.variable_units.z(info.iz_phi_dot_dot,      :)  =  1 / u.s^2;

    ocp.lower_bounds  .t(:) =  10.0;
    ocp.upper_bounds  .t(:) =  240;
    ocp.variable_units.t(:) =  u.s;
    
    ocp.lower_bounds  .reelout_phase_fraction(:) =  0.02;
    ocp.upper_bounds  .reelout_phase_fraction(:) =  0.98;
end

