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

function main
    
    %% Generate model
    if exist('model_lagrangian_inelastic_generated') ~= 2 % Needs to run only once
        % Warning: The model generation takes a very long time (>20 min) in
        % some verions of Matlab due to an issue with matlabFunction().
        % Use Matlab R2019a (No Update) or Matlab R2019a (Update 3) for
        % faster results.
        model_lagrangian(); 
    end
    
    %% Run optimizaiton
    config = struct;
    config.solver                          = 'ipopt_mumps';
    config.ipopt_nlp_scaling_method        = 'gradient-based';
    config.ipopt_bound_push                = 1e-06;
    config.ipopt_mu_max                    = 0.01;
    config.ipopt_tol                       = 1e-05;
    config.wind_speed                      = 11;
    config.tether_diameter                 = 0.003;
    config.secondary_tether_length         = 150;
    config.collocation_n_nodes             = 5;
    config.collocation_n_timesteps         = 100;
    config.ode_weight                      = 1000000;
    config.path_constraint_weight          = 100000;
    config.theta0_upper_bound              = 70/180*pi;
    config.max_curvature                   = 1/30;
    config.initialization_period_shift     = 0.15;
    config.minimum_aircraft_distance       = 60;
    config.n_revolutions                   = 3;
    config.output_directory                = 'output/sweep_test';
    OCP_runall(config);
    
    
    %% Render trajectory animation
    trajectory = load('output/sweep_test/OCP_results2.mat');
    OCP_render(trajectory.ocp2_, trajectory.evaluation2);
    
    
    %% Run parameter a perturbation sweep
    mkdir('output/sweep_perturbation_base');
    copyfile('output/sweep_test/OCP_results2.mat', 'output/sweep_perturbation_base/OCP_results2.mat');
    OCP_solution_perturbation('wind_speed');


    %% Render perturbation animation
    perturbation_output = load('output/sweep_perturbation_base/OCP_wind_speed_perturbation.mat');
    OCP_render_perturbation(perturbation_output);
    
end