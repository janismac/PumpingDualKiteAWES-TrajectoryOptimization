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

function initialization = OCP_initialize_multiple_revolutions(ocp, init_evaluation, n_revolutions)

    T_init = max(init_evaluation.t) * n_revolutions;

    init_trajectory = struct('t',init_evaluation.t,  'u',init_evaluation.u,  'u_dot',init_evaluation.u_dot,  'x',init_evaluation.x,  'z',init_evaluation.z);

    % Create multiple revolutions by copying the periodic trajectory
    init_trajectory = map_struct(init_trajectory, @(e) [repmat(e(:,1:end-1),1,n_revolutions-1)  e]);

    % Fix the timeline
    dt = diff(init_trajectory.t);
    dt(dt < 0) = median(dt);
    init_trajectory.t = [0 cumsum(dt)];
    assert(abs(init_trajectory.t(end) - T_init) < 1e-10);

    % Evaluate on new collocation nodes
    reelout_phase_fraction = 1 - 1/(n_revolutions+1);
    interval_durations = full(ocp.interval_duration_fn(T_init, reelout_phase_fraction));
    interval_duration_sum = cumsum([0; interval_durations]);
    init_t_grid = interval_duration_sum(1:end-1)' + interval_durations' .* ocp.collocation.nodes(1:end-1);
    init_t_grid = init_t_grid(:);
    init_trajectory = map_struct(init_trajectory, @(e) interp1(init_trajectory.t,e',init_t_grid)');

    
    initialization = struct;
    
    if isfield(ocp.lower_bounds, 'phase_switch_speed')
        initialization.phase_switch_speed = 0;
    end
    
    initialization.reelout_phase_fraction = reelout_phase_fraction;
    initialization.t = T_init;
    initialization.u = init_trajectory.u;
    initialization.u_dot = init_trajectory.u_dot;
    initialization.x = init_trajectory.x;
    initialization.z = init_trajectory.z;

    initialization = map_struct2(initialization, ocp.lower_bounds, @(I,B) max(I,B));
    initialization = map_struct2(initialization, ocp.upper_bounds, @(I,B) min(I,B));
end
