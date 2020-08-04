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

function [result,summary] = tracking_validation_lagrangian_inelastic(evaluation, p, u_max, n_loops)

    result = struct;
    
    x_ref = evaluation.x;
    u_ref = evaluation.u;
    t_ref = evaluation.t;
    K     = evaluation.K;
    
    dt = median(diff(t_ref));
    assert(dt < 0.1);
    assert(std(diff(t_ref)) < 1e-13); % has equidistant grid?
    
    % Expand loops / periods
    x_ref = repmat(x_ref(:,2:end), 1, n_loops);
    u_ref = repmat(u_ref(:,2:end), 1, n_loops);
    K     = repmat(K(:,:,2:end), 1, 1, n_loops);
    t_ref = (0:(size(x_ref,2)-1)) * median(diff(t_ref));
    
    
    % Solve ODE with tracking controller
    x0 = x_ref(:, 1);
    x0 = x0 + 1e-6; % Add negligble disturbance
    
    ode_fn = model_lagrangian_inelastic_ode_wrapper;
    info = load('models/model_lagrangian_inelastic_info');
    w = zeros(length(info.w_names), 1);
    
    my_ode = @(t,x) closed_loop_ode_wrapper(...
        t, dt, x, ...
        x_ref,  ...
        u_ref,  ...
        p,w, K, u_max, ode_fn);

    [~,YOUT] = ode45(my_ode, t_ref, x0, odeset('RelTol',1e-4,'AbsTol',1e-4));
    x = YOUT';
    
    
    % Evaluate tracking errors
    pos = convert_lagrangian_inelastic_state_to_cartesian(x,p);
    pos_ref = convert_lagrangian_inelastic_state_to_cartesian(x_ref,p);
    
    abs_errors = map_struct2(pos, pos_ref, @(a,b)a-b);
    
    position0_dist = vecnorm(abs_errors.position0);
    position_dist = vecnorm(abs_errors.position,1,2);
    
    result.x = x;
    result.x_ref = x_ref;
    result.u_ref = u_ref;
    result.t_ref = t_ref;
    result.abs_errors = abs_errors;
    result.position0_dist = position0_dist;
    result.position_dist = position_dist;
    
    % plot(t_ref, [position0_dist; position1_dist; position2_dist])
    
    distances = [position0_dist(:); position_dist(:)];
    
    summary = struct;
    summary.max_distance_deviation = max(distances);
    summary.median_distance_deviation = median(distances);
end

function i=time_index(t,dt,n)
    i = floor(t/dt)+1;
    i = max(1,i);
    i = min(n,i);
end

function dxdt = closed_loop_ode_wrapper(t, dt_ref, x, x_ref, u_ref, p,w, K, u_max, ode_fn)

    I = time_index(t,dt_ref,size(x_ref, 2)-1);
    
    alpha = mod(t, dt_ref)/dt_ref;
    
    x_ref_now = x_ref(:,I) + alpha * (x_ref(:,I+1) - x_ref(:,I));
    u_ref_now = u_ref(:,I) + alpha * (u_ref(:,I+1) - u_ref(:,I));
    K_now     = K(:,:,I)   + alpha * (K(:,:,I+1)   - K(:,:,I));

    u = u_ref_now - K_now * (x - x_ref_now);
    
    u = min(u,  u_max);
    u = max(u, -u_max);

    dxdt = ode_fn(x,u,p,w);
end


