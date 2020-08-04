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

function max_position_error_3sec = OCP_validate_openloop(ocp, evaluation)
    u_trajectory = @(t)(interp1(evaluation.t,evaluation.u',t)');
    
    ode_fn = model_lagrangian_inelastic_ode_wrapper;
    w = zeros(length(ocp.model_info.w_names), 1);
    my_ode = @(t,x) ode_fn(x,u_trajectory(t),ocp.parameters.p,w);
    
    dt = median(diff(evaluation.t));
    
    n_samples = length(evaluation.t);
    n_slice = round(3/dt); % Integrate over 3 second intervals
    
    max_position_error_3sec = 0;
    
    for i = 1:n_slice:n_samples
        slice = (0:n_slice) + i;
        slice = slice(slice <= n_samples);
        
        if length(slice) < 3
            break
        end

        [~,YOUT] = ode45(my_ode,evaluation.t(slice),evaluation.x(:,i),odeset('RelTol',1e-5,'AbsTol',1e-5));
        YOUT = YOUT';
        
        cartesian_ref = convert_lagrangian_inelastic_state_to_cartesian(evaluation.x(:,slice),ocp.parameters.p);
        cartesian_soln = convert_lagrangian_inelastic_state_to_cartesian(YOUT,ocp.parameters.p);
        
        errors = map_struct2(cartesian_ref, cartesian_soln, @(a,b) vecnorm(a-b));
        
        max_position_error_3sec = max(max_position_error_3sec, max(max(max([errors.position0 errors.position]))));
    end
end

