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

function ode_fn = model_lagrangian_inelastic_ode_wrapper
    ode_fn = @(x,u,p,w) ode_wrapper(x,u,p,w);
end

function dxdt = ode_wrapper(x,u,p,w)
    aero_inputs = model_lagrangian_inelastic_aero_inputs_wrapper_generated(x,u,p,w);
    aero_outputs = model_aero_forces(aero_inputs);
    [M,f,C] = model_lagrangian_inelastic_generated(x,u,p,w,aero_outputs.applied_force0,aero_outputs.applied_force);
    assert(max(abs(f)) < 1e14);
    dxdt = C * [x; (M\f); u];
end