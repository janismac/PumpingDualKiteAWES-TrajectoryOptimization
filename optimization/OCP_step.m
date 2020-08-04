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

function checkpoint_new = OCP_step(ocp, checkpoint, initialization_weight)
    
    variables_scaled = map_struct2(checkpoint.variables, ocp.variable_units, @(v,u)v./u);
    lower_bounds_scaled = map_struct2(ocp.lower_bounds, ocp.variable_units, @(v,u)v./u);
    upper_bounds_scaled = map_struct2(ocp.upper_bounds, ocp.variable_units, @(v,u)v./u);
    parameters_scaled = map_struct2(ocp.parameters, ocp.parameter_units, @(v,u)v./u);

    result = ocp.solver.solve(...
        variables_scaled, ...
        lower_bounds_scaled, ...
        upper_bounds_scaled, ...
        parameters_scaled,...
        initialization_weight,...
        checkpoint.solver_output);
    
    solution_SI = map_struct2(result.solution, ocp.variable_units, @(v,u)v.*u);
    checkpoint_new = struct;
    checkpoint_new.variables = solution_SI;
    checkpoint_new.solver_output = result;
    
end

