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

function [dF, vars, vars_dot] = total_time_derivative(F)
    %% Find total derivative with resprect to time.
    % The time derivatives of symbols are automatically created and
    % post-fixed with '_dot'.
    % Example:
    %     syms x real
    %     total_time_derivative(x^2 + sin(x))
    %     Result: 2 * x * x_dot + x_dot * cos(x)


    time_var = sym('time_var_19e3ccff09de09d5703','real'); % arbitrary, unique name, to avoid name collisions
    vars = symvar(F);
    vars_dot = 0*vars;

    for i = 1:length(vars)
        vars_dot(i) = sym([char(vars(i)) '_dot'],'real');
    end
    
    F = subs(F, vars, vars + time_var * vars_dot);
    dF = diff(F, time_var);
    dF = subs(dF, time_var, 0);
end

