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

classdef NLP_Builder < handle
    
    properties
        variables
        parameters
        constraints
        objective
        objective_terms
    end
    
    methods
        function obj = NLP_Builder()
            obj.variables = struct;
            obj.parameters = struct;
            obj.constraints = struct;
            obj.objective_terms = struct;
            obj.objective = casadi.SX(0);
        end
        
        function var = add_variable(obj, name, varargin)
            assert(~isfield(obj.variables, name), 'Variable exists.');
            var = obj.make_SX(name, varargin{:});
            obj.variables.(name) = var;
        end
        
        function var = add_parameter(obj, name, varargin)
            assert(~isfield(obj.parameters, name), 'Parameter exists.');
            var = obj.make_SX(name, varargin{:});
            obj.parameters.(name) = var;
        end
        
        
        function add_constraint(obj, name, g, min_value, max_value, weight)
            % Adds various kinds of nonlinear constraints.
            % Examples:
            % 
            % >> nlp.add_constraint('foo', g, -10, 5, inf)
            % Hard inequality:  -10 <= g(x) <= 5
            % 
            % 
            % >> nlp.add_constraint('foo', g, -inf, 5, inf)
            % Hard inequality:  g(x) <= 5
            % 
            % 
            % >> nlp.add_constraint('foo', g, 0, 0, inf)
            % Hard equation:   g(x) == 0
            % 
            % 
            % >> nlp.add_constraint('foo', g, -10, 5, 1000)
            % Soft inequality:  
            %       -10 - s <= g(x) <= 5 + s
            %             0 <= s
            %      minimize: (...) + 1000 * s
            % 
            % 
            % >> nlp.add_constraint('foo', g, -inf, 5, 1000)
            % Soft inequality:
            %          g(x) <= 5 + s
            %             0 <= s
            %      minimize: (...) + 1000 * s
            % 
            % 
            % >> nlp.add_constraint('foo', g, 0, 0, 1000)
            % Soft equation:
            %            -s <= g(x) <= s
            %             0 <= s
            %      minimize: (...) + 1000 * s
            % 
            
            assert(ischar(name));
            assert(~isfield(obj.constraints, name), 'Constraint exists.');
            assert(isa(g, 'casadi.SX'));
            assert(numel(min_value) == 1);
            assert(numel(max_value) == 1);
            assert((size(weight, 1) == size(g, 1)) || (size(weight, 1) == 1));
            assert((size(weight, 2) == size(g, 2)) || (size(weight, 2) == 1));
            
            assert(isa(min_value, 'casadi.SX') || isfloat(min_value));
            assert(isa(max_value, 'casadi.SX') || isfloat(max_value));
            assert(isa(weight,    'casadi.SX') || isfloat(weight));
            
            if numel(weight) > 1 && size(weight, 1) == 1
                weight = repmat(weight, size(g, 1), 1);
            end
            
            if numel(weight) > 1 && size(weight, 2) == 1
                weight = repmat(weight, 1, size(g, 2));
            end
            
            is_hard_constraint = (isfloat(weight) && all(weight(:) > 1e20));
            has_lower_bound = ~(isfloat(min_value) && min_value < -1e20);
            has_upper_bound = ~(isfloat(max_value) && max_value > 1e20);
            is_equation = isfloat(min_value) && isfloat(max_value) && min_value == 0.0 && max_value == 0.0;
            
            if isfloat(min_value) && isfloat(max_value)
                assert(max_value >= min_value);
            end
            
            assert(~isfloat(weight) || all(weight(:) > 0));
            
            assert(has_upper_bound, 'Constraint must have an upper bound.');
            obj.constraints.(name) = struct(...
                'g',                   g, ...
                'min_value',           min_value, ...
                'max_value',           max_value, ...
                'weight',              weight, ...
                'is_hard_constraint',  is_hard_constraint, ...
                'has_lower_bound',     has_lower_bound, ...
                'is_equation',         is_equation);
        end
        
        function add_objective(obj, name, term)
            assert(~isfield(obj.objective_terms, name), 'Objective term exists.');
            assert(isa(term, 'casadi.SX'));
            assert(numel(term) == 1);
            obj.objective = obj.objective + term;
            obj.objective_terms.(name) = term;
        end
    end
    
    methods (Access = private)
        function var = make_SX(obj, name, varargin)
            assert(ischar(name));
            if nargin == 4
                if isnumeric(varargin{1}) && isnumeric(varargin{2})
                    m = varargin{1};
                    n = varargin{2};
                    vars = arrayfun(@(i)casadi.SX.sym([name '_' num2str(i)], 1, n(1)+1)',(1:m(1)), 'UniformOutput', false);
                    var = [vars{:}]';
                    var = var(1:end,2:end);
                elseif iscell(varargin{1}) && isnumeric(varargin{2})
                    names = varargin{1};
                    n = varargin{2};
                    vars = arrayfun(@(name)casadi.SX.sym(name{1}, 1, n(1)+1)',names', 'UniformOutput', false);
                    var = [vars{:}]';
                    var = var(1:end,2:end);
                else
                    error('-')
                end
            elseif nargin == 3
                if iscell(varargin{1})
                    names = varargin{1};
                    vars = arrayfun(@(name)casadi.SX.sym(name{1})',names', 'UniformOutput', false);
                    var = [vars{:}]';
                elseif isfloat(varargin{1})
                    var = casadi.SX.sym(name, varargin{1});
                else
                    error('-')
                end
            
            elseif nargin == 2
                var = casadi.SX.sym(name);
            else
                error('-')
            end
        end
    end
end

