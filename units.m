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

function u = units
    
    % To avoid numerical problems in optimization, 
    % it is useful to change the system of units. 
    % The conversion constants are defined here. 
    % They are chosen such that typical magnitudes 
    % of various dimensions are logarithmically 
    % close to one.
    
    % Convention: 
    % Multiplication converts from optimization units to SI units.
    %       Division converts from SI units to optimization units.
    
    u = struct;

    u.m = 100;
    u.s = 3;
    u.kg = 30;
    
    u.N = u.kg * u.m / u.s^2;
    u.W = u.kg * u.m^2 / u.s^3;

    u.mps = u.m / u.s;
    u.mps2 = u.m / u.s^2;
    u.mps3 = u.m / u.s^3;
    
    u.kg_m2 = u.kg * u.m^2;
    u.kg_p_m3 = u.kg / u.m^3;
end

