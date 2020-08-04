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

function N = n_aircraft
    % This function acts as a global constant for the number of aircraft
    % and the number of secondary tethers. Allowed values are 1,2,3,...
    % Expect very poor performance for larger values.
    % This number is assumed to be constant within a working directory.
    % The models need to be re-generated after it changes.
    N = 2;
end
