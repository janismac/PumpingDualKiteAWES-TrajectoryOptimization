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

function c = aircraft_color(i)

    if i == 1
        c = [0.3 1 0.3];
    elseif i == 2
        c = [1 0.3 0.3];
    elseif i == 3
        c = [0.5 0.5 1];
    else
        error('not implemented')
    end

end

