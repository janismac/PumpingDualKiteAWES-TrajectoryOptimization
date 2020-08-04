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

function R = rotation_matrix_angle_axis(theta, u)

    ux = u(1);
    uy = u(2);
    uz = u(3);
    
    c = cos(theta);
    s = sin(theta);

    R = [...
    (c + ux^2 * (1-c)        ), (ux * uy * (1-c) - uz * s), (ux * uz * (1-c) + uy * s); ...
    (uy * ux * (1-c) + uz * s), (c + uy^2 * (1-c)        ), (uy * uz * (1-c) - ux * s); ...
    (uz * ux * (1-c) - uy * s), (uz * uy * (1-c) + ux * s), (c + uz^2 * (1-c)        )];

end
