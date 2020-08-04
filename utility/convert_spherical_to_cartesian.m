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

function [p, p_dot, p_dot_dot] = convert_spherical_to_cartesian(...
    theta, phi, length, ...
    theta_dot, phi_dot, length_dot, ...
    theta_dot_dot, phi_dot_dot, length_dot_dot)
    
    p = [(length.*cos(phi).*sin(theta));(length.*sin(phi));(-length.*cos(phi).*cos(theta))];
    p_dot = [(length_dot.*cos(phi).*sin(theta) - length.*phi_dot.*sin(phi).*sin(theta) + length.*theta_dot.*cos(phi).*cos(theta));(length_dot.*sin(phi) + length.*phi_dot.*cos(phi));(length.*phi_dot.*cos(theta).*sin(phi) - length_dot.*cos(phi).*cos(theta) + length.*theta_dot.*cos(phi).*sin(theta))];
    p_dot_dot = [(- length.*cos(phi).*sin(theta).*phi_dot.^2 - 2.*length.*cos(theta).*sin(phi).*phi_dot.*theta_dot - 2.*length_dot.*sin(phi).*sin(theta).*phi_dot - length.*cos(phi).*sin(theta).*theta_dot.^2 + 2.*length_dot.*cos(phi).*cos(theta).*theta_dot + length_dot_dot.*cos(phi).*sin(theta) - length.*phi_dot_dot.*sin(phi).*sin(theta) + length.*theta_dot_dot.*cos(phi).*cos(theta));(- length.*sin(phi).*phi_dot.^2 + 2.*length_dot.*cos(phi).*phi_dot + length_dot_dot.*sin(phi) + length.*phi_dot_dot.*cos(phi));(length.*cos(phi).*cos(theta).*phi_dot.^2 - 2.*length.*sin(phi).*sin(theta).*phi_dot.*theta_dot + 2.*length_dot.*cos(theta).*sin(phi).*phi_dot + length.*cos(phi).*cos(theta).*theta_dot.^2 + 2.*length_dot.*cos(phi).*sin(theta).*theta_dot - length_dot_dot.*cos(phi).*cos(theta) + length.*phi_dot_dot.*cos(theta).*sin(phi) + length.*theta_dot_dot.*cos(phi).*sin(theta))];
end

