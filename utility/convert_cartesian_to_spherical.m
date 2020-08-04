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

function [length, theta, phi, length_dot, theta_dot, phi_dot, length_dot_dot, theta_dot_dot, phi_dot_dot] = ...
    convert_cartesian_to_spherical(p,p_dot,p_dot_dot)
    
%     x  = sym('x', 'real');
%     y  = sym('y', 'real');
%     z  = sym('z', 'real');
%         
%     length = sqrt(sum([x y z].^2));
%     theta = atan2(x, -z);
%     phi = asin(y/length);
%      
%     length_dot = total_time_derivative(length);
%     theta_dot = total_time_derivative(theta);
%     phi_dot = total_time_derivative(phi);
% 
%     length_dot_dot = total_time_derivative(length_dot);
%     theta_dot_dot = total_time_derivative(theta_dot);
%     phi_dot_dot = total_time_derivative(phi_dot);
% 
%     
%     fn = matlabFunction(length, theta, phi, length_dot, theta_dot, phi_dot, length_dot_dot, theta_dot_dot, phi_dot_dot)
    
    x = p(1);
    y = p(2);
    z = p(3);
    
    x_dot = p_dot(1);
    y_dot = p_dot(2);
    z_dot = p_dot(3);
    
    x_dot_dot = p_dot_dot(1);
    y_dot_dot = p_dot_dot(2);
    z_dot_dot = p_dot_dot(3);

    [length, theta, phi, length_dot, theta_dot, phi_dot, length_dot_dot, theta_dot_dot, phi_dot_dot] = deal(sqrt(x.^2+y.^2+z.^2),atan2(x,-z),asin(y.*1.0./sqrt(x.^2+y.^2+z.^2)),((x.*x_dot.*2.0+y.*y_dot.*2.0+z.*z_dot.*2.0).*1.0./sqrt(x.^2+y.^2+z.^2))./2.0,-(z.^2.*(x_dot./z-x.*1.0./z.^2.*z_dot))./(x.^2+z.^2),1.0./sqrt(-y.^2./(x.^2+y.^2+z.^2)+1.0).*(y_dot.*1.0./sqrt(x.^2+y.^2+z.^2)-(y.*(x.*x_dot.*2.0+y.*y_dot.*2.0+z.*z_dot.*2.0).*1.0./(x.^2+y.^2+z.^2).^(3.0./2.0))./2.0),(x.*x_dot.*2.0+y.*y_dot.*2.0+z.*z_dot.*2.0).^2.*1.0./(x.^2+y.^2+z.^2).^(3.0./2.0).*(-1.0./4.0)+(1.0./sqrt(x.^2+y.^2+z.^2).*(x.*x_dot_dot.*2.0+y.*y_dot_dot.*2.0+z.*z_dot_dot.*2.0+x_dot.^2.*2.0+y_dot.^2.*2.0+z_dot.^2.*2.0))./2.0,-(z.^2.*(x_dot_dot./z-x.*1.0./z.^2.*z_dot_dot-x_dot.*1.0./z.^2.*z_dot.*2.0+x.*1.0./z.^3.*z_dot.^2.*2.0))./(x.^2+z.^2)+z.^2.*(x.*x_dot.*2.0+z.*z_dot.*2.0).*1.0./(x.^2+z.^2).^2.*(x_dot./z-x.*1.0./z.^2.*z_dot)-(z.*z_dot.*(x_dot./z-x.*1.0./z.^2.*z_dot).*2.0)./(x.^2+z.^2),1.0./sqrt(-y.^2./(x.^2+y.^2+z.^2)+1.0).*(y_dot_dot.*1.0./sqrt(x.^2+y.^2+z.^2)-(y.*1.0./(x.^2+y.^2+z.^2).^(3.0./2.0).*(x.*x_dot_dot.*2.0+y.*y_dot_dot.*2.0+z.*z_dot_dot.*2.0+x_dot.^2.*2.0+y_dot.^2.*2.0+z_dot.^2.*2.0))./2.0-y_dot.*(x.*x_dot.*2.0+y.*y_dot.*2.0+z.*z_dot.*2.0).*1.0./(x.^2+y.^2+z.^2).^(3.0./2.0)+y.*(x.*x_dot.*2.0+y.*y_dot.*2.0+z.*z_dot.*2.0).^2.*1.0./(x.^2+y.^2+z.^2).^(5.0./2.0).*(3.0./4.0))+(((y.*y_dot.*2.0)./(x.^2+y.^2+z.^2)-y.^2.*(x.*x_dot.*2.0+y.*y_dot.*2.0+z.*z_dot.*2.0).*1.0./(x.^2+y.^2+z.^2).^2).*1.0./(-y.^2./(x.^2+y.^2+z.^2)+1.0).^(3.0./2.0).*(y_dot.*1.0./sqrt(x.^2+y.^2+z.^2)-(y.*(x.*x_dot.*2.0+y.*y_dot.*2.0+z.*z_dot.*2.0).*1.0./(x.^2+y.^2+z.^2).^(3.0./2.0))./2.0))./2.0);
    
end

