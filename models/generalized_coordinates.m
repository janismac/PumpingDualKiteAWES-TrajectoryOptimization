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

function r = generalized_coordinates
    
    r = struct;
    
    r.theta0   = sym('theta0', 'real');
    r.phi0     = sym('phi0', 'real');
    r.length0  = sym('length0', 'real');
    r.strain0  = sym('strain0', 'real');
    
    r.theta0_dot   = total_time_derivative(r.theta0);  
    r.phi0_dot     = total_time_derivative(r.phi0);    
    r.length0_dot  = total_time_derivative(r.length0); 
    r.strain0_dot  = total_time_derivative(r.strain0); 
    
    r.direction0 = [...
         sin(r.theta0) * cos(r.phi0);...
                         sin(r.phi0);...
        -cos(r.theta0) * cos(r.phi0)];
    r.position0 = r.length0 * (1 + r.strain0/1000) * r.direction0;
    r.velocity0 = total_time_derivative(r.position0);
    
    for i = 1:n_aircraft
        i_str = num2str(i);
        
        r.theta(i)   = sym(['theta' i_str], 'real');
        r.phi(i)     = sym(['phi' i_str], 'real');
        r.length(i)  = sym(['length' i_str], 'real');
        r.strain(i)  = sym(['strain' i_str], 'real');
        
        r.theta_dot(i)   = total_time_derivative(r.theta(i));
        r.phi_dot(i)     = total_time_derivative(r.phi(i));
        r.length_dot(i)  = total_time_derivative(r.length(i));
        r.strain_dot(i)  = total_time_derivative(r.strain(i));

        r.direction(:,i) = [...
             sin(r.theta(i)) * cos(r.phi(i));...
                               sin(r.phi(i));...
            -cos(r.theta(i)) * cos(r.phi(i))];

        r.position(:,i) = r.position0 + r.length(i) * (1 + r.strain(i)/1000) * r.direction(:,i);
        r.velocity(:,i) = total_time_derivative(r.position(:,i));
    end
    
    % Generate .m code
    matlabFunction(...
        r.position0, ...
        r.velocity0, ...
        r.direction0, ...
        r.position, ...
        r.velocity, ...
        r.direction, ...
        'File','generalized_coordinates_generated', ...
        'Optimize',false, ...
        'Sparse',false, ...
        'Vars', {
            r.theta0,              r.theta0_dot,   ...
            r.phi0,                r.phi0_dot,     ...
            r.length0,             r.length0_dot,  ...
            r.strain0,             r.strain0_dot,  ...
            r.theta,               r.theta_dot,    ...
            r.phi,                 r.phi_dot,      ...
            r.length,              r.length_dot,   ...
            r.strain,              r.strain_dot,   ...
        }, ...
        'Outputs',{'position0', 'velocity0', 'direction0', 'position', 'velocity', 'direction'});
end
