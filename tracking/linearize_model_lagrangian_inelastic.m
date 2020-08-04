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

function state_space_fn = linearize_model_lagrangian_inelastic
    
    fn = casadi_model_wrapper();
    state_space_fn = @(x, u, z, p) solve_mass_matrix_form(fn, x, u, z, p);
    
end

function [A,B,G] = solve_mass_matrix_form(fn, x, u, z, p)

    % F == 0 == f - M * z
    % dx/dt = Cx * x + Cz * z + Cu * u
    %
    % Implicit function theorem:
    % delta_F == 0 == Fx * delta_x + Fu * delta_u + Fz * delta_z + Fw * delta_w
    % d(delta_x)/dt == Cx * delta_x + Cz * delta_z + Cu * delta_u
    % 
    % -Fz * delta_z == Fx * delta_x + Fu * delta_u + Fw * delta_w
    % delta_z == -inv(Fz) * [ Fx * delta_x + Fu * delta_u + Fw * delta_w ]
    % 
    % 
    % d(delta_x)/dt == Cx * delta_x + Cz * (  -inv(Fz) * [ Fx * delta_x + Fu * delta_u + Fw * delta_w ]  ) + Cu * delta_u
    % 
    % H := -Cz * inv(Fz)
    %
    % d(delta_x)/dt == Cx * delta_x +   H*Fx * delta_x + H*Fu * delta_u + H*Fw * delta_w  + Cu * delta_u
    % d(delta_x)/dt == (Cx + H*Fx) * delta_x       +  (Cu + H*Fu)  * delta_u   +   H*Fw *       delta_w
    %                \ A-state-space /             \  B-state-space /           \ G-state-space /   
    % 

    [C,Fx,Fu,Fz,Fw] = fn(x,u,z,p);
    
    C = sparse(C);
    Fx = sparse(Fx);
    Fu = sparse(Fu);
    Fz = sparse(Fz);
    Fw = sparse(Fw);
    
    
    Cx = C(:, 1:length(x));
    Cz = C(:, length(x) + (1:length(z)));
    Cu = C(:, length(x) + length(z) + (1:length(u)));
    
    H = -Cz/Fz;
    
    A = full(Cx + H*Fx);
    B = full(Cu + H*Fu);
    G = full(H*Fw);
end


function fn = casadi_model_wrapper()
    model_info = load('models/model_lagrangian_inelastic_info');

    x = casadi.SX.sym('x', length(model_info.x_names));
    u = casadi.SX.sym('u', length(model_info.u_names));
    p = casadi.SX.sym('p', length(model_info.p_names));
    z = casadi.SX.sym('z', length(model_info.z_names));
    w = casadi.SX.sym('w', length(model_info.w_names));
    
    aero_inputs = model_lagrangian_inelastic_aero_inputs_wrapper_generated(x,u,p,w);
    aero_outputs = model_aero_forces(aero_inputs);
    [M,f,C] = model_lagrangian_inelastic_generated(x,u,p,w,aero_outputs.applied_force0,aero_outputs.applied_force);
    
    M = M.sparsify();
    f = f.sparsify();
    C = casadi.SX(C).sparsify();
    
    F = f - M*z;
    
    Fx = jacobian(F,x);
    Fu = jacobian(F,u);
    Fz = jacobian(F,z);
    Fw = jacobian(F,w);
    
    Fx = casadi.substitute(Fx,w,0*w);
    Fu = casadi.substitute(Fu,w,0*w);
    Fz = casadi.substitute(Fz,w,0*w);
    Fw = casadi.substitute(Fw,w,0*w);
    
    fn = casadi.Function('model_wrapper',...
        {x,u,z,p},...
        {C,Fx,Fu,Fz,Fw});
end