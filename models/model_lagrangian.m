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

function model_lagrangian
    
    coords = generalized_coordinates;
    strains = [coords.strain0; coords.strain.'];
    
    syms mass0 real            % branching point mass
    syms mass_aircraft real    % aircraft mass
    syms tether_density real   % tether linear density
    syms g0 real               % standard gravity
    
    %% Kinetic and potential energy
    v0_sq = coords.velocity0.' * coords.velocity0;
    
    kinetic_energy = 0.5 * v0_sq * (mass0 + tether_density/3 * coords.length0);
    potential_energy = -g0 * (mass0 + tether_density/2 * coords.length0) * coords.position0(3);

    for i = 1:n_aircraft
        vi_sq = coords.velocity(:,i).' * coords.velocity(:,i);
        v0_vi = coords.velocity0.' * coords.velocity(:,i);
        mass_secondary_tether = tether_density * coords.length(i);
        
        kinetic_energy = kinetic_energy + ...
            0.5 * ((mass_aircraft * vi_sq) ...
            + (mass_secondary_tether / 3 * (v0_sq + v0_vi + vi_sq)));

        potential_energy = potential_energy - g0 * (...
              (mass_secondary_tether/2) * coords.position0(3) ...
            + (mass_aircraft + mass_secondary_tether/2) * coords.position(3,i) ...
        );
    end
    
    % Tether tension force placeholders
    tether_tension0 = sym('tether_tension0','real');
    for i = 1:n_aircraft
        tether_tension(i) = sym(['tether_tension' num2str(i)],'real');
    end
    
    % Aerodynamic force placeholders
    aero_force0 = sym({'aero_force0x';'aero_force0y';'aero_force0z'}, 'real');
    for i = 1:n_aircraft
        aero_force(:,i) = sym({['aero_force' num2str(i) 'x'];['aero_force' num2str(i) 'y'];['aero_force' num2str(i) 'z']}, 'real');
    end
    
    % Total applied force
    F0 = -coords.direction0*tether_tension0 + sum(coords.direction .* tether_tension,2) + aero_force0;
    for i = 1:n_aircraft
        F(:,i) = - coords.direction(:,i) * tether_tension(i) + aero_force(:,i);
    end
    
    % Lagrangian mechanics
    q_full = [coords.theta0; coords.phi0; coords.length0; coords.strain0; coords.theta.'; coords.phi.'; coords.length.'; coords.strain.'];
    q_full_dot = total_time_derivative(q_full);
    lagrangian = kinetic_energy - potential_energy;
    dL_dqdot = jacobian(lagrangian, q_full_dot);
    dL_dq = jacobian(lagrangian, q_full);
    ddt_dL_dqdot = total_time_derivative(dL_dqdot);
    
    generalized_forces = F0.' * jacobian(coords.position0, q_full);
    for i = 1:n_aircraft
        generalized_forces = generalized_forces + ...
            F(:,i).' * jacobian(coords.position(:,i), q_full);
    end
    
    implicit_EOM_full = (ddt_dL_dqdot - dL_dq - generalized_forces).';
    
    for variant_cell = {'inelastic', 'elastic'}
        variant = variant_cell{1};

        % Packed parameter vector
        p = [g0; coords.length.'; mass0; mass_aircraft; tether_density; ...
        sym('wind_ref_speed', 'real'); ...
        sym('wind_ref_height', 'real'); ...
        sym('wind_profile_exponent', 'real'); ...
        sym('air_density_msl', 'real'); ...
        sym('scale_height', 'real'); ...
        sym('tether_drag_diameter', 'real'); ...
        sym('reference_area', 'real'); ...
        sym('CL0', 'real'); ...
        sym('CLalpha', 'real'); ...
        sym('CD0', 'real'); ...
        sym('CDalpha', 'real'); ...
        sym('CDalpha2', 'real'); ...
        ];


        if strcmp(variant, 'inelastic')
            %% Select inelastic EOMs
            q = [coords.theta0; coords.phi0; coords.theta.'; coords.phi.'];
            P = sparse(double(jacobian(q_full, q)));
            zero_vars = [tether_tension0; tether_tension.'; strains; total_time_derivative(strains); total_time_derivative(total_time_derivative(strains))];
            implicit_EOM = subs(P' * implicit_EOM_full, zero_vars, zeros(size(zero_vars)));
            z = total_time_derivative(total_time_derivative(q));

            %% Tether force formulas (not part of the ODE, but reqd for winch power calculation)
            % Because jacobian(implicit_EOM_wrt_tethers, [T0 T1 T2]) == eye(3)
            % we can easily solve for the tether forces:
            % F(T,y) == 0 == eye(3) * T + F(0,y)
            % <=>  T == -F(0,y)
            q_tethers = [coords.length0; coords.length.'];
            P_tethers = sparse(double(jacobian(q_full, q_tethers)));
            zero_vars = [strains; total_time_derivative(strains); total_time_derivative(total_time_derivative(strains))];
            implicit_EOM_wrt_tethers = subs(P_tethers' * implicit_EOM_full, zero_vars, zeros(size(zero_vars)));
            J_wrt_tethers = double(simplify(jacobian(implicit_EOM_wrt_tethers, [tether_tension0 tether_tension])));
            assert(all(all(J_wrt_tethers == eye(1+n_aircraft))));
            tether_forces = -subs(implicit_EOM_wrt_tethers, [tether_tension0 tether_tension], zeros(1,1+n_aircraft));
        elseif strcmp(variant, 'elastic')
            q = [coords.theta0; coords.phi0; coords.strain0; coords.theta.'; coords.phi.'; coords.strain.'];
            P = sparse(double(jacobian(q_full, q)));
            z = total_time_derivative(total_time_derivative(q));
            implicit_EOM = P' * implicit_EOM_full;

            % Apply tether force model
            tether_stiffness = sym('tether_stiffness', 'real'); % The EA product of the tether, unit [N]
            tether_damping = sym('tether_damping', 'real'); % The internal tether friction, unit [N*s]
            tether_forces = (tether_stiffness * strains + tether_damping * total_time_derivative(strains)) / 1000;
            implicit_EOM = subs(implicit_EOM, [tether_tension0; tether_tension.'], tether_forces);

            % Tether parameters
            p = [p; tether_stiffness; tether_damping];
        else
            error('not implemented')
        end

        %% Convert the system into the first order ODE form:
        %    M(x,u,p) * z == f(x,u,p)
        %           dx/dt == C * [x; z; u]

        % Packed state, input and disturbace vectors
        wind_disturbance0  = sym({'wind_disturbance0x';'wind_disturbance0y';'wind_disturbance0z'}, 'real');
        for i = 1:n_aircraft
            angle_of_attack(i) = sym(['angle_of_attack' num2str(i)], 'real');
            roll_angle(i)      = sym(['roll_angle' num2str(i)],      'real');
            
            wind_disturbance(:,i)  = sym({['wind_disturbance' num2str(i) 'x'];['wind_disturbance' num2str(i) 'y'];['wind_disturbance' num2str(i) 'z']}, 'real');
        end
        x = [q; total_time_derivative(q); coords.length0; coords.length0_dot; angle_of_attack.'; roll_angle.'];
        u = total_time_derivative([coords.length0_dot; angle_of_attack.'; roll_angle.']);
        w = [wind_disturbance0(:);wind_disturbance(:)];

        % Transform langrangian EOM to mass-matrix form
        M = jacobian(implicit_EOM, z);
        f = -subs(implicit_EOM, z, zeros(size(z)));
        dxdt = total_time_derivative(x);
        C = jacobian(dxdt, [x;z;u]);

        % Store model info
        info = make_index_struct(struct, 'ix_', x);
        info = make_subvector_slice_struct(info, 'ix_', x);
        info = make_index_struct(info, 'iu_', u);
        info = make_subvector_slice_struct(info, 'iu_', u);
        info = make_index_struct(info, 'ip_', p);
        info = make_subvector_slice_struct(info, 'ip_', p);
        info = make_index_struct(info, 'iw_', w);
        info = make_subvector_slice_struct(info, 'iw_', w);
        info = make_index_struct(info, 'iz_', z);
        info = make_subvector_slice_struct(info, 'iz_', z);
        info.x_names = sym_names(x);
        info.u_names = sym_names(u);
        info.p_names = sym_names(p);
        info.w_names = sym_names(w);
        info.z_names = sym_names(z);
        info.n_x = length(x);
        info.n_u = length(u);
        info.n_p = length(p);
        info.n_w = length(w);
        info.n_z = length(z);
        save(['models/model_lagrangian_' variant '_info'], '-struct', 'info');

        %% Secondary model outputs, for OCP constraints and evaulation
        
        % Path curvature
        for i = 1:n_aircraft
            v = coords.velocity(:,i);
            a = total_time_derivative(v);
            curvature(i) = sqrt(sum(cross(v,a).^2)) / ((sum(v.^2)).^(1.5));
        end
        
        % Aircraft proper acceleration
        for i = 1:n_aircraft
            coordinate_acceleration = total_time_derivative(coords.velocity(:,i));
            proper_acceleration(:,i) = coordinate_acceleration - [0; 0; g0];
        end
                
        outputs = struct;
        outputs.tether_tension0 = tether_forces(1);
        outputs.tether_tension = tether_forces(2:end);
        outputs.position = coords.position;
        outputs.position_dot = coords.velocity;
        outputs.position_dot_dot = total_time_derivative(coords.velocity);
        outputs.position0 = coords.position0;
        outputs.position0_dot = coords.velocity0;
        outputs.position0_dot_dot = total_time_derivative(coords.velocity0);
        outputs.direction0 = coords.direction0;
        outputs.direction = coords.direction;
        outputs.curvature = curvature;
        outputs.proper_acceleration = proper_acceleration;
        outputs.kinetic_energy = kinetic_energy;
        outputs.potential_energy = potential_energy;
        
        output_names = fieldnames(outputs);

        % Parameters are constant, eliminate their derivative
        zero_vars = total_time_derivative(p);
        if strcmp(variant, 'inelastic')
            zero_vars = [zero_vars;strains];
        end
        zero_vars = [zero_vars; total_time_derivative(zero_vars); total_time_derivative(total_time_derivative(zero_vars))];
        M = subs(M, zero_vars, zeros(size(zero_vars)));
        f = subs(f, zero_vars, zeros(size(zero_vars)));
        for i = 1:length(output_names)
            outputs.(output_names{i}) = subs(outputs.(output_names{i}), zero_vars, zeros(size(zero_vars)));
            %outputs.(output_names{i}) = simplify(outputs.(output_names{i}));
        end

        %% Generate Matlab code for the model
        M = simplify(M);
        %f = simplify(f);

        matlabFunction(...
            M,f,C, ...
            'File',['models/model_lagrangian_' variant '_generated'], ...
            'Optimize',true, ...
            'Sparse',false, ...
            'Vars', {x,u,p,w,aero_force0,aero_force}, ...
            'Outputs',{'M','f','C'});

        output_expressions = struct2cell(outputs);
        matlabFunction(...
            output_expressions{:}, ...
            'File',['models/model_lagrangian_' variant '_outputs_generated'], ...
            'Optimize',true, ...
            'Sparse',false, ...
            'Vars', {x,u,p,w,z,aero_force0,aero_force}, ...
            'Outputs',output_names);
        
        generate_struct_wrapper_function(...
            ['models/model_lagrangian_' variant '_outputs_wrapper_generated.m'], ...
            output_names, ...
            ['model_lagrangian_' variant '_outputs_generated'], ...
            ['model_lagrangian_' variant '_outputs_wrapper_generated']);
        
        %% Generate inputs for the aerodynamic model
        aero_model_inputs = struct;
        aero_model_inputs.position0              = coords.position0;
        aero_model_inputs.position               = coords.position;
        aero_model_inputs.position0_dot          = coords.velocity0;
        aero_model_inputs.position_dot           = coords.velocity;
        aero_model_inputs.direction0             = coords.direction0;
        aero_model_inputs.direction              = coords.direction;
        aero_model_inputs.length0                = coords.length0;
        aero_model_inputs.angle_of_attack        = angle_of_attack;
        aero_model_inputs.roll_angle             = roll_angle;
        aero_model_inputs.secondary_length       = coords.length;
        aero_model_inputs.wind_ref_speed         = sym('wind_ref_speed',        'real');
        aero_model_inputs.wind_ref_height        = sym('wind_ref_height',       'real');
        aero_model_inputs.wind_profile_exponent  = sym('wind_profile_exponent', 'real');
        aero_model_inputs.wind_disturbance0       = wind_disturbance0;
        aero_model_inputs.wind_disturbance       = wind_disturbance;
        aero_model_inputs.air_density_msl        = sym('air_density_msl',       'real');
        aero_model_inputs.scale_height           = sym('scale_height',          'real');
        aero_model_inputs.tether_drag_diameter   = sym('tether_drag_diameter',  'real');
        aero_model_inputs.reference_area         = sym('reference_area',        'real');
        aero_model_inputs.CL0                    = sym('CL0',                   'real');
        aero_model_inputs.CLalpha                = sym('CLalpha',               'real');
        aero_model_inputs.CD0                    = sym('CD0',                   'real');
        aero_model_inputs.CDalpha                = sym('CDalpha',               'real');
        aero_model_inputs.CDalpha2               = sym('CDalpha2',              'real');
        
        aero_model_input_names = fieldnames(aero_model_inputs);
        
        % Parameters are constant, eliminate their derivative
        zero_vars = total_time_derivative(p);
        if strcmp(variant, 'inelastic')
            zero_vars = [zero_vars;strains];
        end
        zero_vars = [zero_vars; total_time_derivative(zero_vars); total_time_derivative(total_time_derivative(zero_vars))];
        for i = 1:length(aero_model_input_names)
            aero_model_inputs.(aero_model_input_names{i}) = subs(aero_model_inputs.(aero_model_input_names{i}), zero_vars, zeros(size(zero_vars)));
        end

        aero_model_inputs_expressions = struct2cell(aero_model_inputs);
        matlabFunction(...
            aero_model_inputs_expressions{:}, ...
            'File',['models/model_lagrangian_' variant '_aero_inputs_generated'], ...
            'Optimize',true, ...
            'Sparse',false, ...
            'Vars', {x,u,p,w}, ...
            'Outputs',aero_model_input_names);
        
        generate_struct_wrapper_function(...
            ['models/model_lagrangian_' variant '_aero_inputs_wrapper_generated.m'], ...
            aero_model_input_names, ...
            ['model_lagrangian_' variant '_aero_inputs_generated'], ...
            ['model_lagrangian_' variant '_aero_inputs_wrapper_generated']);
        
    end

end
