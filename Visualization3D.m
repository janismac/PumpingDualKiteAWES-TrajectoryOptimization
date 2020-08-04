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

classdef Visualization3D < handle
    
    properties
        ax
        aircraft_patch
        aircraft_transform
        sphere_surf
        sphere_transform
        tether0_surf
        tether0_transform
        tether_surf
        tether_transform
        
        trail_plot
        trail_tick_plot
        
        reference_mode
        
        slider_handles
        slider_triangles
        
        scale_factor
    end
    
    methods
        function obj = Visualization3D(inputs)
            
            if isfield(inputs, 'axes')
                obj.ax = inputs.axes;
            else
                obj.ax = Visualization3D_axes(inputs.enable_slider_overlay);
            end
            
            obj.reference_mode = false;
            if isfield(inputs, 'reference_mode')
                obj.reference_mode = inputs.reference_mode;
            end
            
            obj.scale_factor = 8;
            aircraft_3d_model = load('Glider_3D_model/files/Glider.mat');
            
            aircraft_color_modifier = @(c) c;
            if obj.reference_mode
                aircraft_color_modifier = @(c) (0.3 * c);
            end

            for i = 1:n_aircraft
                obj.aircraft_patch{i} = patch(...
                    'Faces',           aircraft_3d_model.faces, ...
                    'Vertices',        aircraft_3d_model.vertices, ...
                    'FaceColor',       aircraft_color_modifier(aircraft_color(i)*0.999), ...
                    'EdgeColor',       'none', ...
                    'FaceLighting',    'gouraud', ...
                    'SpecularStrength',   0.0, ...
                    'AmbientStrength', 1);

                obj.aircraft_transform{i} = hgtransform('Parent', obj.ax.axis_handle);
                obj.aircraft_patch{i}.Parent = obj.aircraft_transform{i};
                obj.aircraft_transform{i}.Matrix = [eye(3) [0 0 -100]';[0 0 0 1]];
            end
            
            [sphere_x,sphere_y,sphere_z] = sphere(32);
            obj.sphere_surf = surf(...
                sphere_x,sphere_y,sphere_z, ...
                'FaceColor',       [1 1 1]*0.999, ...
                'EdgeColor',       'none', ...
                'FaceLighting',    'gouraud', ...
                'SpecularStrength',   0.0, ...
                'AmbientStrength', 1);
            obj.sphere_transform = hgtransform('Parent', obj.ax.axis_handle);
            obj.sphere_surf.Parent = obj.sphere_transform;
            
            
            [cylinder_X,cylinder_Y,cylinder_Z] = cylinder;
            
            obj.tether0_surf = surf(...
                cylinder_X,cylinder_Y,cylinder_Z, ...
                'FaceColor',       [1 1 1]*0.999, ...
                'EdgeColor',       'none', ...
                'FaceLighting',    'gouraud', ...
                'SpecularStrength',   0.0, ...
                'AmbientStrength', 1);
            obj.tether0_transform = hgtransform('Parent', obj.ax.axis_handle);
            obj.tether0_surf.Parent = obj.tether0_transform;
            
            for i = 1:n_aircraft
                obj.tether_surf{i} = surf(...
                    cylinder_X,cylinder_Y,cylinder_Z, ...
                    'FaceColor',       [1 1 1]*0.999, ...
                    'EdgeColor',       'none', ...
                    'FaceLighting',    'gouraud', ...
                    'SpecularStrength',   0.0, ...
                    'AmbientStrength', 1);
                obj.tether_transform{i} = hgtransform('Parent', obj.ax.axis_handle);
                obj.tether_surf{i}.Parent = obj.tether_transform{i};
            end
            
            
            % Cube for the winch station
            cube = [1 -1 -1 1;1 1 -1 -1; 1 1 1 1];
            cube = cat(3,cube, circshift(cube,1,1), circshift(cube,2,1));
            cube = cat(3,cube, -cube);
            patch(...
                'XData',squeeze(cube(1,:,:)) * 2 * obj.scale_factor,...
                'YData',squeeze(cube(2,:,:)) * 2 * obj.scale_factor,...
                'ZData',squeeze(cube(3,:,:)) * obj.scale_factor, ...
                'FaceColor',       [1 1 1]*0.999, ...
                'EdgeColor',       'none', ...
                'FaceLighting',    'gouraud', ...
                'SpecularStrength',   0.0, ...
                'AmbientStrength', 1);
            
            % Aircraft trails / streamers
            for i = 1:n_aircraft
                obj.trail_plot{i} = plot3(0,0,0,'Color',aircraft_color_modifier(aircraft_color(i)),'LineWidth',1.4);
                obj.trail_tick_plot{i} = plot3(0,0,0,'o','MarkerFaceColor',aircraft_color_modifier(aircraft_color(i)),'MarkerEdgeColor','none','MarkerSize',4);
            end
            
            
            if inputs.enable_slider_overlay
                obj.slider_triangles = patch(obj.ax.axis_overlay_handle, 0, 0, 0);
                obj.slider_triangles.EdgeColor = 'none';
                if obj.reference_mode
                    obj.slider_triangles.FaceColor = [1 1 1]*0.5;
                else
                    obj.slider_triangles.FaceColor = 'k';
                end
            end
                
            
            % Hide axis lines and ticks
            obj.ax.axis_handle.XTickLabel = [];
            obj.ax.axis_handle.YTickLabel = [];
            obj.ax.axis_handle.ZTickLabel = [];
            
            obj.ax.axis_handle.XAxis.Visible = 'off';
            obj.ax.axis_handle.YAxis.Visible = 'off';
            obj.ax.axis_handle.ZAxis.Visible = 'off';
            
            obj.ax.axis_handle.XColor = 'none';
            obj.ax.axis_handle.YColor = 'none';
            obj.ax.axis_handle.ZColor = 'none';
            
            obj.ax.axis_handle.GridColor = [0 0 0];
            obj.ax.axis_handle.GridAlpha = 0.5;
            
            xlabel('')
            ylabel('')
            zlabel('')
            
            xlim([-50 600])
            ylim([-600 600])
            zlim([-600 100])
            
            campos([20 20 -120])
            camtarget([0 0 -100])
            camva(40)
            
            obj.ax.axis_handle.Clipping = 'off';
            grid off
            
            % Disable toolbars, since they are also visible in the video
            obj.ax.axis_handle.Toolbar.Visible = 'off';
        end
        
        function update(obj, varargin)
            warning('off','MATLAB:hg:DiceyTransformMatrix');
            
            assert(mod(length(varargin),2) == 0);
            inputs = struct;
            for i = 1:2:length(varargin)
                inputs.(varargin{i}) = varargin{i+1};
            end
            
            obj.sphere_transform.Matrix = [eye(3)*0.2*obj.scale_factor inputs.position0; 0 0 0 1];
            obj.tether0_transform.Matrix = cylinder_pose_matrix(obj,inputs.position0,[0;0;0]);
            
            for i = 1:n_aircraft
                obj.aircraft_transform{i}.Matrix = ...
                    inputs.aircraft_pose_matrix(:,:,i) ...
                    * diag([[1 1 1]*obj.scale_factor 1]);

                obj.tether_transform{i}.Matrix = cylinder_pose_matrix(obj,inputs.position0,inputs.aircraft_pose_matrix(1:3,4,i));

                obj.trail_plot{i}.XData = inputs.trail_data(1,i,:);
                obj.trail_plot{i}.YData = inputs.trail_data(2,i,:);
                obj.trail_plot{i}.ZData = inputs.trail_data(3,i,:);

                obj.trail_tick_plot{i}.XData = inputs.trail_ticks(1,i,:);
                obj.trail_tick_plot{i}.YData = inputs.trail_ticks(2,i,:);
                obj.trail_tick_plot{i}.ZData = inputs.trail_ticks(3,i,:);
            end
            
            
            
            obj.ax.axis_handle.CameraUpVector = [0 0 -1];
            
            %% Sliders
            if isfield(inputs, 'sliders')
                if isempty(obj.slider_handles)
                    for i = 1:length(inputs.sliders)
                        h_align = 'left';
                        if strcmp(inputs.sliders(i).side, 'left')
                            h_align = 'right';
                        end
                        obj.slider_handles(i).label = text(obj.ax.axis_overlay_handle, 0, 0, 'hello','Interpreter','Latex','FontSize',14,'HorizontalAlignment',h_align);

                        % Slider lines
                        plot(obj.ax.axis_overlay_handle,...
                            inputs.sliders(i).position([1 3]), inputs.sliders(i).position([2 4]),...
                            'k','LineWidth',2);
                    end
                end

                slider_triangle_data = zeros(2,3,length(inputs.sliders));

                for i = 1:length(inputs.sliders)
                    p0 = inputs.sliders(i).position(1:2)';
                    p1 = inputs.sliders(i).position(3:4)';
                    L = norm(p1-p0);
                    D = (p1-p0)./L;
                    N = [0 -1;1 0] * D;
                    if ~strcmp(inputs.sliders(i).side, 'left')
                        N = -N;
                    end
                    value = inputs.sliders(i).data(inputs.slider_index);

                    max_value = inputs.sliders(i).max;
                    min_value = inputs.sliders(i).min;

                    relative_position = (value-min_value)/(max_value-min_value);
                    abs_position = p0 + relative_position * (p1-p0);

                    
                    if obj.reference_mode
                        obj.slider_handles(i).label.String = '';
                        f = 0.5;
                    else
                        value_str = sprintf(inputs.sliders(i).format, value);
                        label_str = [inputs.sliders(i).label ': ' value_str];
                        obj.slider_handles(i).label.String = label_str;
                        f = 1;
                    end
                    obj.slider_handles(i).label.Position = (abs_position + 0.1 * N)';
                    %obj.slider_handles(i).label = text(obj.ax.axis_overlay_handle, 0, 0, 'hello');
                    
                    
                    slider_triangle_data(:,:,i) = abs_position + 0.02 * [0*D f*D+2.5*N -f*D+2.5*N] + 0.01 * N;
                end

                obj.slider_triangles.XData = squeeze(slider_triangle_data(1,:,:));
                obj.slider_triangles.YData = squeeze(slider_triangle_data(2,:,:));
            end
        end
        
        function P = cylinder_pose_matrix(obj,p1,p2)
            ez = p2 - p1;
            dist = norm(ez);
            ez = ez / dist;
            ex = [0 0 1]' - (([0 0 1]*ez)*ez);
            ex = ex / norm(ex);
            ey = cross(ez,ex);
            R = [ex ey ez];
            R = R * diag([[1 1]*0.05*obj.scale_factor dist]);
            P = [R p1; 0 0 0 1];
        end
    end
end

