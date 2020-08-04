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

function OCP_render(ocp, evaluation)

    close all
    n_frames = floor(evaluation.t(end) * 30);
    t_grid_video = ((1:n_frames)-1) / 30;
    assert(t_grid_video(end) < evaluation.t(end))
        
    vis = Visualization3D(...
        struct('enable_slider_overlay', true));
    vid = MyVideoWriter(false);
    
    fprintf('Rendering Video:\n ');
    lineLength = 0;
    
    frame_values = map_struct(evaluation, @(e) interp1_ND(evaluation.t,e,t_grid_video));
    
    position_center = mean(mean(cat(2,frame_values.position0, frame_values.position), 3), 2);
    
    min_aoa = ocp.lower_bounds.x(ocp.model_info.ix_angle_of_attack1,1)*180/pi;
    max_aoa = ocp.upper_bounds.x(ocp.model_info.ix_angle_of_attack1,1)*180/pi;
    L0_max = ocp.upper_bounds.x(ocp.model_info.ix_length0,1);
    L0_dot_min = ocp.lower_bounds.x(ocp.model_info.ix_length0_dot,1);
    L0_dot_max = ocp.upper_bounds.x(ocp.model_info.ix_length0_dot,1);
    

    patch(vis.ax.axis_overlay_handle,[-0.5 -0.45 -0.45 -0.5],[0.7 0.7 1.3 1.3],aircraft_color(1),'EdgeColor','none')
    sliders(1) = struct('label','$V_{A,1}$',   'format','%.1f m/s',      'min',0,            'max',ocp.parameters.max_airspeed,                'position',[-0.5 0.7 -0.5 1.3],    'side','left',   'data',squeeze(frame_values.aero_speed(:,1,:)));
    sliders(end+1) = struct('label','$\alpha_1$',  'format','%.1f$^\\circ$', 'min',min_aoa,      'max',max_aoa,                                    'position',[-0.45 0.7 -0.45 1.3],  'side','right',  'data',frame_values.angle_of_attack1*180/pi);
    sliders(end+1) = struct('label','$T_0$',       'format','%.f N',         'min',0,            'max',ocp.parameters.max_tension_main_tether,     'position',[0 -1.2 0 0],           'side','left',   'data',frame_values.tether_tension0);
    sliders(end+1) = struct('label','$T_1$',       'format','%.f N',         'min',0,            'max',ocp.parameters.max_tension_secondary_tether,'position',[0 0 -0.5 0.6],         'side','left',   'data',frame_values.tether_tension(1,:));
    sliders(end+1) = struct('label','$\dot{L}_0$',  'format','%.1f m/s',     'min',L0_dot_min,   'max',L0_dot_max,                           'position',[0 -1.2 0 0],           'side','right',  'data',frame_values.length0_dot);

%     if isfield(frame_values, 'max_real_pole')
%         sliders(end+1) = struct('label','$Re(p)_{max}$',  'format','%.1f $$s^{-1}$$',   'min',0,            'max',2,                                    'position',[-1.0 -1.3 -1.0 -0.5],  'side','right',  'data',frame_values.max_real_pole);
%     end
        
    if n_aircraft >= 2
        patch(vis.ax.axis_overlay_handle,[0.5 0.45 0.45 0.5],[0.7 0.7 1.3 1.3],aircraft_color(2),'EdgeColor','none')
        sliders(end+1) = struct('label','$V_{A,2}$',   'format','%.1f m/s',      'min',0,            'max',ocp.parameters.max_airspeed,                'position',[0.5 0.7 0.5 1.3],      'side','right',  'data',squeeze(frame_values.aero_speed(:,2,:)));
        sliders(end+1) = struct('label','$\alpha_2$',  'format','%.1f$^\\circ$', 'min',min_aoa,      'max',max_aoa,                                    'position',[0.45 0.7 0.45 1.3],    'side','left',   'data',frame_values.angle_of_attack2*180/pi);
        sliders(end+1) = struct('label','$T_2$',       'format','%.f N',         'min',0,            'max',ocp.parameters.max_tension_secondary_tether,'position',[0 0 0.5 0.6],          'side','right',  'data',frame_values.tether_tension(2,:));
    end
    
    if n_aircraft >= 3
        patch(vis.axis_overlay_handle,[-0.67 -0.72 -0.72 -0.67],[-0.4 -0.4 0.2 0.2],aircraft_color(3),'EdgeColor','none')
        sliders(end+1) = struct('label','$V_{A,3}$',   'format','%.1f m/s',      'min',0,            'max',ocp.parameters.max_airspeed,                'position',[-0.67 -0.4 -0.67 0.2],      'side','right',  'data',squeeze(frame_values.aero_speed(:,3,:)));
        sliders(end+1) = struct('label','$\alpha_3$',  'format','%.1f$^\\circ$', 'min',min_aoa,      'max',max_aoa,                                    'position',[-0.72 -0.4 -0.72 0.2],    'side','left',   'data',frame_values.angle_of_attack3*180/pi);
        sliders(end+1) = struct('label','$T_3$',       'format','%.f N',         'min',0,            'max',ocp.parameters.max_tension_secondary_tether,'position',[0 0 -0.2 0.6],          'side','right',  'data',frame_values.tether_tension(3,:));
    end
    
    if n_aircraft > 3
        error('not implemented')
    end
    
    
        
    for frame_index = 2:2:n_frames

        cam_angle = 2*pi*frame_index/n_frames;
        %cam_angle = 1;
        camtarget(position_center')
        campos(position_center' + [320*cos(cam_angle) 320*sin(cam_angle) -50])
        camva(70)
        
        frame_value = map_struct(frame_values, @(e) pick_frame(e, frame_index, n_frames));
        
        trail_data = circshift(frame_values.position,-frame_index,3);
        m = 30;
        n = floor(n_frames/m);
        trail_data = trail_data(:, :, (end-m*(n-3)+1):(end-m));
        trail_ticks = circshift(trail_data,frame_index,3);
        trail_ticks = trail_ticks(:, :, 1:m:end);
        
        for i = 1:n_aircraft
            aircraft_pose_matrix(:,:,i) = [frame_value.orientation_matrix(:,:,i) frame_value.position(:,i); 0 0 0 1];
        end
        
        vis.update(...
            'position0',    frame_value.position0, ...
            'aircraft_pose_matrix',        aircraft_pose_matrix, ...
            'trail_data',  trail_data, ...
            'trail_ticks', trail_ticks, ...
            'sliders',      sliders, ...
            'slider_index', frame_index ...
        );
        pause(1e-5)
        drawnow
        vid.writeFrame();
        fprintf(repmat('\b',1,lineLength));
        lineLength = fprintf('%9.2f %%',100*frame_index/n_frames);
    end
    
    fprintf('\nDone\n');
    
end


function e = pick_frame(e, frame_index, n_frames)
    
    D = find(size(e) == n_frames);
    assert(numel(D)==1);
    if D == 1
        e = e(frame_index,:);
    elseif D == 2
        e = e(:,frame_index);
    elseif D == 3
        e = e(:,:,frame_index);
    elseif D == 4
        e = e(:,:,:,frame_index);
    else
        error('');
    end
    
end
