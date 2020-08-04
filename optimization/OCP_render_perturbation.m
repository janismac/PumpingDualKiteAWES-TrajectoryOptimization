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

function OCP_render_perturbation(perturbation_output)

    % Find camera target
    N_perturbation = length(perturbation_output.perturbation_values);
    for i = 1:N_perturbation
        P = cat(2, perturbation_output.perturbation_evaluation(i).position0, perturbation_output.perturbation_evaluation(i).position);
        mean_pos(:,i) = mean(mean(P, 2), 3);
    end
    cam_target = mean(mean_pos,2);
    
    % Animation sequence
    perturbation_indices = round(1:0.5:N_perturbation);
    perturbation_indices = [perturbation_indices fliplr(perturbation_indices)];
    perturbation_indices = repmat(perturbation_indices, 1, 8);
    n_frames = length(perturbation_indices);
    
    vis = Visualization3D(...
        struct('enable_slider_overlay', false));
    vid = MyVideoWriter(false);
    lineLength = 0;
    
    text_perturbation_value = annotation(...
        'textarrow',[0.7 0.1],[0.95 0.95],...
        'String','value: 123',...
        'FontSize',18,...
        'HeadLength',0,...
        'TextColor',[1 1 1]*0.999,...
        'HeadWidth',0,...
        'LineStyle','none',...
        'Interpreter','none');
    
    for j = 1:n_frames
        i_perturbation = perturbation_indices(j);

        evaluation = perturbation_output.perturbation_evaluation(i_perturbation);
        frame_value = map_struct(evaluation, @(e) pick_frame(e, 1, length(evaluation.t)));

        trail_data = evaluation.position;
        trail_ticks = interp1_ND(evaluation.t,trail_data,0:1:evaluation.t(end));

        for i = 1:n_aircraft
            aircraft_pose_matrix(:,:,i) = [frame_value.orientation_matrix(:,:,i) frame_value.position(:,i); 0 0 0 1];
        end

        vis.update(...
            'position0',    frame_value.position0, ...
            'aircraft_pose_matrix',        aircraft_pose_matrix, ...
            'trail_data',  trail_data, ...
            'trail_ticks', trail_ticks ...
        );
    
        text_perturbation_value.String = ...
            [perturbation_output.perturbation_config.parameter_name ': ' ...
            sprintf('%6.1f', perturbation_output.perturbation_values(i_perturbation))];
    
    
        camtarget(cam_target);
        cam_angle = 2*pi*j/n_frames;
        campos(cam_target' + [390*cos(cam_angle) 390*sin(cam_angle) -50])
        camva(70)
        pause(1e-5)
        drawnow
        vid.writeFrame();
        fprintf(repmat('\b',1,lineLength));
        lineLength = fprintf('%9.2f %%',100*j/n_frames);
    end


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
