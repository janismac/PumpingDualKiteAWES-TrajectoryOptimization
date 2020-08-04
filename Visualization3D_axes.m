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

classdef Visualization3D_axes < handle
    
    properties
        figure_handle
        axis_handle
        axis_overlay_handle
    end
    
    methods
        function obj = Visualization3D_axes(enable_slider_overlay)
            
            close all
            setup_3D_plot();
            
            obj.figure_handle = gcf;
            
            obj.figure_handle.Units = 'pixel';
            drawnow 
            obj.figure_handle.Position = [100 100 1280 720];
            
            obj.axis_handle = gca;
            
            if enable_slider_overlay
                obj.axis_overlay_handle = axes;

                obj.figure_handle.CurrentAxes = obj.axis_overlay_handle;
                obj.axis_overlay_handle.Color = [1 1 1]*0.999;
                split = 0.55;
                margin = 0.01;
                obj.axis_overlay_handle.Position = [(split+margin) margin*16/9 (1-split-2*margin) 1-2*margin*16/9];
                obj.axis_overlay_handle.XLim = [-1 1]*1.1;
                obj.axis_overlay_handle.Visible = 'on';
                obj.axis_overlay_handle.XTickLabel = [];
                obj.axis_overlay_handle.YTickLabel = [];
                obj.axis_overlay_handle.ZTickLabel = [];
                obj.axis_overlay_handle.XAxis.Visible = 'off';
                obj.axis_overlay_handle.YAxis.Visible = 'off';
                obj.axis_overlay_handle.ZAxis.Visible = 'off';
                hold(obj.axis_overlay_handle,'on');
                axis(obj.axis_overlay_handle,'equal')
                obj.axis_overlay_handle.YLim = diff(obj.axis_overlay_handle.YLim)/2 * [-1 1];
                
                obj.axis_overlay_handle.Toolbar.Visible = 'off';
                
                obj.figure_handle.CurrentAxes = obj.axis_handle;
                obj.axis_handle.Position = [0 0 split 1];
            end
            

            light('Position',[0 0 -1])
            light('Position',[0 0 -1],'Color',[1 1 1]*0.4)

            obj.figure_handle.Color = [77 113 175]/255;
            obj.figure_handle.InvertHardcopy = 'off';
            obj.axis_handle.Color = [77 113 175]/255;
            obj.axis_handle.Visible = 'off';
            obj.axis_handle.AmbientLightColor = 0.2*[1 1 1];

            
                
            % Green checkerboard for the ground
            [gird_x,gird_y]=meshgrid(-5:5,-5:5);
            gird_x = gird_x(:);
            gird_y = gird_y(:);
            gird_x = gird_x + [0 0 0.5 0.5];
            gird_y = gird_y + [0 0.5 0.5 0];
            gird_x = [gird_x;gird_x+0.5]';
            gird_y = [gird_y;gird_y+0.5]';
            
            patch(...
                'XData',1000*gird_x,...
                'YData',1000*gird_y,...
                'ZData',0*gird_y, ...
                'FaceColor',       [69 92 64]/255, ...
                'EdgeColor',       'none', ...
                'FaceLighting',    'gouraud', ...
                'SpecularStrength',   0.0, ...
                'AmbientStrength', 1);
            
            patch(...
                'XData',1000*(gird_x+0.5),...
                'YData',1000*gird_y,...
                'ZData',0*gird_y, ...
                'FaceColor',       [121 126 96]/255, ...
                'EdgeColor',       'none', ...
                'FaceLighting',    'gouraud', ...
                'SpecularStrength',   0.0, ...
                'AmbientStrength', 1);
            
            
            % Ground texture
            % texture_size_in_meters = 2119.4;
            % surf(...
            %     [-1 1; -1 1]*texture_size_in_meters/2,...
            %     [1 1; -1 -1]*texture_size_in_meters/2,...
            %     -0.1*[1 1;1 1],...
            %     imread('ground_texture.jpg'), ... % image must be square
            %     'FaceColor','texturemap',...
            %     'EdgeColor','none' ...
            % );
            
            % Hide axis lines and ticks
            obj.axis_handle.XTickLabel = [];
            obj.axis_handle.YTickLabel = [];
            obj.axis_handle.ZTickLabel = [];
            
            obj.axis_handle.XAxis.Visible = 'off';
            obj.axis_handle.YAxis.Visible = 'off';
            obj.axis_handle.ZAxis.Visible = 'off';
            
            obj.axis_handle.XColor = 'none';
            obj.axis_handle.YColor = 'none';
            obj.axis_handle.ZColor = 'none';
            
            obj.axis_handle.GridColor = [0 0 0];
            obj.axis_handle.GridAlpha = 0.5;
            
            xlabel('')
            ylabel('')
            zlabel('')
            
            xlim([-50 600])
            ylim([-600 600])
            zlim([-600 100])
            
            campos([20 20 -120])
            camtarget([0 0 -100])
            camva(40)
            
            obj.axis_handle.Clipping = 'off';
            grid off
            
            % Scale label
            textarrow = annotation('textarrow',[1 1]*0.18,[.95 .5],'String',...
                'The aircraft are shown 8x larger than their true size.');
            textarrow.Color = 'none';
            textarrow.TextColor = [1 1 1];
            textarrow.FontSize = 14;
            
            % Disable toolbars, since they are also visible in the video
            obj.axis_handle.Toolbar.Visible = 'off';
        end
        
    end
end

