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

function setup_3D_plot(fig_handle)

    if nargin == 0
        fig_handle = gcf;
    end
    set(0, 'CurrentFigure', fig_handle)
    clf
    hold on
    grid on
    box on
    camproj perspective
    daspect([1 1 1])
    cameratoolbar('Show')
    cameratoolbar('SetMode', 'orbit')
    
    set(gca, 'ZDir','reverse')
    set(gca, 'YDir','reverse')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    set(fig_handle,'units','normalized')
    sz = get(fig_handle,'outerposition');
    if sz(3) < 0.5
        set(fig_handle,'outerposition',[0.1 0.1 0.8 0.8])
        %fig_handle.WindowState = 'maximized';
    end
    
end

