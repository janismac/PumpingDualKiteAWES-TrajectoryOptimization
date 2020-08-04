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

classdef MyVideoWriter < handle
    
    
    properties
        cleanupObj
        videoWriter
        width
        height
        filename
        frame_count
        high_quality
    end
    
    methods
        function obj = MyVideoWriter(high_quality)
            
            
            if nargin < 1
                high_quality = false;
            end
            
            obj.high_quality = high_quality;

            obj.filename = ['output/video_' datestr(now,'yyyy_mm_dd_HH_MM_SS')];
            
            if high_quality
                obj.filename = [obj.filename '_HD'];
            end
            
            obj.videoWriter = VideoWriter(obj.filename,'MPEG-4');
            obj.videoWriter.Quality = 99;
            obj.videoWriter.open();
            obj.frame_count = 0;

            obj.width = 1280;
            obj.height = 720;
        end
        
        function img = getFrame(obj)
            drawnow
            
            set(gcf,'units','pixel');
            old_sz = get(gcf,'position');
            if old_sz(3) ~= obj.width || old_sz(4) ~= obj.height
                new_sz = old_sz .* [1 1 0 0] + [0 0 obj.width obj.height];
                set(gcf,'position',new_sz);
                drawnow
            end
            
            sz = get(gcf,'position');
            assert(sz(3) == obj.width);
            assert(sz(4) == obj.height);
            
            if obj.high_quality
                frame_file = [getenv('TEMP') '/frame' num2str(feature('getpid')) '.png'];
                print(frame_file,'-opengl','-dpng','-r288');
                img = imread(frame_file);
                assert(size(img,1) == 1080*2);
                assert(size(img,2) == 1920*2);
                img = imresize(img,0.5);
                delete(frame_file);
            else
                frame = getframe(gcf);
                img = frame.cdata;
            end
        end
        
        
        
        
        function writeFrame(obj, img)
            
            if nargin < 2
                img = obj.getFrame();
            end
            
            writeVideo(obj.videoWriter,img);
            obj.frame_count = obj.frame_count + 1;
        end
    end
end

