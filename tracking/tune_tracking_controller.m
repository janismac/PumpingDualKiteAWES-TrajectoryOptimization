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

function K = tune_tracking_controller(A,B,dt,Q,R)
    N = size(A,3);

    reps = 5;

    A = repmat(A,1,1,reps);
    B = repmat(B,1,1,reps);

    K = lqr_continuous_finite_horizon(A,B,dt,Q,R);
    K = K(:,:,1:N);
end

