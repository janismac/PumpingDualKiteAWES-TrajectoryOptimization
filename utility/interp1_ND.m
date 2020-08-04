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

function xq = interp1_ND(t,x,tq)
    
    % 1-D interpolation of an ND array.
    % The interpolation is carried out along the unique dimension
    % for which size(x, i) == length(t).
    % The interpolation is parallelized/vectorized with respect to the
    % other dimensions.

    assert(numel(t) == length(t));
    assert(numel(tq) == length(tq));
    N = numel(t);
    
    sz = size(x);
    I = find(sz==N);
    assert(numel(I) == 1);
    
    dim_swap = 1:length(sz);
    dim_swap([I 1]) = dim_swap([1 I]);
    x_swapped = permute(x, dim_swap);
    sz_swapped = size(x_swapped);
    x_flat = reshape(x_swapped, N, numel(x)/N);
    xq_flat = interp1(t(:), x_flat, tq(:));
    sz_swapped(1) = numel(tq);
    xq_swapped = reshape(xq_flat, sz_swapped);
    xq = permute(xq_swapped, dim_swap);
    
end

