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

function S = make_subvector_slice_struct(S, prefix, sym_vec)
    assert(length(sym_vec) == numel(sym_vec));
    
    for i = 1:length(sym_vec)
        sym_names{i} = char(sym_vec(i));
    end
    
    similarity_matrix = false(length(sym_vec), length(sym_vec));
    
    pos_digits = {'1','2','3','4','5','6','7','8','9'};
    for i = 1:length(sym_vec)
        for j = (i+1):length(sym_vec)
            if length(sym_names{i}) > 2 && length(sym_names{i}) == length(sym_names{j})
                diff_slice = sym_names{i} ~= sym_names{j};
                hamming_distance = sum(diff_slice);
                assert(hamming_distance > 0);
                if hamming_distance == 1 ...
                && any(strcmp(pos_digits, sym_names{i}(diff_slice))) ...
                && any(strcmp(pos_digits, sym_names{j}(diff_slice))) 
                    similarity_matrix(i,j) = true;
                    similarity_matrix(j,i) = true;
                end
            end
        end
        similarity_matrix(i,i) = true;
    end
    
    for i = 1:length(sym_vec)
        I = find(similarity_matrix(:,i));
        if numel(I) > 1
            name_slice = sym_names{I(1)} == sym_names{I(2)};
            name = [prefix sym_names{I(1)}(name_slice)];
            S.(name) = I';
        elseif numel(I) == 1 && any(sym_names{I(1)} == '1')
            name = [prefix sym_names{I(1)}(~(sym_names{I(1)} == '1'))];
            S.(name) = I';
        end
    end
    
end