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

function mystruct_out = map_struct2(mystruct1, mystruct2, operator)
    fields1 = fieldnames(mystruct1);
    fields2 = fieldnames(mystruct2);
    mystruct_out = struct;
    assert(numel(fields1) == numel(fields2));
    for i = 1:length(fields1)
        f1 = fields1{i};
        f2 = fields2{i};
        assert(strcmp(f1,f2))
        mystruct_out.(f1) = operator(mystruct1.(f1), mystruct2.(f2));
    end
end

