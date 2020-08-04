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

function c = collocation_constants(n_nodes)

    N = n_nodes - 1;

    % Legendre-Gauss-Lobatto nodes
    P = zeros(N+1, N+1);
    P(1,1) = 1;
    P(2,2) = 1;
    for k=2:N
        P(:,k+1)=( (2*k-1) * [0;P(1:end-1,k)]  - (k-1) * P(:,k-1) ) /k;
    end
    dP = P(2:end,end)' .* (1:N);
    c = -dP(1:end-1)/dP(end);
    A = diag(ones(length(c)-1,1),1);
    A(end,:) = c;
    nodes = sort([eig(A); -1; 1]);
    
    % Legendre-Gauss-Lobatto weights
    P = zeros(N+1, N+2);
    P(:,1)=1;
    P(:,2)=nodes;
    for k=2:(N+1)
        P(:,k+1)=(( (2*k-1)*nodes.*P(:,k)-(k-1)*P(:,k-1) )/k);
    end
    weights = 2./(N*(N+1)*P(:,(N+1)).^2);

    
    
    integration_matrix = zeros(N+1,N+1);
    for j = 2:(N+1)
        integration_matrix = integration_matrix + P(:,j)' .* (P(:,j+1) - P(:,j-1));
    end
    integration_matrix = (integration_matrix + nodes + 1) .* weights' / 2;
    integration_matrix = integration_matrix(2:end,:);
    
    % Switch to [0, 1] interval
    nodes = (nodes + 1) / 2;
    weights = weights / 2;
    integration_matrix = integration_matrix / 2;
    
    % Differentiation matrix
    deltas = nodes - nodes' + eye(length(nodes));
    P = prod(deltas, 2);
    interpolation_weights = 1 ./ P;
    Q = P ./ deltas;
    D = interpolation_weights' .* Q + diag(-2 + interpolation_weights .* sum(Q, 2));

    % Interpolation matrix
    t_interp = linspace(0, 1, 20*n_nodes + 1);
    deltas = (t_interp - nodes)';
    interpolation_matrix = 0*deltas;
    for j = 1:length(nodes)
        deltas_copy = deltas;
        deltas_copy(:,j) = [];
        interpolation_matrix(:,j) = interpolation_weights(j) * prod(deltas_copy,2);
    end

    c = struct;
    c.D                    = D;
    c.nodes                = nodes;
    c.n_nodes              = n_nodes;
    c.integration_weights  = weights;
    c.integration_matrix   = integration_matrix;
    c.interpolation_nodes  = t_interp;
    c.interpolation_matrix = interpolation_matrix;
    
end

