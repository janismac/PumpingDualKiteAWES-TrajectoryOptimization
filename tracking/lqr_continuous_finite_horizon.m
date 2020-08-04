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

function K = lqr_continuous_finite_horizon(A,B,dt,Q,R)
% Source "Linear, time-varying approximations to nonlinear dynamical
% systems: with applications in control and optimization" by Maria
% Tomas-Rodriguez, Stephen P. Banks, page 105

    Pf = Q;
    yf = Pf(:);
    n = size(A,1);
    N = size(A,3);
    
    Rinv = inv(R);
    
    my_ode = @(t,y) riccati_ode(t,dt,y,A,B,Q,Rinv,n);
    
    t_grid = fliplr((0:(N-1))*dt);
    
    [~,YOUT] = ode45(my_ode, t_grid, yf);
    YOUT = YOUT';
    YOUT = fliplr(YOUT);
    t_grid = fliplr(t_grid);
    
    P = reshape(YOUT, n, n, length(t_grid));
    K = nan(size(B,2), size(B,1), length(t_grid));
    
    for i = 1:length(t_grid)
        K(:,:,i) = Rinv * B(:,:,i)' * P(:,:,i);
    end
    
end



function dydt = riccati_ode(t,dt,y,A,B,Q,Rinv,n)

    % Linear interpolation of A,B
    i1 = floor(t/dt+1);
    i1 = max(1, i1);
    i1 = min(size(A, 3), i1);
    
    i2 = ceil(t/dt+1);
    i2 = max(1, i2);
    i2 = min(size(A, 3), i2);
    
    if i1 < i2
        tau2 = t/dt - (i1-1);
        tau1 = 1-tau2;
        
        A_now = tau1 * A(:,:,i1) + tau2 * A(:,:,i2);
        B_now = tau1 * B(:,:,i1) + tau2 * B(:,:,i2);
    else
        A_now = A(:,:,i1);
        B_now = B(:,:,i1);
    end
    
    % Riccati ODE
    P = reshape(y, n, n);
    dPdt = P*A_now + A_now'*P + Q - P*B_now*Rinv*B_now'*P;
    dPdt = -dPdt;
    dydt = dPdt(:);
end



