function [g,t_feas, P_feas, P_dot_feas] = fct_calc_nom_LTV_directInt_public(LTVsys, T_grid, P_rho_CL)

g = inf;

P = [];

t = []; % discrete time steps of the ODE45 solver

[A, B, C, D, M] = fct_get_SysMat_LTV(P_rho_CL);


%% Bisection Bounds and Tolerance

gLower = 0;
gUpper = 10;

tol = 1.0e-4;

%% Upper Bound Phase

cnt             = 0;
cntMax          = 10;
UpperBndFound   = false;

warning off

while ~UpperBndFound && cnt <= cntMax
    
    gTry = gUpper; % initial guess
    
    % solve LTV Ricatti Equation
    
    tic
    [t, P, P_dot] = fct_solve_DARE_LTV_2(LTVsys, gTry, T_grid, A, B, C, D, M);
    toc
    
    if t(1) > T_grid(1) % if no solution the time vector is not flipped to 0 -> T
        % P did not converge
        gUpp = 2 * gUpper; % try bigger initial guess (way higher norm bound)
        
    else
        % P did converge
        UpperBndFound = true;
        
    end
    
    cnt = cnt + 1;
    
end

if ~UpperBndFound % no solution found LTV system is not stable
    
    g = inf;
    warning on;
    return
    
end

%% Bisection Phase

while (gUpper - gLower > tol * gUpper)
    
    gTry = (gUpper + gLower)/2;
    
    % solve the LTV DARE
    
    [t, P, P_dot] = fct_solve_DARE_LTV_2(LTVsys, gTry, T_grid, A, B, C, D, M);
    
    if t(1) > T_grid(1)
        % P did not converge
        gLower = gTry;
        disp('------------------------FAILURE solution!----------------------')
    else
        % P converged
        gUpper = gTry;
        
        % to store last feasible solution
        t_feas      = t;
        P_feas      = P;
        P_dot_feas  = P_dot;
        disp('------------------------SUCCESS solution!----------------------')
        
    end
    
    
end

% store final results

g = gUpper;

warning on


%% Local function (may solve that via a function handle)
% solves the ricatti equetion for T = t_i
function [t, P, P_dot] = fct_solve_DARE_LTV_2(LTVsys, gTry, T_grid ,A, B, C, D, M)

xNr = size(A, 1);

uNr = size(B, 2); % nr of cloumns equals number of inputs

ODEfh = @(t, P) fct_build_local_P_dot(t, P, LTVsys, gTry, T_grid ,A, B, C, D, M);
options = odeset('stats', 'on', 'RelTol',1e-1,'AbsTol',1e-3);


P0      = zeros(xNr, xNr);


tspan   = [T_grid(end) T_grid(1)];

[t, P]  = ode15s(ODEfh, tspan, P0); % nr 1


% Reverse time back from tau to t and then reshape bacj to matrix array (to start from 0 again)

t   = flipud(t); % flip the vector upside down
tNr  = numel(t);

P   = flipud(P);
%P   = reshape(shiftdim(P, -1), [xNr, xNr, tNr]); % shift dimensions to the right and pads with singletons

P = P';

P = reshape(P, [xNr, xNr, tNr]);

% return P_dot

P_dot = zeros(xNr, xNr, tNr);
I = eye(uNr);

for i = 1 : 1 : tNr
    
    P_dot(:, :, i) = reshape(ODEfh(t(i), P(:, :, i)), [xNr, xNr]);
    
end


%% Local function
function  P_dot = fct_build_local_P_dot(t, P, LTVsys, gTry, T_grid, A, B, C, D, M)

xNr = size(A,1);
uNr = size(B,2);
yNr = size(C, 1);
mdlNr = size(A,3);
% Get system matrices
[A_int,B_int,C_int,D_int] = LTVsys(t, T_grid, A, B, C, D, M,xNr, uNr, yNr, mdlNr);
I = eye(uNr);

% Convert P from column to matrix
P = reshape(P,[xNr xNr]);

% Compute Pdot
%P_dot = -P*A_int + -A_int'*P + (P*B_int+C_int'*D_int)*((D_int'*D_int-gTry^2*I)\(P*B_int+C_int'*D_int)') - C_int'*C_int;
P_dot = -P*A_int + -A_int'*P + (P*B_int+C_int'*D_int)*(inv(D_int'*D_int-gTry^2*I)*(P*B_int+C_int'*D_int)') - C_int'*C_int;

% Convert Pdot from matrix back to column
P_dot = P_dot(:);



