function [gamma_grid, P, P_dot, tDRE] = fct_get_gamma_grid_public(LTVsys_IQC, lambda, gammaLB, gammaUB, P_rho_CL_nom, IQClist, b, T, DRE_event)

T_grid = T;

[SystemMatrices, SysSize]= fct_get_SysMat_LTV_DM_public(P_rho_CL_nom, IQClist(1).IQC, T_grid); % geht gerade nur fuer ein b ... muss angepasst werden

nrIQC = numel(IQClist.IQC); % number of IQCs equivalent to the number of lambdas


%options = odeset('stats', 'on','RelTol',1e-2,'AbsTol',1e-6, 'BDF', 'off'); % adjusted to see the evaluation of the ODE
options = odeset('RelTol',1e-2,'AbsTol',1e-6, 'BDF', 'off', 'Event', DRE_event);

% just works for 2 columns at the moment
% [A,B] = meshgrid(lambda(:,1),lambda(:, 2));
% c=cat(2,A',B');
% x=reshape(c,[],2);

gamma_grid{size(lambda, 1), size(lambda, 1)} = [];
P_dot{size(lambda, 1), size(lambda, 1)} = [];
P{size(lambda, 1), size(lambda, 1)} = [];
tDRE{size(lambda, 1), size(lambda, 1)} =[];


for ii = 1 : 1 : size(lambda, 1)
    
    for jj = 1 : 1 : size(lambda, 2)
        
        [gamma_grid{ii, jj}, P{ii, jj}, P_dot{ii, jj}, tDRE{ii, jj}] = mySearch(lambda{ii, jj}, LTVsys_IQC, T, b, gammaLB, gammaUB, SystemMatrices, IQClist, SysSize, options, DRE_event);
        gamma_grid{ii, jj}
    
    end
    
end


%% Local function
function [gamma_cost, P, P_dot, tDRE] = mySearch(x, LTVsys_IQC, T, b, gammaLBinit, gammaUBinit, SystemMatrices, IQClist, SysSize, options, DRE_event)

lambdaTry = x;
tDRE    = [];
P       = [];

nx = SysSize(1);
% bisection starts here

gammaLB = gammaLBinit;
gammaUB = gammaUBinit;

FirstRun  = 1;
while (gammaUB - gammaLB) > gammaUB * 5*1e-6 % -6 bis jetzt beste resulate /-7 auch gut
    
    if isequal(FirstRun, 1)
    
        gammaTry = gammaUB;
        FirstRun = 0;
        
    else
        
        gammaTry = (gammaLB + gammaUB)/2;
    
    end
    [P_rhs, BXBtrans, CTransC] = fct_build_RDE_components_public(SystemMatrices,IQClist, gammaTry, lambdaTry);
    
    tspan = [T(end) T(1)];
    
    P_0 = zeros(nx, nx);
    %tic
    odefh = @(t, P) fct_build_local_P_dot_IQC(t, P,LTVsys_IQC, SystemMatrices.T_grid,P_rhs, BXBtrans, CTransC, SysSize);
    %toc
    tic % turn on to time evaluation of ode15s
    [tDRE_temp, P_temp] = ode15s(odefh, tspan, P_0, options);
    toc
    
    tDRE_temp = flipud(tDRE_temp);
    
    
    if tDRE_temp(1) > T(1)
        % P did not converge
        gammaLB = gammaTry;
        disp('------------------------FAILURE solution!----------------------') % turn on to see if integration has valid result
        
        
    else
        % P did converge
        disp('------------------------SUCCESS solution!----------------------') % turn on to see if integration has valid result
        gammaUB = gammaTry
        
        % store the new time vector
        tDRE    = tDRE_temp;
        P       = P_temp;
        nt = numel(tDRE);
        
    end
 
   
end

if isempty(tDRE) || isempty(P)
    
    P = zeros(nx, nx);
    tDRE = 0;
    nt = 1;
    P_dot = zeros(nx, nx);
    
else
    
    P = flipud(P);
    P = P';
    P = reshape( P, [nx nx nt]);
    P_dot = zeros(nx, nx);
    
    for i = 1 : 1 : nt
        
        P_dot(:, :, i) = reshape(odefh(tDRE(i), P(:, :, i)), [nx, nx]);
        
    end
    
end
gamma_cost = gammaUB; % scale the cost results in better results of the optimization




%% Local function
function  P_dot = fct_build_local_P_dot_IQC(t, P, LTVsys, T_grid,P_rhs, BXBtrans, CTransC, SysSize)

% Get system matrices

[P_rhs_t, P_lhs_t, CtransC_t, BXBtrans_t] = LTVsys(t, T_grid,P_rhs, BXBtrans, CTransC);

xNr = SysSize(1);

% Convert P from column to matrix

P = reshape(P,[xNr xNr]);

P_dot = CtransC_t + P_lhs_t*P + P*P_rhs_t - P*BXBtrans_t*P; % P_rhs -> A P_lhs -> A^T

P_dot = P_dot(:);

