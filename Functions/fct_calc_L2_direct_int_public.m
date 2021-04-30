function [cost, x, P, P_dot, tDRE, fval, exitflag, info ] = fct_calc_L2_direct_int_public(LTVsys_IQC, lambda_init, gammaLB, gammaUB, P_rho_CL_nom, IQClist, b, T, DRE_event)

T_grid = T;

[SystemMatrices, SysSize]= fct_get_SysMat_LTV_DM_public(P_rho_CL_nom, IQClist(1).IQC, T_grid);

nrIQC = numel(IQClist.IQC); % number of IQCs equivalent to the number of lambdas

x = lambda_init;

% modify the odeset options
options = odeset('RelTol',1e-2,'AbsTol',1e-6, 'BDF', 'off', 'Event', DRE_event); 

% turn warnings off
warning off

% build the optimization object

myObj = @(x) mySearch(x, LTVsys_IQC,T,gammaLB, gammaUB, SystemMatrices, IQClist, SysSize, options, DRE_event);

% Set the constraints for lambda_i during the optimization
%lb = 0.01 * ones(nrIQC, 1);
lb = diag([0.001, 1e3]) * ones(nrIQC, 1);
ub = inf * ones(nrIQC, 1);

% build the optiprob structure (intermediate structure)

prob = optiprob('fun', myObj, 'bounds', lb, ub, 'x0', x);

% setup the solver options

opts = optiset('solver', 'nlopt', 'maxtime', 0.20*3600, 'tolrfun', 1e-1, 'tolafun', 1e-1, 'solverOpts',nloptset('algorithm', 'LN_COBYLA', 'submaxtime', 3*1800, 'subtolrfun', 1e-2, 'subtolafun', 1e-2));

% build the OPTI object

Opt = opti(prob, opts);

% solve the problem

[x, fval, exitflag, info] = solve(Opt, x);

fval = fval/100000 % rescale the cost of the optimization

% run mysearch again to determine P, P_dot, tDRE and check the cost
[cost, P, P_dot, tDRE] = mySearch(x, LTVsys_IQC, T, gammaLB, fval+100, SystemMatrices, IQClist, SysSize, options, DRE_event); 

cost = cost / 100000 % rescale the cost function

warning on


%% Local function 
function [gamma_cost, P, P_dot, tDRE] = mySearch(x, LTVsys_IQC, T, gammaLBinit, gammaUBinit, SystemMatrices, IQClist, SysSize, options, DRE_event)

lambdaTry = x;
tDRE    = [];
P       = []; 

nx = SysSize(1);

% bisection starts here

gammaLB = gammaLBinit;
gammaUB = gammaUBinit;


FirstRun = 1;
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
    %tic % turn on to time evaluation of ode15s
    [tDRE_temp, P_temp] = ode15s(odefh, tspan, P_0, options);
    %toc
    
    tDRE_temp = flipud(tDRE_temp);
    
    
    if tDRE_temp(1) > T(1)
        
        % P did not converge
        gammaLB = gammaTry;
        %disp('------------------------FAILURE solution!----------------------') % turn on to see if integration has valid result

        
    else
        
        % P did converge
       
        %disp('------------------------SUCCESS solution!----------------------') % turn on to see if integration has valid result
        gammaUB = gammaTry;
        
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

gamma_cost = gammaUB*100000; % scale the cost results in better results of the optimization




%% Local function
function  P_dot = fct_build_local_P_dot_IQC(t, P, LTVsys, T_grid,P_rhs, BXBtrans, CTransC, SysSize)

% Get system matrices

[P_rhs_t, P_lhs_t, CtransC_t, BXBtrans_t] = LTVsys(t, T_grid,P_rhs, BXBtrans, CTransC);

xNr = SysSize(1);

% Convert P from column to matrix

P = reshape(P,[xNr xNr]);

P_dot = CtransC_t + P_lhs_t*P + P*P_rhs_t - P*BXBtrans_t*P; % P_rhs -> A P_lhs -> A^T

P_dot = P_dot(:);

