function [cost, x, P, P_dot, tDRE, fval, exitflag, info ] = fct_calc_L2_direct_int_public_MHs(LTVsys_IQC, lambda_init, gammaLB, gammaUB, P_rho_CL_nom, IQClist, b, T, DRE_event,bno)

T_grid = T;

[SystemMatrices, SysSize]= fct_get_SysMat_LTV_DM_public(P_rho_CL_nom, IQClist(1).IQC, T_grid);

nrIQC = numel(IQClist.IQC); % number of IQCs equivalent to the number of lambdas

x = lambda_init;

% modify the odeset options
options = odeset('RelTol',1e-2,'AbsTol',1e-6, 'BDF', 'off', 'Event', DRE_event); 
% 
% % turn warnings off
% warning off
% 

%%
%%%%%%%%%%%%%% original implementation of the benchmark %%%%%%%%%%%%%%%%%%%
% Remarks:
%   - related to the authors' paper https://ieeexplore.ieee.org/document/8511591
%   - it provides the lower and upper bounds for gamma in the bisection
%   - to run this section the latest version of the OPTI-toolbox needs to be installed into the main directory:
%       -> https://www.inverseproblem.co.nz/OPTI/index.php/DL/DownloadOPTI
%       -> https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
%          Copyright (C) 2011-2018, Jonathan Currie, All rights reserved.





% Set the constraints for lambda_i during the optimization
lb = diag([0.001, 1e3]) * ones(nrIQC, 1);
ub = inf * ones(nrIQC, 1);
myObj = @(x) mySearch(x, LTVsys_IQC,T,gammaLB, gammaUB, SystemMatrices, IQClist, SysSize, options, DRE_event);

% build the optiprob structure (intermediate structure)


try

prob = optiprob('fun', myObj, 'bounds', lb, ub, 'x0', x);

catch
    
    warning('Calculation not possible. The OPTI-toolbox is not installed please visit: https://www.inverseproblem.co.nz/OPTI/index.php/DL/DownloadOPTI for installation')
    
    cost     = [];
    x        = [];
    P        = [];
    P_dot    = [];
    tDRE     = [];
    fval     = [];
    exitflag = [];
    info     = [];
    
    return
    
end

% setup the solver options

opts = optiset('solver', 'nlopt', 'maxtime', 0.20*3600, 'tolrfun', 1e-1, 'tolafun', 1e-1, 'solverOpts',nloptset('algorithm', 'LN_COBYLA', 'submaxtime', 3*1800, 'subtolrfun', 1e-2, 'subtolafun', 1e-2));

% build the OPTI object

Opt = opti(prob, opts);

%solve the problem
tic
[x, fval, exitflag, info] = solve(Opt, x);
time=toc;
save(['OrignalbNo' num2str(bno)],'x','fval','exitflag','info','time')

fval = fval; % rescale the cost of the optimization

%run mysearch again to determine P, P_dot, tDRE and check the cost
[cost,fobj,gcon, P, P_dot, tDRE] = mySearch(x, LTVsys_IQC, T, gammaLB, fval+100, SystemMatrices, IQClist, SysSize, options, DRE_event); 

cost = cost ; % rescale the cost function



%%
%%%%%% MHs
gammaLBmh=0.1*gammaLB;
gammaUBmh=0.5*gammaUB;
myObjMH = @(x) mySearchMH(x, LTVsys_IQC,T,gammaLBmh, gammaUBmh, SystemMatrices, IQClist, SysSize, options, DRE_event);
myObjMH02 = @(x,gammaUBmh) mySearchMH(x, LTVsys_IQC,T,gammaLBmh, gammaUBmh, SystemMatrices, IQClist, SysSize, options, DRE_event);
myObjMH03 = @(x,fold) mySearchMH02(x, LTVsys_IQC,T,gammaLBmh, gammaUBmh, SystemMatrices, IQClist, SysSize, options, DRE_event,fold);



nbit=5;
nvar=nrIQC;

algo=[{'SODE'}
    {'SOJADE'}
    {'SOSHADE'}
    {'SOLSHADE'}
    {'SOSCA'}    %%%5
    {'SODA'}        
    {'SOGWO'}
    {'SOMFO'}
    {'SOWOA'}
    {'SOISCA'}     %%%%10
    {'SOm_SCA'}
    {'SOSinDE'}
    {'SOLSHADEND'}
    {'SOSPS_LSHADE_EIG'}
    {'SOSCAadbL02'}   %%%% 15
    {'SODEobj3'}   
    {'SOJADEobj3'}
    {'SOSHADEobj3'}
    {'SOLSHADEobj3'}   
    {'SOSCAobj3'}   %%%20
    {'SODAobj3'}        
    {'SOGWOobj3'}
    {'SOMFOobj3'}
    {'SOWOAobj3'}
    {'SOSCAadbL02obj3'}  %%%% 25
    {'SOISCAobj3'}    
    {'SOm_SCAobj3'}
    {'SOLSHADENDobj03'}  
    {'SOSPS_LSHADE_EIGobj03'}
    {'SOSinDEobj03'}    %%%%30
];

for i=1:size(algo,1)
    if i==15||i==25
        lb =0*ones(nrIQC, 1);
        ub =1e2*ones(nrIQC, 1);
    else
        lb =0*ones(nrIQC, 1);
        ub =1e8*ones(nrIQC, 1);
    end
    parfor j=1:5
        [b,i,j]
        if i==15
            nloop=25;  %%% defult is 50x50
            nsol=100;
            feval(char(algo(i,:)),myObjMH02,['bNo' num2str(bno) 'Al' num2str(i) 'Run' num2str(j)],nloop,nsol,nvar,nbit,lb,ub,gammaUBmh)
        elseif i==25
            nloop=25;  %%% defult is 50x50
            nsol=100;
            feval(char(algo(i,:)),myObjMH03,['bNo' num2str(bno) 'Al' num2str(i) 'Run' num2str(j)],nloop,nsol,nvar,nbit,lb,ub,gammaUBmh)
        elseif ismember(i,[16:24,26:30])
            nloop=100;  %%% defult is 50x50
            nsol=25;
            feval(char(algo(i,:)),myObjMH03,['bNo' num2str(bno) 'Al' num2str(i) 'Run' num2str(j)],nloop,nsol,nvar,nbit,lb,ub)
        else
            nloop=50;  %%% defult is 50x50
            nsol=50;
            feval(char(algo(i,:)),myObjMH,['bNo' num2str(bno) 'Al' num2str(i) 'Run' num2str(j)],nloop,nsol,nvar,nbit,lb,ub)
        end
    end
end

warning on


%% Local function 
function [gamma_cost,fobj,gcon, P, P_dot, tDRE] = mySearchMH(x, LTVsys_IQC, T, gammaLBinit, gammaUBinit, SystemMatrices, IQClist, SysSize, options, DRE_event)

% a = 0 * ones(length(x), 1);
% b = 1e6 * ones(length(x), 1);

% % a = diag([-3,3]) * ones(length(x), 1);
% % b = 8 * ones(length(x), 1);
% x=a+(b-a).*x;
% x=10.^x;

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
    
    gamma_cost = gammaUB*1e10;
    fobj=gammaUB;
    gcon=0;
    
else
    P = flipud(P);
    P = P';
    P = reshape( P, [nx nx nt]);
    P_dot = zeros(nx, nx);
    
    for i = 1 : 1 : nt
        
        P_dot(:, :, i) = reshape(odefh(tDRE(i), P(:, :, i)), [nx, nx]);
        
    end
    
    gamma_cost = gammaUB;
    fobj=gamma_cost;
    gcon=0;
    
end


%% Local function 
function [gamma_cost,fobj,gcon, P, P_dot, tDRE] = mySearchMH02(x, LTVsys_IQC, T, gammaLBinit, gammaUBinit, SystemMatrices, IQClist, SysSize, options, DRE_event,fold)

if fold<gammaUBinit
    gammaUBinit=fold;
end
% a = 0 * ones(length(x), 1);
% b = 1e6 * ones(length(x), 1);

% % a = diag([-3,3]) * ones(length(x), 1);
% % b = 8 * ones(length(x), 1);
% x=a+(b-a).*x;
% x=10.^x;

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
    
    gamma_cost = gammaUB*1e10;
    fobj=gammaUB;
    gcon=0;
    
else
    P = flipud(P);
    P = P';
    P = reshape( P, [nx nx nt]);
    P_dot = zeros(nx, nx);
    
    for i = 1 : 1 : nt
        
        P_dot(:, :, i) = reshape(odefh(tDRE(i), P(:, :, i)), [nx, nx]);
        
    end
    
    gamma_cost = gammaUB;
    fobj=gamma_cost;
    gcon=0;
    
end


%% Local function 
function [gamma_cost,fobj,gcon, P, P_dot, tDRE] = mySearch(x, LTVsys_IQC, T, gammaLBinit, gammaUBinit, SystemMatrices, IQClist, SysSize, options, DRE_event)

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

gamma_cost = gammaUB;
fobj=gamma_cost;
gcon=0;

%% Local function
function  P_dot = fct_build_local_P_dot_IQC(t, P, LTVsys, T_grid,P_rhs, BXBtrans, CTransC, SysSize)

% Get system matrices

[P_rhs_t, P_lhs_t, CtransC_t, BXBtrans_t] = LTVsys(t, T_grid,P_rhs, BXBtrans, CTransC);

xNr = SysSize(1);

% Convert P from column to matrix

P = reshape(P,[xNr xNr]);

P_dot = CtransC_t + P_lhs_t*P + P*P_rhs_t - P*BXBtrans_t*P; % P_rhs -> A P_lhs -> A^T

P_dot = P_dot(:);

