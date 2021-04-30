%% This script analysis the vanguard space launcher's pitch channel using IQCs
%Author: Nantiwat Pholdee (Khon Kaen University), Felix Biertümpfel (Technische Universität Dresden)
%Last Edited:  30/04/2021
%Email: nantiwat@kku.ac.th, felix.biertuempfel@tu-dresden.de
%Reference: Unpublished work
%DOI:
%% ----------------------------- Init Data --------------------------------
% Init Paths:
clc
clear all

CurrentFolder = cd;
addpath(genpath(CurrentFolder))

LauncherData;

%% ----------------------------- Settings ---------------------------------

% t_grid : - set the time grid for the LTV model
%              -> minimum t = 15s;  Data (11.347s)
%              -> maximum t = 100s; Data (146.35s)
%              -> grid density of the original data is 2.7s, points in between are interpolated

t_grid = [15: 0.5 :100];

% t_ctr: - sets the point in time along the envelope for synthesis of the
%          controller

t_ctr = 50; 

% b: - norm bound on the the LTI uncertainty Delta 

% b = [0.000001, 0.01];
b = [0.01 0.03 0.05 0.06 0.07 0.075 0.08 0.085];


% GammaOverLambda: - runs the algorithm without the optimization for fixed
%                    lambdas
%   1 - this does run a bisection over gamma for a fixed set of lambdas
%   0 - no execution (default as it is time consuming)

GammaOverLambda = 0;


%% ------------------------------ IQC Data --------------------------------

% IQC Filter Psi

IQC.Psi{1} = eye(2); % just a feedthrough block
IQC.Psi{2} = blkdiag(tf(1, [1 1]), tf(1, [1 1])); % PT1

for i = 1 : 1 : length(b)

    IQC.M{i} = blkdiag(b(i)^2, -1); 

end

% IQC list, build structure containing all the IQCs for later use

IQClist = fct_build_IQClist_public(IQC);

% Select the order of the basis function used for the LMI approach

FBasisOrder = 8;


%% ------------ Build the Launchers Pitch Model Components ----------------

% ----- unaugmented launcher
% builds the launchers LTV model (it's not an LPV toolbox model)

[LauncherSS, LauncherSSwind, LauncherSS50s]    = fct_Launcher_Pitch_ss(Launcher, Env, Path, t_grid);

% ----- Gimbal System

[ActuatorSS]    = fct_Actuator();

% ----- Controller

[ControllerSS, ControllerSSdata, ControllerSSdata_DM]  = fct_controller(LauncherSS, LauncherSS50s);

%% ------------------- Build the nominal Closed Loop ----------------------
% ---- nominal closed loop with wind input build like the reference system
% - inputs: mu_ref, wind (alpha disturbance)
% - output: alpha 

ClosedLoopSS_nom_ref = fct_build_closed_CL_nominal_ref(LauncherSSwind, ActuatorSS, ControllerSSdata);

ClosedLoopSS_nom     = fct_build_closed_CL_nominal(LauncherSSwind, ActuatorSS, ControllerSSdata_DM);


%% --------------- Build the system for the IQC Anlysis -------------------

OpenLoopSS_pert_IQC = fct_build_OL_pert_IQC(LauncherSSwind, ActuatorSS, ControllerSSdata_DM);

%% ------------------ Basis Function for LMI Analysis ---------------------

[Fbasis, Fgrad] = fct_build_BasisFct_LTV(OpenLoopSS_pert_IQC, IQClist, FBasisOrder, t_grid);


%% ------------------- Nominal L2 gain integration ------------------------
LTV_sys_nom = @fct_build_LTV_DARE_nom;

[g_nom, t_nom, P_nom, P_dot_nom] = fct_calc_nom_LTV_directInt_public(LTV_sys_nom, t_grid, ClosedLoopSS_nom);
    

%% ------------------- Optional gamma over lambda -------------------------

if isequal(GammaOverLambda, 1)

nrIQC = numel(IQC.Psi); % # number of iqcs defines equivalent to the number of lambdas

LTV_sys_IQC = @fct_build_LTV_DRE_IQC_public;
DRE_event = @fct_event_DRE_public;

% select lambda values
lambda = unique([0.00001, 0.0001, 0.001, 0.005, 0.01 : 0.005 :0.1, 0.01 , 1 : 1 : 10, 10 : 10 : 100, 100 : 50 : 1000, 1000 : 200 : 10000, 10000 : 5000 : 100000, 100000 : 100000 : 1000000, 1000000 : 5000000 : 3e7]');

% build lambda grid
for i = 1 : 1 : length(lambda)
    
    for j = 1 : 1 : length(lambda)
        
        lambda_grid{i, j} = [lambda(i), lambda(j)];
    
    end
    
end

% in the moment just implemented for a single b(k), has to be selected
% accordingly

[gamma_Grid, P_grid, P_dot_grid, tDRE_Grid] = fct_get_gamma_grid(LTV_sys_IQC, lambda_grid, gamma_LB, gamma_UB, OpenLoopSS_pert_IQC, IQClist(1), b(1), t_grid, DRE_event);

end



%% ------------------ Worst Case L2 gain integration ----------------------
nrIQC = numel(IQC.Psi); % # number of iqcs defines equivalent to the number of lambdas

LTV_sys_IQC = @fct_build_LTV_DRE_IQC_public;
DRE_event = @fct_event_DRE_public;

% pre-allocation
g_IQC{length(b)} = [];
x_IQC{length(b)} = [];
P_IQC{length(b)} = [];
P_dot_IQC{length(b)} = [];
tDRE_IQC{length(b)} = [];
fval_IQC{length(b)} = [];
P_IQC{length(b)} = [];
exitflag_IQC{length(b)} = [];
info_IQC{length(b)} = [];

for i = 1 : 1 : length(b) 

    
    if i > 1
        
        lambda_init     = x_IQC{i-1};
        gamma_LB_init   = 0.8*fval_IQC{i-1};
        gamma_UB_init   = 10*fval_IQC{i-1};
    
        
        
    elseif i == 1
        
        
        lambda_init     = 1e5 * ones(1, nrIQC);
        gamma_LB_init   = g_nom;
        gamma_UB_init   = 10*g_nom;
        
        
    end
    
    [g_IQC{i}, x_IQC{i}, P_IQC{i}, P_dot_IQC{i}, tDRE_IQC{i}, fval_IQC{i}, exitflag_IQC{i}, info_IQC{i}] = fct_calc_L2_direct_int_public_MHs(LTV_sys_IQC, lambda_init, gamma_LB_init, gamma_UB_init, OpenLoopSS_pert_IQC, IQClist(i), b(i), t_grid, DRE_event,i);
    
    
end





