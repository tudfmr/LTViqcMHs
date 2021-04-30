% scripts which:
%   - loads the aerodynamic derivatives
%   - calculates and provides the inertias
%   - thrust characteristic
%   - environmental data
%   - prepares data for the look up tables
%   - gets the desired path angle (gamma_d)


%% ------------------ Load the aerodynamic derivatives --------------------
% data is from pg. 291 Fig. 5.18 of the launcher control book
% collected via the grabbit tool

load('D_Launcher.mat')
load('M_alpha_Launcher.mat')
load('Z_alpha_Launcher.mat')
load('M_q_Launcher.mat')

%% ----------------------- Caculate the inertias --------------------------

% ----- create the time vector  -> [11.347, 97.747] step size is 2.7s
%                               -> [100.45, 146.35] step size is 2.7s

t       = [11.347 : 2.7 : 97.747, 100.45 : 2.7 : 146.35];

% ----- calculate the inertias

m_p     = 7727.3 * (1 - 0.0054 * (t - 11.347));
J_yy_p  = 215000 * (1 - 0.006 *(t - 11.347));

% ----- prepare data for look up tables

J_yy    = [t.', J_yy_p.'];
m       = [t.', m_p.'];
eps     = 8.2317; % [m] distance c.g. and nozzle

%% ------------------------ Thrust characteristic -------------------------
% trhust is maintained constant

T = 133202.86; % N

%% -------------------------- Environmental data --------------------------

g_0 = 9.806; %m/s^2

r_0 = 6378.14; % km (earth radius)

%% --------------------------- Trajetory data -----------------------------

v = [0.088604; 0.1065; 0.12437; 0.14287; 0.16166; 0.18074; 0.20021; 0.21983;...
     0.24004; 0.2608; 0.28225; 0.30443; 0.32742; 0.35114; 0.37494; 0.39865;...
     0.42268; 0.44717; 0.47264; 0.49934; 0.52757; 0.55753; 0.58946; 0.62356;...
     0.65985; 0.6983 ; 0.73892; 0.7817; 0.82663; 0.87368; 0.92285; 0.97414;...
     1.0276; 1.0832; 1.1411; 1.2014; 1.2641; 1.3295; 1.3975; 1.4685; 1.5425;...
     1.6196; 1.7002; 1.7843; 1.8722; 1.964; 2.0601; 2.1607; 2.2662; 2.377; 2.4937]*1000; % [m/s]
 

 v_d     = [t', v];
 
 
 gamma = [1.5349; 1.5237; 1.5115; 1.4985; 1.4848; 1.4704; 1.4555; 1.4401; 1.4242;...
          1.408; 1.3915; 1.3747; 1.3577; 1.3407; 1.3235; 1.3063;1.289; 1.2717;...
          1.2544; 1.2372; 1.22; 1.2031; 1.1863; 1.1698; 1.1536; 1.1378; 1.1223;...
          1.1072; 1.0926; 1.0784; 1.0646; 1.0512; 1.0383; 1.0258; 1.0138; 1.0021;...
          0.99093; 0.98013; 0.96973; 0.95973; 0.95013; 0.94091; 0.93206; 0.92359;...
          0.91548; 0.90772; 0.90031; 0.89324; 0.8865; 0.88009; 0.874];
 
 gamma_d = [t', gamma];

      
%% ---------------------------- Data Storage ------------------------------

Launcher.D          = D_Launcher;
Launcher.M_alpha    = M_alpha_Launcher;
Launcher.M_q        = M_q_Launcher;
Launcher.Z_alpha    = Z_alpha_Launcher;
Launcher.J_yy       = J_yy;
Launcher.m          = m;
Launcher.T          = T;
Launcher.eps        = eps;

Path.gamma_d        = gamma_d;
%Path.r_d            = r_d;
Path.v_d            = v_d;
Path.t_d            = t;

Env.g_0             = g_0;
Env.r_0             = r_0;










