function [ControllerSS, ControllerSSdata, ControllerSSdata_DM]  = fct_controller(LauncherSS, LauncherSS50s)

% Controller

% the controller is taken from the exaple given on pg. 292ff
% it's designed for the max pressure point @ t = 50s

A = [-0.0356, -0.0233, 1.000  ;...
    0  ,     0  ,    1   ;...
    1.1961,     0  , -0.0341];

B = [0.0522; 0; 7.0084];

C = [0, 1, 0];

D = 0;

% ----- regulator design by LQR

[k, S, E] = lqr(A, B, [0.01, 0, 0; 0, 0.1, 0;  0, 0, 0.1], 1); %Regulator design by LQR

% Observer (Kalman Filter) design by LQR

Lp   = (lqr(A', C', eye(3), 1))'; % L = Lp'


% Controller SS
% Observer: x_hat_dot = (A - LC)x_hat + Be + Ly;
% Controller is just a pure prop of the Observer Output;

A_ctr = A - Lp * C - B*k;
B_ctr = [B*k(:, 2), Lp];
C_ctr = -k;
D_ctr = [k(:, 2), 0];

ControllerSS            = ss(A_ctr, B_ctr, C_ctr, D_ctr);
ControllerSS.StateName  = {'x_ctr1', 'x_ctr2', 'x_ctr3'};


ControllerSS.OutputName = {'mu_cmd'};
ControllerSS.InputName  = {'mu_d', 'q'};

%% Controller via data from data set @50s


A_data = LauncherSS50s.A;
B_data = LauncherSS50s.B;
C_data = [0, 1, 0];
D_data = 0;

% ----- regulator design by LQR

[k_data, S_data, E_data] = lqr(A_data, B_data, [0.01, 0, 0; 0, 0.1, 0;  0, 0, 0.1], 1); %Regulator design by LQR

% Observer (Kalman Filter) design by LQR

Lp_data   = (lqr(A_data', C_data', eye(3), 1))'; % L = Lp'


%% Controller SS according to reference system
% Observer: x_hat_dot = (A - LC)x_hat + Be + Ly;
% Controller is just a pure prop of the Observer Output;

A_ctr_data = A_data - Lp_data * C_data - B_data*k_data;
B_ctr_data = [B_data*k_data(:, 2), Lp_data];
C_ctr_data = -k_data;
D_ctr_data = [k_data(:, 2), 0];

ControllerSSdata            = ss(A_ctr_data, B_ctr_data, C_ctr_data, D_ctr_data);
ControllerSSdata.StateName  = {'x_ctr1', 'x_ctr2', 'x_ctr3'};


ControllerSSdata.OutputName = {'mu_cmd'};
ControllerSSdata.InputName  = {'mu_d', 'theta'};

%% Controller SS without reference signal, used for the robustness analysis

A_ctr_data_DM = A_data - Lp_data * C_data - B_data*k_data;
B_ctr_data_DM = [Lp_data];
C_ctr_data_DM = -k_data;
D_ctr_data_DM = [0];

ControllerSSdata_DM            = ss(A_ctr_data_DM, B_ctr_data_DM, C_ctr_data_DM, D_ctr_data_DM);
ControllerSSdata_DM.StateName  = {'x_ctr1', 'x_ctr2', 'x_ctr3'};


ControllerSSdata_DM.OutputName = {'mu_cmd'};
ControllerSSdata_DM.InputName  = {'theta'};









