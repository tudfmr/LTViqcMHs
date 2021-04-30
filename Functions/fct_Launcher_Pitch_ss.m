function [LauncherSS, LauncherSSwind, LauncherSS50s] = fct_Launcher_Pitch_ss(Launcher, Env, Path, t_grid)

% builds the Launchers_Pitch_ss
% - it creates 


% the simplified pitch state space model looks like:
% - M_alpha_dot ? 0

%                   |  Z_alpha       g*sin(theta_d)         |       
%   alpha_dot       | ---------   - ----------------    1   | alpha
%                   |  m * v_d            v_d               |
%                   |                                       |
%   theta_dot   =   |     0                0            1   | theta
%                   |                                       |
%                   |   M_alpha                        M_q  |    
%   q_dot           |  ---------           0          ----- |   q
%                   |    J_yy                          J_yy |
%
%
%                               |   T   |
%                               | ----- |
%                               | m*v_d |
%                               |       |
%                               |   0   | mu
%                               |       |
%                               | T*eps |
%                               | ----- |
%                               |  J_yy |


%% ------------- generate data specific to the desired grid ---------------

Z_alpha = interp1(Launcher.Z_alpha(:, 1), Launcher.Z_alpha(:, 2), t_grid);

M_alpha = interp1(Launcher.M_alpha(:, 1), Launcher.M_alpha(:, 2), t_grid); % [Nm/ rad]

M_q     = interp1(Launcher.M_q(:, 1), Launcher.M_q(:, 2), t_grid); % [s]

J_yy    = interp1(Launcher.J_yy(:, 1), Launcher.J_yy(:, 2), t_grid); % [kg m^2]

m       = interp1(Launcher.m(:, 1), Launcher.m(:, 2), t_grid);

gamma_d = interp1(Path.gamma_d(:, 1), Path.gamma_d(:, 2), t_grid);

v_d     = interp1(Path.v_d(:, 1), Path.v_d(:, 2), t_grid);


%% -------------------- build the state space models ----------------------
% the implementation is based on the simulink model given on pg. 295 Fig.
% 5.22

% check the units
% there is an error of factor 1000 in the alpha channel


for i = 1 : 1 : length(t_grid)

% ----- A
    
    A(:, :, i) = [Z_alpha(i)/(m(i) * v_d(i)), -(Env.g_0/ v_d(i)) .* sin(gamma_d(i)),      1        ;...
                                   0        ,                0                 ,      1        ;...
                      M_alpha(i)/J_yy(i)    ,                0                 , M_q(i) / J_yy(i)  ];
    
% ----- B

    B(:, :, i) = [     Launcher.T/(m(i)*v_d(i))      ;...
                                 0                   ;...
                  (Launcher.T * Launcher.eps)/J_yy(i) ];


% ----- C
    
    C(:, :, i) = eye(3, 3);
    
% ----- D
    
    D(:, :, i) = zeros(3, 1);


end


LauncherSS = ss(A, B, C, D);


LauncherSS.InputName    = 'mu';
LauncherSS.OutputName   = {'alpha', 'theta', 'q'};
LauncherSS.StateName    = {'alpha', 'theta', 'q'};



%% Launcher model with wind input
% wind input is a alpha disturbance, therefore the wind input is modeled
% the same way as the alpha channel

for i = 1 : 1 : length(t_grid)

% ----- A
    
    A_wind(:, :, i) = [Z_alpha(i)/(m(i) * v_d(i)), -(Env.g_0/ v_d(i)) .* sin(gamma_d(i)),      1        ;...
                                   0        ,                0                 ,      1        ;...
                      M_alpha(i)/J_yy(i)    ,                0                 , M_q(i) / J_yy(i)  ];
    
% ----- B

    B_wind(:, :, i) = [     Launcher.T/(m(i)*v_d(i)), Z_alpha(i)/(m(i) * v_d(i));...
                                       0            ,            0              ;...
                       (Launcher.T * Launcher.eps)/J_yy(i) , M_alpha(i)/J_yy(i)];


% ----- C
    
     C_wind(:, :, i) = eye(3, 3);
    %C_wind(:, :, i) = [0, 0, 0; 0, 1, 0; 0, 0, 0];
    
% ----- D
    
    D_wind(:, :, i) = zeros(3, 2);


end

LauncherSSwind = ss(A_wind, B_wind, C_wind, D_wind);

LauncherSSwind.InputName    = {'mu', 'wind'};
LauncherSSwind.OutputName   = {'alpha', 'theta', 'q'};
LauncherSSwind.StateName    = {'alpha', 'theta', 'q'};


%% ------------------------ generate data @50s ----------------------------

Z_alpha_50s = interp1(Launcher.Z_alpha(:, 1), Launcher.Z_alpha(:, 2), 50);

M_alpha_50s = interp1(Launcher.M_alpha(:, 1), Launcher.M_alpha(:, 2), 50); % [Nm/ rad]

M_q_50s     = interp1(Launcher.M_q(:, 1), Launcher.M_q(:, 2), 50); % [s]

J_yy_50s    = interp1(Launcher.J_yy(:, 1), Launcher.J_yy(:, 2), 50); % [kg m^2]

m_50s       = interp1(Launcher.m(:, 1), Launcher.m(:, 2), 50);

gamma_d_50s = interp1(Path.gamma_d(:, 1), Path.gamma_d(:, 2), 50);

v_d_50s     = interp1(Path.v_d(:, 1), Path.v_d(:, 2), 50);

%% ----------------------- Launcher Model @50s ----------------------------

% check the units
% there is an error of factor 1000 in the alpha channel


for i = 1 : 1 : length(t_grid)

% ----- A
    
    A_50s = [Z_alpha_50s/(m_50s * v_d_50s), -(Env.g_0/ v_d_50s) .* sin(gamma_d_50s),      1        ;...
                                   0        ,                0                 ,      1        ;...
                      M_alpha_50s/J_yy_50s    ,                0                 , M_q_50s / J_yy_50s  ];
    
% ----- B

    B_50s = [     Launcher.T/(m_50s*v_d_50s)      ;...
                                 0                   ;...
                  (Launcher.T * Launcher.eps)/J_yy_50s ];


% ----- C
    
    C_50s = eye(3, 3);
    
% ----- D
    
    D_50s = zeros(3, 1);


end


LauncherSS50s = ss(A_50s, B_50s, C_50s, D_50s);


LauncherSS50s.InputName    = 'mu';
LauncherSS50s.OutputName   = {'alpha', 'theta', 'q'};
LauncherSS50s.StateName    = {'alpha', 'theta', 'q'};





