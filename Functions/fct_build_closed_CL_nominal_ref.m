function [ClosedLoopSS_nominal] = fct_build_closed_CL_nominal_ref(LauncherSS, ActuatorSS, ControllerSS)

% build the closed loop via sysIC

ControllerSS = ControllerSS;
LauncherSS   = LauncherSS;
ActuatorSS   = ActuatorSS;

systemnames = 'ControllerSS LauncherSS ActuatorSS';

inputvar    = ' [mu_d; wind]';
outputvar   = ' [LauncherSS(1:3)] ';

input_to_ActuatorSS     = ' [ControllerSS(1)] ';
input_to_LauncherSS     = ' [ActuatorSS(1); wind(1)] ';
input_to_ControllerSS   = ' [mu_d; LauncherSS(2)] ';
cleanupsysic = 'yes';

ClosedLoopSS_nominal = sysic;

