function ClosedLoopSS_pert_LTI_WC = fct_build_CL_pert_LTI_WC(LauncherSS, ActuatorSS, ControllerSS)

Delta = ultidyn('Delta', [1, 1], 'bound', 1);

systemnames = 'ControllerSS LauncherSS ActuatorSS Delta';

inputvar    = ' [wind]';
outputvar   = ' [LauncherSS(1)] ';

input_to_ActuatorSS     = ' [Delta(1)+ ControllerSS(1) + Delta(1)] ';
input_to_Delta          = ' [ControllerSS(1) + Delta(1)] ';
input_to_LauncherSS     = ' [ActuatorSS(1); wind(1)] ';
input_to_ControllerSS   = ' [LauncherSS(2)] ';
cleanupsysic = 'yes';

ClosedLoopSS_pert_LTI_WC = sysic;