function ClosedLoopSS_pert_IQC = fct_build_OL_pert_IQC(LauncherSS, ActuatorSS, ControllerSS)


systemnames = 'ControllerSS LauncherSS ActuatorSS';

inputvar    = ' [ w{1} ;d{1} ]'; % wind, it's the disturbance input d
outputvar   = ' [ ControllerSS(1) + w; LauncherSS(1) ] '; % alpha, it's the disturbance output e 

input_to_ActuatorSS     = ' [w + ControllerSS(1) + w] ';
input_to_LauncherSS     = ' [ActuatorSS(1); d] ';
input_to_ControllerSS   = ' [LauncherSS(2)] ';
cleanupsysic = 'yes';

ClosedLoopSS_pert_IQC = sysic;