function [ActuatorSS] = fct_Actuator()
% is a function for easy extension to more complex actuators

% easy PT1 gimbal

ActuatorSS = ss(-50, 50, 1, 0);

ActuatorSS.StateName    = {'x_Act'};
ActuatorSS.InputName    = {'mu_cmd'};
ActuatorSS.OutputName   = {'mu'};


