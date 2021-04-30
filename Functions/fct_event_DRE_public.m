function [Value, isterminal, direction] = fct_event_DRE_public(t, y)

% value(i) mathematical expression describing the ith event. An event
% occurs when value(i) is equal to zero

% isterminal(i) = 1 if the integration is to terminate when the ith event
% occurs, otherwise it is 0

% direction(i) = 0 if all zeros to be detected, +1 if event function is
% increasing (zeros where the event fct is increasing)

isterminal = 1; % stop the integration
direction = 0;
nx = sqrt(length(y));
y = reshape(y, [nx, nx]);
if max(abs(eig(y))) > 10^5
    
   Value = 0;
    
else
    
    Value = 1;
    
end