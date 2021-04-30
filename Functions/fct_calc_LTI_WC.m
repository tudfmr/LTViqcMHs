function wcglti = fct_calc_LTI_WC(LauncherSS_grid, ActuatorSS, ControllerSS, b)

% number of models

mdls = length(LauncherSS_grid.A);

%at each grid point compute worst case LTI gain using wcgain.m
for jj = 1:length(b)
    kappa=b(jj);
    clear delta
    Delta = ultidyn('Delta',[1,1],'Bound',kappa);
    
    for ii = 1:mdls
    
    LauncherSS = LauncherSS_grid(:, :, ii);

    systemnames = 'ControllerSS LauncherSS ActuatorSS Delta';

    inputvar    = ' [wind]';
    outputvar   = ' [LauncherSS(1)] ';

    input_to_ActuatorSS     = ' [Delta(1)+ ControllerSS(1) + Delta(1)] ';
    input_to_Delta          = ' [ControllerSS(1) + Delta(1)] ';
    input_to_LauncherSS     = ' [ActuatorSS(1); wind(1)] ';
    input_to_ControllerSS   = ' [LauncherSS(2)] ';
    cleanupsysic = 'yes';    

    ClosedLoopSS_pert_LTI_WC = sysic;

%         plant =P(:,:,ii);
%         controller = C(:,:,ii);
%         
%         systemnames = 'plant controller delta tdelay';
%         inputvar = '[ d{1} ]';
%         outputvar = '[ -plant(1) ]';
%         input_to_plant =  '[tdelay+d(1)]';
%         input_to_tdelay ='[delta+controller]';
%         input_to_controller = '[-plant(1) ]';
%         input_to_delta = '[controller]';
%         sysoutname = 'sensii';
%         cleanupsysic = 'yes';
%         sysic
        
        [wcg{ii},wcu{ii},info] = wcgain(ClosedLoopSS_pert_LTI_WC);
        
        wcgjj(ii) = wcg{ii}.UpperBound;
        
    end
    
    %worst case gain for a given kappa over all grid points
    wcglti(jj+1) = max(wcgjj);
end


%
% for jj = 1:length(b)
%     kappa=b(jj);
%     clear delta
%     delta = ultidyn('delta',[1,1],'Bound',kappa);
%     for ii = 1:length(yvec)
%
%         plant =P(:,:,ii);
%         controller = C(:,:,ii);
%
%         systemnames = 'plant controller delta tdelay';
%         inputvar = '[ d{1} ]';
%         outputvar = '[ -plant(1) ]';
%         input_to_plant =  '[tdelay+d(1)]';
%         input_to_tdelay ='[delta+controller]';
%         input_to_controller = '[-plant(1) ]';
%         input_to_delta = '[controller]';
%         sysoutname = 'sensii';
%         cleanupsysic = 'yes';
%         sysic
%
%         [wcg{ii},wcu{ii},info] = wcgain(sensii);
%
%         wcgjj(ii) = wcg{ii}.UpperBound;
%
%     end
%
%     %worst case gain for a given kappa over all grid points
%     wcglti(jj+1) = max(wcgjj);
% end





