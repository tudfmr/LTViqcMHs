function P_rho_CL = fct_build_LPV_norm_CL(G_rho, K_rho, G_rho_id, K_rho_id, P_rho_id, Settings)


if isequal(Settings.IOgain, 'psi_ref2psi_err')
    
    G_rho_temp = G_rho.(G_rho_id);
    K_rho_temp = K_rho.(K_rho_id);
    
    systemnames      = 'G_rho_temp K_rho_temp';
    
    inputvar         = '[ d{1} ]';   % w - output of the uncertainty
    % e - performance output for L2 gain
    
    outputvar        = '[ d(1) - G_rho_temp(1) ]';
    
    input_to_G_rho_temp   = '[ K_rho_temp ]';
    input_to_K_rho_temp   = '[ d(1) - G_rho_temp(1); -G_rho_temp(2); -G_rho_temp(3) ]';
    
    
    sysoutname       = 'P_rho_CL_temp';
    cleanupsysic     = 'yes';
    sysic
    
    % P_rho_id = [G_rho_id, '_', K_rho_id];
    
    P_rho_CL.(P_rho_id) = P_rho_CL_temp;
    
    
elseif isequal(Settings.IOgain, 'zeta_ref2zeta_err')
    
    if isequal(Settings.IncludeWeight, 0)
        
        G_rho_temp = G_rho.(G_rho_id);
        K_rho_temp = K_rho.(K_rho_id);
        
        systemnames      = 'G_rho_temp K_rho_temp';
        
        inputvar         = '[ d{1} ]';   % w - output of the uncertainty
        % e - performance output for L2 gain
        
        outputvar        = '[ d(1) - G_rho_temp(2) ]';
        
        input_to_G_rho_temp   = '[ K_rho_temp ]';
        input_to_K_rho_temp   = '[ -G_rho_temp(1); d(1) - G_rho_temp(2); -G_rho_temp(3) ]';
        
        
        sysoutname       = 'P_rho_CL_temp';
        cleanupsysic     = 'yes';
        sysic
        
        % P_rho_id = [G_rho_id, '_', K_rho_id];
        
        P_rho_CL.(P_rho_id) = P_rho_CL_temp;
        
    else
        
        G_rho_temp = G_rho.(G_rho_id);
        K_rho_temp = K_rho.(K_rho_id);
        W_P        = ss(tf([1/1.5, 0.0481], [1 0.0481 * 10^-3]));
        systemnames      = 'G_rho_temp K_rho_temp W_P';
        
        inputvar         = '[ d{1} ]';   % w - output of the uncertainty
        % e - performance output for L2 gain
        
        outputvar        = '[ W_P ]';
        
        input_to_W_P          = '[ d(1) - G_rho_temp(2) ]';
        input_to_G_rho_temp   = '[ K_rho_temp ]';
        input_to_K_rho_temp   = '[ -G_rho_temp(1); d(1) - G_rho_temp(2); -G_rho_temp(3) ]';
        
        
        sysoutname       = 'P_rho_CL_temp';
        cleanupsysic     = 'yes';
        sysic
        
        % P_rho_id = [G_rho_id, '_', K_rho_id];
        
        P_rho_CL.(P_rho_id) = P_rho_CL_temp;
        
        
    end
    
elseif isequal(Settings.IOgain, 'psi_ref2Qalpha')
    
    
    G_rho_temp = G_rho.(G_rho_id);
    K_rho_temp = K_rho.(K_rho_id);
    
    systemnames      = 'G_rho_temp K_rho_temp';
    
    inputvar         = '[ d{1} ]';   % w - output of the uncertainty
    % e - performance output for L2 gain
    
    outputvar        = '[ G_rho_temp(4) ]';
    
    input_to_G_rho_temp   = '[ K_rho_temp ]';
    input_to_K_rho_temp   = '[ d(1) - G_rho_temp(1); -G_rho_temp(2); -G_rho_temp(3) ]';
    
    
    sysoutname       = 'P_rho_CL_temp';
    cleanupsysic     = 'yes';
    sysic
    
    % P_rho_id = [G_rho_id, '_', K_rho_id];
    
    P_rho_CL.(P_rho_id) = P_rho_CL_temp;
    
elseif isequal(Settings.IOgain, 'zeta_ref2Qalpha')
    
    G_rho_temp = G_rho.(G_rho_id);
    K_rho_temp = K_rho.(K_rho_id);
    
    systemnames      = 'G_rho_temp K_rho_temp';
    
    inputvar         = '[ d{1} ]';   % w - output of the uncertainty
    % e - performance output for L2 gain
    
    outputvar        = '[ G_rho_temp(4) ]';
    
    input_to_G_rho_temp   = '[ K_rho_temp ]';
    input_to_K_rho_temp   = '[ -G_rho_temp(1); d(1)-G_rho_temp(2); -G_rho_temp(3) ]';
    
    
    sysoutname       = 'P_rho_CL_temp';
    cleanupsysic     = 'yes';
    sysic
    
    % P_rho_id = [G_rho_id, '_', K_rho_id];
    
    P_rho_CL.(P_rho_id) = P_rho_CL_temp;
    
elseif isequal(Settings.IOgain, 'psi_ref2zeta_dot_err')
    
    G_rho_temp = G_rho.(G_rho_id);
    K_rho_temp = K_rho.(K_rho_id);
    
    systemnames      = 'G_rho_temp K_rho_temp';
    
    inputvar         = '[ d{1} ]';   % w - output of the uncertainty
    % e - performance output for L2 gain
    
    outputvar        = '[ -G_rho_temp(3) ]';
    
    input_to_G_rho_temp   = '[ K_rho_temp ]';
    input_to_K_rho_temp   = '[ d(1) - G_rho_temp(1); -G_rho_temp(2); -G_rho_temp(3) ]';
    
    
    sysoutname       = 'P_rho_CL_temp';
    cleanupsysic     = 'yes';
    sysic
    
    % P_rho_id = [G_rho_id, '_', K_rho_id];
    
    P_rho_CL.(P_rho_id) = P_rho_CL_temp;
    
    
    elseif isequal(Settings.IOgain, 'psi_ref2zeta_err')
    
    G_rho_temp = G_rho.(G_rho_id);
    K_rho_temp = K_rho.(K_rho_id);
    
    systemnames      = 'G_rho_temp K_rho_temp';
    
    inputvar         = '[ d{1} ]';   % w - output of the uncertainty
    % e - performance output for L2 gain
    
    outputvar        = '[ -G_rho_temp(2) ]';
    
    input_to_G_rho_temp   = '[ K_rho_temp ]';
    input_to_K_rho_temp   = '[ d(1) - G_rho_temp(1); -G_rho_temp(2); -G_rho_temp(3) ]';
    
    
    sysoutname       = 'P_rho_CL_temp';
    cleanupsysic     = 'yes';
    sysic
    
    % P_rho_id = [G_rho_id, '_', K_rho_id];
    
    P_rho_CL.(P_rho_id) = P_rho_CL_temp;
    
    
else
    
    error(['This in-/ output combination (',Settings.IOgain ,') is not implemented yet', char(10)])
    
end


end
