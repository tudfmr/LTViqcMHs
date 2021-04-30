function [A_int, B_int, C_int, D_int] = fct_build_LTV_DARE_nom(t, T_grid, A, B, C, D, M,nx, nd, ny,nmdl)

% Get system matrices

% A = P_rho_CL.(P_rho_id).Data.A;
% B = P_rho_CL.(P_rho_id).Data.B;
% C = P_rho_CL.(P_rho_id).Data.C;
% D = P_rho_CL.(P_rho_id).Data.D;

% T = P_rho_CL.(P_rho_id).Parameter.T.GridData; % gets the time grid -> x data for the interpolation

T_shift = T_grid; % is hard coded must be more flexible for different time spans of the LTV model is used to integrate to zero, but I think it's not necessary

% permute to get interpolated dimension first (it's the 3rd related to the time steps of the LTV)

% % A_per = permute(A, [3 1 2]);
% % B_per = permute(B, [3 1 2]);
% % C_per = permute(C, [3 1 2]);
% % D_per = permute(D, [3 1 2]);
% % 
% % %interpolate for time t and permute to the correct dimensions
% % 
% % A_int = permute(interp1(T_shift, A_per, t, 'linear'), [2, 3, 1]);
% % B_int = permute(interp1(T_shift, B_per, t, 'linear'), [2, 3, 1]);
% % C_int = permute(interp1(T_shift, C_per, t, 'linear'), [2, 3, 1]);
% % D_int = permute(interp1(T_shift, D_per, t, 'linear'), [2, 3, 1]);


% A_int =  fct_lin_interp(t , T_shift, A);
% B_int =  fct_lin_interp(t , T_shift, B);
% C_int =  fct_lin_interp(t , T_shift, C);
% D_int =  fct_lin_interp(t , T_shift, D);

M_int = fct_lin_interp(t , T_shift, M);
%M_int = fct_interp1qr_1(t , T_shift, M);

A_int = M_int(1:nx, 1:nx);
B_int = M_int(1:nx, (nx+1) : end);
C_int = M_int((nx+1) : end, 1 : nx);
D_int = M_int((nx+1) : end, (nx+1) : end);



end