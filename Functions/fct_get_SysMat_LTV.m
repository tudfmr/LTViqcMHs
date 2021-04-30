function [A, B, C, D, M] = fct_get_SysMat_LTV(P_rho_CL)

A = P_rho_CL.A;
B = P_rho_CL.B;
C = P_rho_CL.C;
D = P_rho_CL.D;

%T_grid = P_rho_CL.(P_rho_id).Parameter.T.GridData; % gets the time grid -> x data for the interpolation

nx      = size(A, 1);
nd      = size(B, 2);
ny      = size(C, 1);
nmdl    = size(A, 3);

M = zeros((nx+ny), (nx+nd), nmdl);
M(1:nx, 1:nx, :) = A;
M(1:nx, (nx+1) : end, :) = B;
M((nx+1) : end, 1 : nx, :) = C;
M((nx+1) : end, (nx+1) : end, :) = D;



end