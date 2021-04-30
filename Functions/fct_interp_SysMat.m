function [P_rhs, P_lhs, CtransC, BXBtrans] = fct_interp_SysMat(t, SystemMatrices,M, M_IQC, SysSize, lambda, IQClist)

T_grid  = SystemMatrices.T_grid;
T_shift = T_grid; % is hard coded must be more flexible for different time spans of the LTV model is used to integrate to zero, but I think it's not necessary

nx = SysSize(1);
nw = SysSize(2);


nrIQC   = numel(SystemMatrices.C1);
C1_int{nrIQC}  = [];
D11_int{nrIQC} = [];
D12_int{nrIQC} = [];

for i = 1 : 1 : nrIQC

    M_IQC_int = fct_lin_interp(t , T_shift, M_IQC{i});
    
    nz = size(SystemMatrices.C1{i}, 1);
    
    
     C1_int{i}  = M_IQC_int(1 : nz, 1 : nx);
     D11_int{i} = M_IQC_int(1 : nz, (nx+1) : (nx + nw));
     D12_int{i} = M_IQC_int(1 : nz, (nx + nw + 1) : end);
    
end
    
% build the matrices like in my example for TUC

M_int   = fct_lin_interp(t , T_shift, M);
A_int   = M_int(1 : nx, 1 : nx);
B1_int  = M_int(1 : nx, (nx + 1) : (nx + nw));
B2_int  = M_int(1 : nx, (nx + nw + 1) : end);
C2_int  = M_int((nx + 1) : end, 1 : nx);
D21_int = M_int((nx + 1) : end, (nx + 1) : (nx + nw));
D22_int = M_int((nx + 1) : end, (nx + nw + 1) : end);

% can use them to build the new formulation

O = lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D11_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D11_int{2};
R = (lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D11_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D11_int{2})';

Q = lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D12_int{2};
S = (lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D12_int{2})';

U = -lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*C1_int{1} - lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*C1_int{2};

X = inv([(lambda(1)*D11_int{1}'*IQClist.IQC{1}{2}*D11_int{1} + lambda(2)*D11_int{2}'*IQClist.IQC{2}{2}*D11_int{2}),   (lambda(1)*D11_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*D11_int{2}'*IQClist.IQC{2}{2}*D12_int{2});...
         (lambda(1)*D12_int{1}'*IQClist.IQC{1}{2}*D11_int{1}+ lambda(2)*D12_int{2}'*IQClist.IQC{2}{2}*D11_int{2}),   (lambda(1)*D12_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*D12_int{2}'*IQClist.IQC{2}{2}*D12_int{2} - eye(1))]);
    
B = [B1_int, B2_int];

BXBtrans = -B * X * B';

P_rhs = B * X * [R ; S] - A_int;
P_lhs = (B * X * [R; S]- A_int)';

CtransC = [O, Q] * X * [R; S] + U;
    





end