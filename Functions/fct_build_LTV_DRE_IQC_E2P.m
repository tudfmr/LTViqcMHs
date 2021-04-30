function [P_rhs_t, P_lhs_t, CtransC_t, BXBtrans_t, C2_t] = fct_build_LTV_DRE_IQC_E2P(t, T_grid,P_rhs, BXBtrans, CTransC, C2)

T_shift = T_grid;

P_rhs_t     = fct_lin_interp(t , T_shift, P_rhs);
P_lhs_t     = P_rhs_t'; 
CtransC_t   = fct_lin_interp(t, T_shift, CTransC);
BXBtrans_t  = fct_lin_interp(t, T_shift, BXBtrans);
C2_t        = fct_lin_interp(t, T_shift, C2);

% nx = SysSize(1);
% nw = SysSize(2);
% nz = SysSize(4); % fuer andere parametrisierungen mal schauen
% 
% nrIQC   = numel(SystemMatrices.C1);
% C1_int{nrIQC}  = [];
% D11_int{nrIQC} = [];
% D12_int{nrIQC} = [];
% M_IQC_int{nrIQC} = [];
% 
% for i = 1 : 1 : nrIQC
%     
%     M_IQC_int{i} = fct_lin_interp(t , T_shift, M_IQC{i});
%     
%     %nz = size(SystemMatrices.C1{i}, 1);
%             
%     C1_int{i}  = M_IQC_int{i}(1 : nz, 1 : nx);
%     D11_int{i} = M_IQC_int{i}(1 : nz, (nx+1) : (nx + nw));
%     D12_int{i} = M_IQC_int{i}(1 : nz, (nx + nw + 1) : end);
%     
% end
% 
% % build the matrices as shown in the corresponding paper
% 
% M_int   = fct_lin_interp(t , T_shift, M);
% A_int   = M_int(1 : nx, 1 : nx);
% B1_int  = M_int(1 : nx, (nx + 1) : (nx + nw));
% B2_int  = M_int(1 : nx, (nx + nw + 1) : end);
% C2_int  = M_int((nx + 1) : end, 1 : nx);
% D21_int = M_int((nx + 1) : end, (nx + 1) : (nx + nw));
% D22_int = M_int((nx + 1) : end, (nx + nw + 1) : end);
% 
% % use them to build the new formulation
% 
% O = lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D11_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D11_int{2} + C2_int'*D21_int;
% %R = (lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D11_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D11_int{2} + C2_int'*D21_int)';
% 
% Q = lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D12_int{2} + C2_int'*D22_int;
% %S = (lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*D12_int{2} + C2_int'*D22_int)';
% 
% %U = -lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*C1_int{1} - lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*C1_int{2} - C2_int'*C2_int;
% 
% X = inv([(D21_int'*D21_int + lambda(1)*D11_int{1}'*IQClist.IQC{1}{2}*D11_int{1} + lambda(2)*D11_int{2}'*IQClist.IQC{2}{2}*D11_int{2}),   (D21_int'*D22_int + lambda(1)*D11_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*D11_int{2}'*IQClist.IQC{2}{2}*D12_int{2});...
%     (D22_int'*D21_int + lambda(1)*D12_int{1}'*IQClist.IQC{1}{2}*D11_int{1}+ lambda(2)*D12_int{2}'*IQClist.IQC{2}{2}*D11_int{2}),   (D22_int'*D22_int + lambda(1)*D12_int{1}'*IQClist.IQC{1}{2}*D12_int{1} + lambda(2)*D12_int{2}'*IQClist.IQC{2}{2}*D12_int{2} - eye(1)*gamma^2)]);
% 
% %B = [B1_int, B2_int];
% 
% %BXBtrans = -B * X * B';
% BXBtrans = -[B1_int, B2_int] * X * [B1_int, B2_int]';
% 
% P_rhs = [B1_int, B2_int] * X * [O' ; Q'] - A_int;
% %P_lhs = (B * X * [O'; Q']- A_int)';
% 
% %CtransC = [O, Q] * X * [O'; Q'] + U;
% CtransC = [O, Q] * X * [O'; Q'] -lambda(1)*C1_int{1}'*IQClist.IQC{1}{2}*C1_int{1} - lambda(2)*C1_int{2}'*IQClist.IQC{2}{2}*C1_int{2} - C2_int'*C2_int;

end