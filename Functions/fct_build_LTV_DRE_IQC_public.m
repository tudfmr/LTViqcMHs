function [P_rhs_t, P_lhs_t, CtransC_t, BXBtrans_t] = fct_build_LTV_DRE_IQC_public(t, T_grid,P_rhs, BXBtrans, CTransC)

T_shift = T_grid;

P_rhs_t     = fct_lin_interp(t , T_shift, P_rhs);
P_lhs_t     = P_rhs_t'; 
CtransC_t   = fct_lin_interp(t, T_shift, CTransC);
BXBtrans_t  = fct_lin_interp(t, T_shift, BXBtrans);

end