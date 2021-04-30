function   [M_int] = fct_lin_interp(t , t_grid, M_grid)

% find the position in the Grid

idx_min = find(t_grid <= t, 1, 'last');
idx_max = find(t_grid >= t, 1, 'first');

M_min = M_grid(:, :, idx_min);
M_max = M_grid(:, :, idx_max);

t_min = t_grid(idx_min);
t_max = t_grid(idx_max);

dMdt = (M_max - M_min)/(t_max - t_min);

dMdt(isnan(dMdt)) = 0;
%test1(isnan(test1))= 0

M_int = dMdt*(t - t_min) + M_min;
