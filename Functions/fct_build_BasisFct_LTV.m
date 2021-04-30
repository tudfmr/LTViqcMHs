function [Fbasis, Fgrad] = fct_build_BasisFct_LTV(P_rho, IQClist, FBasisOrder,t_grid)

NrGridPts = length(P_rho);

for i = 1 : 1 : NrGridPts % das sollte sich ueber vertcat() loesen lassen
    
    ParamVec(i) = t_grid(i); % geht noch flexibler fuer beliebige Parameter
    
end


if isequal(FBasisOrder, 1)
    
    % ----------------------- Linear Basis Function -----------------------
    % rho enters affine
    
    % get the scheduling/ variable parameter vector
    
    
    % basis function for the parameter dependent lyapunov function
    % --> P affine defined as follows:
    %           P(rho) = P_0 + rho * P_1
    
    Fbasis = ones(2,1,NrGridPts); %First term is constant
    Fbasis(2,1,:) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2,1,:) = ones(1,1,NrGridPts);
    
    Fbasis_id = 'LinBasisFct';
    
    %
    %         T = pgrid('T', T);
    %
    %         b0      = basis(1, 0);
    %         b1      = basis(T, 1);
    %         Basis   = [b0, b1];
    
elseif isequal(FBasisOrder, 2)
    
    % ------------------- Quadratic Basis Function --------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2,1,:) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3,1,:) = reshape(ParamVec'.^2,1,1,NrGridPts); %third term quadratic in rho
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2,1,:) = ones(1,1,NrGridPts);
    Fgrad(3,1,:) = reshape(2*ParamVec',1,1,NrGridPts);
    
    Fbasis_id = 'QuadBasisFct';
    
    
    
elseif isequal(FBasisOrder, 3)
    
    % ---------------------- Cubic Basis Function ---------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2,1,:) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3,1,:) = reshape(ParamVec'.^2,1,1,NrGridPts); %third term quadratic in rho
    Fbasis(4,1,:) = reshape(ParamVec'.^3,1,1,NrGridPts); %third term cubic in rho
    
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2,1,:) = ones(1,1,NrGridPts);
    Fgrad(3,1,:) = reshape(2*ParamVec',1,1,NrGridPts);
    Fgrad(4,1,:) = reshape((3*(ParamVec.^2))',1,1,NrGridPts);
    
    Fbasis_id = 'CubicBasisFct';
    
elseif isequal(FBasisOrder, 4)
    
    % ---------------------- Cubic Basis Function ---------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2, 1, :) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3, 1, :) = reshape(ParamVec'.^2,1,1,NrGridPts); % third term quadratic in rho
    Fbasis(4, 1, :) = reshape(ParamVec'.^3,1,1,NrGridPts); % fourth term cubic in rho
    Fbasis(5, 1, :) = reshape(ParamVec'.^4,1,1,NrGridPts); % fith term 4 order in rho
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2,1,:) = ones(1,1,NrGridPts);
    Fgrad(3,1,:) = reshape(2*ParamVec',1,1,NrGridPts);
    Fgrad(4,1,:) = reshape((3*(ParamVec.^2))',1,1,NrGridPts);
    Fgrad(5,1,:) = reshape((4*(ParamVec.^3))',1,1,NrGridPts);
    
    
    Fbasis_id = 'Ord4BasisFct';
    
    
elseif isequal(FBasisOrder, 5)
    
    % ---------------------- Cubic Basis Function ---------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2, 1, :) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3, 1, :) = reshape(ParamVec'.^2,1,1,NrGridPts); % third term quadratic in rho
    Fbasis(4, 1, :) = reshape(ParamVec'.^3,1,1,NrGridPts); % fourth term cubic in rho
    Fbasis(5, 1, :) = reshape(ParamVec'.^4,1,1,NrGridPts); % fith term 4 order in rho
    Fbasis(6, 1, :) = reshape(ParamVec'.^5,1,1,NrGridPts); % sixth term 5 order in rho
    
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2,1,:) = ones(1,1,NrGridPts);
    Fgrad(3,1,:) = reshape(2*ParamVec',1,1,NrGridPts);
    Fgrad(4,1,:) = reshape((3*(ParamVec.^2))',1,1,NrGridPts);
    Fgrad(5,1,:) = reshape((4*(ParamVec.^3))',1,1,NrGridPts);
    Fgrad(6,1,:) = reshape((5*(ParamVec.^4))',1,1,NrGridPts);
    
    
    Fbasis_id = 'Ord5BasisFct';
    
    
elseif isequal(FBasisOrder, 6)
    
    % ---------------------- Cubic Basis Function ---------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2, 1, :) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3, 1, :) = reshape(ParamVec'.^2,1,1,NrGridPts); % third term quadratic in rho
    Fbasis(4, 1, :) = reshape(ParamVec'.^3,1,1,NrGridPts); % fourth term cubic in rho
    Fbasis(5, 1, :) = reshape(ParamVec'.^4,1,1,NrGridPts); % fith term 4 order in rho
    Fbasis(6, 1, :) = reshape(ParamVec'.^5,1,1,NrGridPts); % sixth term 5 order in rho
    Fbasis(7, 1, :) = reshape(ParamVec'.^6,1,1,NrGridPts); % fith term 6 order in rho
    
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2,1,:) = ones(1,1,NrGridPts);
    Fgrad(3,1,:) = reshape(2*ParamVec',1,1,NrGridPts);
    Fgrad(4,1,:) = reshape((3*(ParamVec.^2))',1,1,NrGridPts);
    Fgrad(5,1,:) = reshape((4*(ParamVec.^3))',1,1,NrGridPts);
    Fgrad(6,1,:) = reshape((5*(ParamVec.^4))',1,1,NrGridPts);
    Fgrad(7,1,:) = reshape((6*(ParamVec.^5))',1,1,NrGridPts);
    
    
    Fbasis_id = 'Ord6BasisFct';
    
    
elseif isequal(FBasisOrder, 7)
    
    % ---------------------- Cubic Basis Function ---------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2, 1, :) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3, 1, :) = reshape(ParamVec'.^2,1,1,NrGridPts); % third term quadratic in rho
    Fbasis(4, 1, :) = reshape(ParamVec'.^3,1,1,NrGridPts); % fourth term cubic in rho
    Fbasis(5, 1, :) = reshape(ParamVec'.^4,1,1,NrGridPts); % fith term 4 order in rho
    Fbasis(6, 1, :) = reshape(ParamVec'.^5,1,1,NrGridPts); % sixth term 5 order in rho
    Fbasis(7, 1, :) = reshape(ParamVec'.^6,1,1,NrGridPts); % fith term 6 order in rho
    Fbasis(8, 1, :) = reshape(ParamVec'.^7,1,1,NrGridPts); % fith term 7 order in rho
    
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2,1,:) = ones(1,1,NrGridPts);
    Fgrad(3,1,:) = reshape(2*ParamVec',1,1,NrGridPts);
    Fgrad(4,1,:) = reshape((3*(ParamVec.^2))',1,1,NrGridPts);
    Fgrad(5,1,:) = reshape((4*(ParamVec.^3))',1,1,NrGridPts);
    Fgrad(6,1,:) = reshape((5*(ParamVec.^4))',1,1,NrGridPts);
    Fgrad(7,1,:) = reshape((6*(ParamVec.^5))',1,1,NrGridPts);
    Fgrad(8,1,:) = reshape((7*(ParamVec.^6))',1,1,NrGridPts);
    
    
    
    
    Fbasis_id = 'Ord7BasisFct';
    
    
elseif isequal(FBasisOrder, 8)
    
    % ---------------------- Cubic Basis Function ---------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2, 1, :) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3, 1, :) = reshape(ParamVec'.^2,1,1,NrGridPts); % third term quadratic in rho
    Fbasis(4, 1, :) = reshape(ParamVec'.^3,1,1,NrGridPts); % fourth term cubic in rho
    Fbasis(5, 1, :) = reshape(ParamVec'.^4,1,1,NrGridPts); % fith term 4 order in rho
    Fbasis(6, 1, :) = reshape(ParamVec'.^5,1,1,NrGridPts); % sixth term 5 order in rho
    Fbasis(7, 1, :) = reshape(ParamVec'.^6,1,1,NrGridPts); % fith term 6 order in rho
    Fbasis(8, 1, :) = reshape(ParamVec'.^7,1,1,NrGridPts); % fith term 7 order in rho
    Fbasis(9, 1, :) = reshape(ParamVec'.^8,1,1,NrGridPts); % fith term 8 order in rho
    
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2, 1, :) = ones(1,1,NrGridPts);
    Fgrad(3, 1, :) = reshape(2*ParamVec',1,1,NrGridPts);
    Fgrad(4, 1, :) = reshape((3*(ParamVec.^2))',1,1,NrGridPts);
    Fgrad(5, 1, :) = reshape((4*(ParamVec.^3))',1,1,NrGridPts);
    Fgrad(6, 1, :) = reshape((5*(ParamVec.^4))',1,1,NrGridPts);
    Fgrad(7, 1, :) = reshape((6*(ParamVec.^5))',1,1,NrGridPts);
    Fgrad(8, 1, :) = reshape((7*(ParamVec.^6))',1,1,NrGridPts);
    Fgrad(9, 1, :) = reshape((8*(ParamVec.^7))',1,1,NrGridPts);
    
    
    
    
    Fbasis_id = 'Ord8BasisFct';
    
    
    
elseif isequal(FBasisOrder, 9)
    
    % ---------------------- Cubic Basis Function ---------------------
    
    % basis function for the parameter dependent lyapunov function
    
    Fbasis = ones(3,1,NrGridPts); %First term is constant
    Fbasis(2, 1, :) = reshape(ParamVec',1,1,NrGridPts); %second term linear in rho
    Fbasis(3, 1, :) = reshape(ParamVec'.^2,1,1,NrGridPts); % third term quadratic in rho
    Fbasis(4, 1, :) = reshape(ParamVec'.^3,1,1,NrGridPts); % fourth term cubic in rho
    Fbasis(5, 1, :) = reshape(ParamVec'.^4,1,1,NrGridPts); % fith term 4 order in rho
    Fbasis(6, 1, :) = reshape(ParamVec'.^5,1,1,NrGridPts); % sixth term 5 order in rho
    Fbasis(7, 1, :) = reshape(ParamVec'.^6,1,1,NrGridPts); % fith term 6 order in rho
    Fbasis(8, 1, :) = reshape(ParamVec'.^7,1,1,NrGridPts); % fith term 7 order in rho
    Fbasis(9, 1, :) = reshape(ParamVec'.^8,1,1,NrGridPts); % fith term 8 order in rho
    Fbasis(10, 1, :) = reshape(ParamVec'.^9,1,1,NrGridPts); % fith term 9 order in rho
    
    
    %gradient of the basis function with respect to the parameters
    
    Fgrad = zeros(2,1,NrGridPts);
    Fgrad(2, 1, :) = ones(1,1,NrGridPts);
    Fgrad(3, 1, :) = reshape(2*ParamVec',1,1,NrGridPts);
    Fgrad(4, 1, :) = reshape((3*(ParamVec.^2))',1,1,NrGridPts);
    Fgrad(5, 1, :) = reshape((4*(ParamVec.^3))',1,1,NrGridPts);
    Fgrad(6, 1, :) = reshape((5*(ParamVec.^4))',1,1,NrGridPts);
    Fgrad(7, 1, :) = reshape((6*(ParamVec.^5))',1,1,NrGridPts);
    Fgrad(8, 1, :) = reshape((7*(ParamVec.^6))',1,1,NrGridPts);
    Fgrad(9, 1, :) = reshape((8*(ParamVec.^7))',1,1,NrGridPts);
    Fgrad(10, 1, :) = reshape((9*(ParamVec.^8))',1,1,NrGridPts);
    
    
    
    
    Fbasis_id = 'Ord9BasisFct';
    
    
end



