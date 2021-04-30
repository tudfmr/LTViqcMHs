function [P_rhs, BXBtrans, CTransC] = fct_build_RDE_components_08052018(SystemMatrices, IQClist, gamma, lambda)

nModels = size(SystemMatrices.A, 3);


for i = 1 : 1 : nModels

    O(:, :, i) = lambda(1)*SystemMatrices.C1{1}(:, :, i)'*IQClist.IQC{1}{2}*SystemMatrices.D11{1}(:, :, i) + lambda(2)*SystemMatrices.C1{2}(:, :, i)'*IQClist.IQC{2}{2}*SystemMatrices.D11{2}(:, :, i) + SystemMatrices.C2(:, :, i)'*SystemMatrices.D21(:, :, i);
    Q(:, :, i) = lambda(1)*SystemMatrices.C1{1}(:, :, i)'*IQClist.IQC{1}{2}*SystemMatrices.D12{1}(:, :, i) + lambda(2)*SystemMatrices.C1{2}(:, :, i)'*IQClist.IQC{2}{2}*SystemMatrices.D12{2}(:, :, i) + SystemMatrices.C2(:, :, i)'*SystemMatrices.D22(:, :, i);
    
    X(:, :, i) = inv([(SystemMatrices.D21(:, :, i)'*SystemMatrices.D21(:, :, i) + lambda(1)*SystemMatrices.D11{1}(:, :, i)'*IQClist.IQC{1}{2}*SystemMatrices.D11{1}(:, :, i) + lambda(2)*SystemMatrices.D11{2}(:, :, i)'*IQClist.IQC{2}{2}*SystemMatrices.D11{2}(:, :, i)), (SystemMatrices.D21(:, :, i)'*SystemMatrices.D22(:, :, i) + lambda(1)*SystemMatrices.D11{1}(:, :, i)'*IQClist.IQC{1}{2}*SystemMatrices.D12{1}(:, :, i) + lambda(2)*SystemMatrices.D11{2}(:, :, i)'*IQClist.IQC{2}{2}*SystemMatrices.D12{2}(:, :, i));...
                      (SystemMatrices.D22(:, :, i)'*SystemMatrices.D21(:, :, i) + lambda(1)*SystemMatrices.D12{1}(:, :, i)'*IQClist.IQC{1}{2}*SystemMatrices.D11{1}(:, :, i) + lambda(2)*SystemMatrices.D12{2}(:, :, i)'*IQClist.IQC{2}{2}*SystemMatrices.D11{2}(:, :, i)), (SystemMatrices.D22(:, :, i)'*SystemMatrices.D22(:, :, i) + lambda(1)*SystemMatrices.D12{1}(:, :, i)'*IQClist.IQC{1}{2}*SystemMatrices.D12{1}(:, :, i) + lambda(2)*SystemMatrices.D12{2}(:, :, i)'*IQClist.IQC{2}{2}*SystemMatrices.D12{2}(:, :, i) - eye(1)*gamma^2)]);

    
    BXBtrans(:, :, i) = -[SystemMatrices.B1(:, :, i), SystemMatrices.B2(:, :, i)] * X(:, :, i) * [SystemMatrices.B1(:, :, i), SystemMatrices.B2(:, :, i)]';
    P_rhs(:, :, i)  = [SystemMatrices.B1(:, :, i), SystemMatrices.B2(:, :, i)] * X(:, :, i) * [O(:, :, i)' ; Q(:, :, i)'] - SystemMatrices.A(:, :, i);
%   P_lhs(:, :, i)  = P_rhs(:, :, i)';
    CTransC(:, :, i) = [O(:, :, i), Q(:, :, i)] * X(:, :, i) * [O(:, :, i)' ; Q(:, :, i)'] -lambda(1)*SystemMatrices.C1{1}(:, :, i)'*IQClist.IQC{1}{2}*SystemMatrices.C1{1}(:, :, i) - lambda(2)*SystemMatrices.C1{2}(:, :, i)'*IQClist.IQC{2}{2}*SystemMatrices.C1{2}(:, :, i) - SystemMatrices.C2(:, :, i)'*SystemMatrices.C2(:, :, i);

   
end
