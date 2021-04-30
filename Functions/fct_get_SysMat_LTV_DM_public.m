function [SystemMatrices, SysSize] = fct_get_SysMat_LTV_DM_public(G,IQC, T_grid)


%% Dimensions
szG = size(G);
if numel(szG) == 2
    szG = [szG 1];
end
nmod = szG(3); % # of model in the G array.
nve = szG(1);  % # of outputs of G = nv + ne
nwd = szG(2);  % # of inputs of G = nw + nd

Niqc = size(IQC,2);    % # of IQCs
szPsi = size(IQC{1}{1},1);
nv = floor(szPsi/2);   % # of inputs to Delta
nw = nv;              % # of outputs of Delta (Delta assumed to be square)
nd = nwd-nw;   % # of disturbances
ne = nve-nv;   % # of errors


%% Compute hard factorization for each IQC multiplier
Psi = cell(Niqc,1);
M = cell(Niqc,1);
for i=1:Niqc
    
    Psi{i} = IQC{i}{1};
    M{i} = IQC{i}{2};
    
end

%%
% Construct extended system Pext including dynamics of plant G and
% all auxiliary systems Psi. The I/O relationship for Pext is:
%     [z1; ...; zNiqc; e] = Pext*[w; d]
% where zj is the output of the auxiliary system Psi for the j^th IQC.

% Stack all auxiliary systems:

%      [z1; ...; zNiqc] = AllPsi*[v; w]
% zidx is a cell array of the output indices for zi.
AllPsi = vertcat(Psi{:});
nz = size(AllPsi,1);

zptr = 0;
zidx = cell(Niqc,1);
for i=1:Niqc
    zidx{i} = zptr + (1:2*nv); % for nv = 2 --> 1:1:4 --> 1,2,3,4
    zptr = zptr + 2*nv;        % the point shifts through for nv = 1 --> 2 outputs first loop zidx= 1,2 second lop is zidx = 3,4
end                            %

% Form mapping [z;e] = AllPsiI*[v;w;e]
AllPsiI = blkdiag(AllPsi,eye(ne)); % keeps it more flexible so that the number of IQCs connected to the OL Plant doesn't matter

% Form mapping [v;w;e]=Gw*[w;d]
GI = [G; eye(nwd)]; % just adds new outputs ... pure feedthrough, depending on the number w + d (inputs of the system G)
vweidx = [1:nv, (nv+ne)+(1:nw), nv+(1:ne)];
Gw = GI(vweidx,:); % reorders the outputs keeps all inputs removes in d(1), v(1) and w(1) case one output
% - puts the initially first placed outputs associated with the v(i) in first position (1 : nv)
% - puts the initially second placed outputs associated with the w(i) (1 : nw) in the
%   last position (nv+ne --> after the v's and e's)
% - puts the initially third placed outputs associated with the e (1 : ne) in the
%   second position (after the v)
% -

% Form [z;e]=Pext*[w;d] from the product AllPsi*M ([w & d are the inputs of G and the ext. system])
% z (output psi) & e (output G) are the outputs of the ext. system

% all the reordering is just done to keep it as flexible as possible for
% all combination of input and output sizes

% Pext = ss( zeros(nz+ne,nwd,nmod) );
for i=1:nmod
    Gwi = Gw(:,:,i);
    Pext(:,:,i) = AllPsiI*Gwi;
end

% Partition state matrices of Pext based on inputs [w;d] and outputs [z;e]
[A,B,C,D] = ssdata(Pext);
nx = size(A,1);
Ce = C(end-ne+1:end,:,:);
De = D(end-ne+1:end,:,:);

%% ----------------- Build and evaluate the Inequality --------------------


for k1 = 1 : 1 : nmod
    
    Ak      = A(:, :, k1); % state matrix of the ext system
    Bk      = B(:,:,k1);
    Ck      = C(:, :, k1);
    Dk      = D(:, :, k1);
    Cek     = Ce(:,:,k1);
    Dek     = De(:,:,k1);
     
    Bwk     = Bk(:, (1 : nw));
    Bdk     = Bk(:, (nw +(1 : nd)));
    
    Dewk    = Dek(:, (1 : nw));
    Dedk    = Dek(:, (nw + (1 : nd)));
    
    for i = 1 : 1 : Niqc
    
        FieldName = ['IQC',num2str(i)];
        
        Cwk1.(FieldName)   = Ck(zidx{i}, :);
        Dwk11.(FieldName)  = Dk(zidx{i}, (1 : nw));
        Dwk12.(FieldName)  = Dk(zidx{i}, (nw +(1 : nd)));
 
    end
 

    
    SystemMatrices.A(:, :, k1)      = Ak;
    SystemMatrices.B1(:, :, k1)     = Bwk;
    SystemMatrices.B2(:, :, k1)     = Bdk;
    
    SystemMatrices.C2(:, :, k1)     = Cek;
    SystemMatrices.D21(:, :, k1)    = Dewk;
    SystemMatrices.D22(:, :, k1)    = Dedk;
    
    for i = 1 : 1 : Niqc
    
        FieldName = ['IQC',num2str(i)];
        
        SystemMatrices.C1{i}(:, :, k1)     = Cwk1.(FieldName);
        SystemMatrices.D11{i}(:, :, k1)    = Dwk11.(FieldName);
        SystemMatrices.D12{i}(:, :, k1)    = Dwk12.(FieldName);
        
        BlkMat_IQC{i}(:, :, k1) = [Cwk1.(FieldName), Dwk11.(FieldName), Dwk12.(FieldName)];
    
    end
    
    
    
end

SystemMatrices.T_grid = T_grid;

% allcocate all the non IQC terms/ matrices in one matrix for the later interpolation

BlkMat = [[SystemMatrices.A; SystemMatrices.C2], [[SystemMatrices.B1, SystemMatrices.B2]; [SystemMatrices.D21, SystemMatrices.D22]]];

SysSize(1) = size(SystemMatrices.A, 1);     % number of x
SysSize(2) = size(SystemMatrices.B1, 2);    % number of w
SysSize(3) = size(SystemMatrices.B2, 2);    % number of d
SysSize(4) = size(SystemMatrices.C1{1}, 1); % number of z1 (make it more flexible later on to get the results)
SysSize(5) = size(SystemMatrices.C2, 1);    % number of e

% Allocate all the IQC terms/ matrices in one matrix fot the later
% interpolation






