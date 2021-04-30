
% under construction
% these are actually just quick shots
% proper implementation follows

function [g,info] = fct_calc_wc_L2_LTV_LMI(G,IQC,Fbasis,Fgrad)

% function [g,info] =  simplenormboundedgain(G,PI,Fbasis,Fgrad,RateBounds)
%
% Compute bound on gain of Fu(G,Delta) where G is an LTI system and Delta
% is an operator described by IQCs PI.
%
% Inputs
%  G: Linear parameter varying system.
%  IQC: Cell array containing the IQC factorization {Psi_i,M_i} for
%       the i^th multiplier for Delta. This function assumes Delta
%       is square m-by-m.
% Fbasis: Basis functions used in storage function Nbasis-by-1-by-nmod
% Fgrad: Gradient of the basis functions wrt. the parameters:
%        Nbasis-by-Nparameters-by-nmod
% RateBounds: Upper and lower bounds on the parameter rates: Nparameter-by-2

% Outputs
%  g: Bound on the gain of Fu(G,Delta)
%  info:  Structure with information regarding solution of dissipation
%         inequality. Structure includes the fields:
%    -P:  Symmetric matrix used to construct the storage function.
%    -lamopt: Coefficients for optimal conic combination of IQC multipliers
%

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


nbasis = size(Fbasis,1); % # of basis functions


% Check for non-rate bounded case
ratebndflg = true;
if (nbasis==1 && npar==0)
    ratebndflg = false;
end

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
    zidx{i} = zptr + (1:2*nv);
    zptr = zptr + 2*nv;
end

% Form mapping [z;e] = AllPsiI*[v;w;e]
AllPsiI = blkdiag(AllPsi,eye(ne));

% Form mapping [v;w;e]=Gw*[w;d]
GI = [G; eye(nwd)];
vweidx = [1:nv, (nv+ne)+(1:nw), nv+(1:ne)];
Gw = GI(vweidx,:);

% Form [z;e]=Pext*[w;d] from the product AllPsi*M
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

%%
% Initialize LMI and define variables
setlmis([])
X = zeros(nbasis,1);
for i1 = 1:nbasis
    X(i1) = lmivar(1,[nx 1]);
end

for i2 = 1:Niqc
    [lam(i2),nlam(i2)] = lmivar(1,[1 1]);
end
[gsq,ndec] = lmivar(1,[1 1]);
gsqI = lmivar(3,blkdiag(zeros(nx+nw),ndec*eye(nd)));

%%
% Specify LMIs

% LMI for Dissipation Inequality - Loop through array of models
cnt = 1;
for k1 = 1:nmod
    Ak = A(:,:,k1);
    Bk = B(:,:,k1);
    Cek = Ce(:,:,k1);
    Dek = De(:,:,k1);
    
    if ratebndflg
        Fgk1 = Fgrad(:,:,k1);
    else
        Fgk1 = [];
    end
    Fbk1 = Fbasis(:,:,k1);
    
    % L2 Bound LMI in block 3-by-3 form with gamma
    
    % der for loop kann raus un der cnt = cnt + 1 auch --> da das nur noch
    % eine einzelne LMI ist
    
    
    % ---- the original deltaP part
    
    %for q = 1:2^npar  % number of variables appearing in X's basis fcns
    
    for ibasis=1:nbasis
        lmiterm([cnt 1 1 X(ibasis)],[eye(nx); zeros(nw+nd,nx)],[Ak Bk]*Fbk1(ibasis,1),'s');
    end
    
    %L2 performance term
    lmiterm([-cnt 1 1 gsqI],1,1);
    lmiterm([cnt 1 1 0],[Cek Dek]'*[Cek Dek]);
    for i1=1:Niqc
        
        % Grab i1^th IQC
        Mi = M{i1};
        Cyi = C(zidx{i1},:,k1);
        Dyi = D(zidx{i1},:,k1);
        Mtil = [Cyi Dyi]'*Mi*[Cyi Dyi];
        % lmiterm([cnt 1 1 lam(i1)],Mtil,1,'s'); % wouldn't this be 2 times
        lmiterm([cnt 1 1 lam(i1)],Mtil,1);
        
    end
    
    %if ratebndflg
    %rbvec = RateBounds(:,1);
    %tmp = dec2bin(q-1,npar);
    %idx = find(tmp=='1');
    %rbvec(idx) = RateBounds(idx,2);
    for k2=1:nbasis
        
        % this term would mean that the P_dot are on the rhs therefore
        % negative on the lhs --> doesn't work
        %         lmiterm([-cnt 1 1 X(k2)],[eye(nx); zeros(nwd,nx)],...
        %             0.5*Fgk1(k2, 1)*...                             % pure time differentiation
        %             [eye(nx) zeros(nx,nwd)],'s');
        
        lmiterm([cnt 1 1 X(k2)],[eye(nx); zeros(nwd,nx)],...
                0.5*Fgk1(k2, 1)*...                             % pure time differentiation
                [eye(nx) zeros(nx,nwd)],'s');
        
        
    end
    %end
    
    
    %end
    
    cnt = cnt+1;
end

% change to the
% LMI to enforce P(T)> = 0 for end time!

if ratebndflg         % || k1 == 1 % check for what this is needed
    for n1=1:nbasis
        lmiterm([-cnt 1 1 X(n1)],Fbasis(n1, :, nmod),1);
    end
    cnt = cnt+1;
end


% LMI to enforce scalings to be non-negative
for i1 = 1:Niqc
    
    lmiterm([-(cnt) 1 1 lam(i1)],1,1);
    cnt = cnt+1;
    
end

%% Solve SDP: min gamsq subject to LMI constraints
lmisys = getlmis;
c = zeros(ndec,1);
c(end) = 1;

% opt = [0 0 0 0 1]; % org
opt = [0 0 0 0 0];
[copt,xopt] = mincx(lmisys,c,opt);

%% Store solution
g = [];
info.P = [];
info.lamopt = [];
if ~isempty(xopt)
    
    evals = evallmi(lmisys, xopt);
    
    info.eval.lmisys = lmisys;
    info.eval.xopt   = xopt;
    
    info.eval.lhs = zeros((nx+nw+nd), (nx+nw+nd), nmod);
    info.eval.lhs = zeros((nx+nw+nd), (nx+nw+nd), nmod);
    
    info.eval.lambda = zeros(1,1, Niqc);
    info.eval.P_T      = zeros(nx, nx, 1);
    
    for i  = 1 : 1 : (cnt-1)
        
        [lhs, rhs] = showlmi(evals,i);
        
        if isequal(size(rhs, 1), size(info.eval.lhs, 1))
            
            [info.eval.lhs(:, :, i), info.eval.rhs(:, :, i)] = showlmi(evals, i);
            
            info.eval.lhs(:, :, i) = lhs;
            info.eval.rhs(:, :, i) = rhs;
            
        elseif isequal(size(rhs, 1), size(info.eval.lambda, 1))
            
            info.eval.lambda(:, :, i) = rhs;
            
        elseif isequal(size(rhs, 1), size(info.eval.P_T, 1))
            
            info.eval.P_T(:, :, i) = rhs;
            
        end
    end
    
    
    g = sqrt(copt);
    for i1 = 1:nbasis
        info.P(:,:,i1) =  dec2mat(lmisys,xopt,X(i1));
    end
    lamopt = zeros(Niqc,1);
    for i = 1:Niqc
        lamopt(i) = dec2mat(lmisys,xopt,lam(i));
    end
    info.lamopt = lamopt;
    
elseif isempty(xopt)
    
    info.P  = zeros(nx, nx, nbasis);
    
    lamopt = zeros(Niqc,1);
    
    info.lamopt = lamopt;
    
    info.eval = [];
    
end
