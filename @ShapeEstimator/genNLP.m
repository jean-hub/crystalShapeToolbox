function obj = genNLP(obj,n_prt) 
%
% Generate NLP matrices and callbacks from data and partition 
% polytope
% ========================================================
%

% Coefs for relative weighting of views
a_12 = obj.a12;
a_13 = obj.a13;

obj.estNLP.pVw = obj.pVw;

obj.estNLP.nprt = n_prt;

n12 = obj.e_nD12;
n13 = obj.e_nD13;
nV = obj.crysModel.part{n_prt}.nV;
m12 = size(obj.Ad12,1);
m13 = size(obj.Ad13,1);
m = size(obj.crysModel.A,1);

obj.estNLP.n12 = n12;
obj.estNLP.n13 = n13;
obj.estNLP.nV = nV;
obj.estNLP.m12 = m12;
obj.estNLP.m13 = m13;
obj.estNLP.pDim = obj.pDim;

d12 = obj.e_d12;
d13 = obj.e_d13;

Ad12 = obj.Ad12;
Ad13 = obj.Ad13;
bd12 = obj.bd12;
bd13 = obj.bd13;

%% Set vertex matrices for partition n_prt
obj.estNLP.V = obj.crysModel.part{n_prt}.V;
obj.estNLP.v = obj.crysModel.part{n_prt}.v;

%% Linear constraints
% Constraint on scaling parameter t
obj.estNLP.bigH = obj.crysModel.part{n_prt}.A;
obj.estNLP.bigH = [obj.estNLP.bigH,zeros(size(obj.estNLP.bigH,1),4*obj.estNLP.nV+3*(n12+n13))];
obj.estNLP.bigh = obj.crysModel.part{n_prt}.b;
% Containment in data polytopes
obj.estNLP.bigH = [obj.estNLP.bigH;...
                   zeros(nV*m12,obj.pDim),kron(eye(nV),Ad12),zeros(nV*m12,2*nV+3*(n12+n13));...
                   zeros(nV*m13,obj.pDim+2*nV),kron(eye(nV),Ad13),zeros(nV*m13,3*(n12+n13))];
obj.estNLP.bigh = [obj.estNLP.bigh;...
                   repmat(bd12,nV,1);...
                   repmat(bd13,nV,1)];
% Containment in model
obj.estNLP.bigH = [obj.estNLP.bigH;...
                   repmat(-obj.crysModel.B,n12,1),zeros(m*n12,4*nV),kron(eye(n12),obj.crysModel.A),zeros(m*n12,3*n13);...
                   repmat(-obj.crysModel.B,n13,1),zeros(m*n13,4*nV+3*n12),kron(eye(n13),obj.crysModel.A)];
obj.estNLP.bigh = [obj.estNLP.bigh;...
                   zeros(m*(n12+n13),1)];

%% Sparsity structure of H
obj.estNLP.spars.bigH = obj.crysModel.part{n_prt}.A; 
obj.estNLP.spars.bigH(find(obj.estNLP.spars.bigH)) = 1;
obj.estNLP.spars.bigH = [obj.estNLP.spars.bigH,...
                         zeros(size(obj.estNLP.spars.bigH,1),4*nV+3*(n12+n13))];

Ad12_ = obj.Ad12; Ad12_(find(Ad12_)) = 1;
Ad13_ = obj.Ad13; Ad13_(find(Ad13_)) = 1;
obj.estNLP.spars.bigH = [obj.estNLP.spars.bigH;...
                         zeros(nV*m12,obj.pDim),kron(eye(nV),Ad12_),zeros(nV*m12,2*nV+3*(n12+n13));...
                         zeros(nV*m13,obj.pDim+2*nV),kron(eye(nV),Ad13_),zeros(nV*m13,3*(n12+n13))];

A_ = obj.crysModel.A; A_(find(A_)) = 1;
B_ = obj.crysModel.B; B_(find(B_)) = 1;
obj.estNLP.spars.bigH = [obj.estNLP.spars.bigH;...
                         repmat(B_,n12,1),zeros(m*n12,4*nV),kron(eye(n12),A_),zeros(m*n12,3*n13);...
                         repmat(B_,n13,1),zeros(m*n13,4*nV+3*n12),kron(eye(n13),A_)];


%% Objective (Re-projection error)
obj.estNLP.bigb = [zeros(4*nV,1);...
                   d12(:)/sqrt(n12);...
                   d13(:)/sqrt(n13)];
               
obj.estNLP.d12 = d12;
obj.estNLP.d13 = d13;              

% Dimension of linear part of optimiser
obj.estNLP.dimX = obj.pDim+4*nV+3*(n12+n13);

%% Set symbolic expressions for hessians

setupOpt(obj.estNLP,a_12,a_13);

obj.estNLP.hssM2Dat12 = str2func(['hessM2Dat_12_',num2str(obj.pDim)]);
obj.estNLP.hssD2Mod12 = str2func(['hessD2Mod_12_',num2str(obj.pDim)]);
obj.estNLP.hssM2Dat13 = str2func(['hessM2Dat_13_',num2str(obj.pDim)]);
obj.estNLP.hssD2Mod13 = str2func(['hessD2Mod_13_',num2str(obj.pDim)]);

%% Set wms points on (q',X')'
obj.estNLP.wms = zeros(4+obj.estNLP.dimX,1);
obj.estNLP.wms(1) = 1;

end



function setupOpt(nlp,a_12,a_13)
%
% Symbolic expressions (objective, gradient, hessian)
% ===================================================
%


nV = nlp.nV;
nD12 = nlp.n12;
nD13 = nlp.n13;
m = nlp.pDim;

fname_M2D_12 = ['hessM2Dat_12_',num2str(m)];
fname_M2D_13 = ['hessM2Dat_13_',num2str(m)];
fname_D2M_12 = ['hessD2Mod_12_',num2str(m)];
fname_D2M_13 = ['hessD2Mod_13_',num2str(m)];


if (exist(fname_M2D_12)~=2)||(exist(fname_M2D_13)~=2)...
        ||(exist(fname_D2M_12)~=2)||(exist(fname_D2M_13)~=2)
    
    % Compute a function to give value and gradient of
    %  f(q,th,x) = |Rx*V*t - x|
    %
    t = sym('t',[m 1]);
    t = sym(t, 'real'); % Shape param
    
    q  = sym('q',[4 1]);
    q  = sym(q, 'real'); % Quaternion
    R = getRot(q);
    
    V = sym('V',[3 m]);
    V = sym(V, 'real'); % Vertex matrix
    
    % ||Rx*V*t - xx||
    
    xx_12 = sym('xx_12',[2 1]);
    xx_12 = sym(xx_12,'real'); % Data point on projection 12
    
    err_12 = (R([1 2],:)*V*t-xx_12)/sqrt(nV);
    cost_12 = err_12'*err_12;
    cost_12 = a_12*cost_12;
    
    fprintf('==> Generate hessian, M2D_12 <==\n');
    
    hCost_12 = hessian(cost_12,[q;t;xx_12]);
    
    % Generate Matlab function
    matlabFunction(hCost_12,...
        'outputs', {'hCost_12'},...
        'vars', {q,t,xx_12,V},...
        'file', sprintf('hessM2Dat_12_%i',m));
    
    xx_13 = sym('xx_13',[2 1]);
    xx_13 = sym(xx_13,'real'); % Data point on projection 13
    
    err_13 = (nlp.pVw*R*V*t-xx_13)/sqrt(nV);
    cost_13 = err_13'*err_13;
    cost_13 = a_13*cost_13;
    
    fprintf('==> Generate hessian, M2D_13 <==\n');
    
    hCost_13 = hessian(cost_13,[q;t;xx_13]);
    
    % Gen. Matlab function
    matlabFunction(hCost_13,...
        'outputs',{'hCost_13'},...
        'vars',{q,t,xx_13,V},...
        'file',sprintf('hessM2Dat_13_%i',m));
    
    % ||Rx*yy-d||
    yy_12 = sym('yy_12', [3 1]);
    yy_12 = sym(yy_12, 'real'); % Model point on projection 12
    d_12 = sym('d_12', [2 1]);
    d_12 = sym(d_12, 'real'); % Data point on projection 12
    
    err_12 = (R([1 2],:)*yy_12-d_12)/sqrt(nD12);
    cost_12 = err_12'*err_12;
    cost_12 = a_12*cost_12;
    
    
    fprintf('==> Generate hessian, D2M_12 <==\n');
    
    hCost_12 = hessian(cost_12,[q;yy_12]);
    
    % Gen. Matlab function
    matlabFunction(hCost_12,...
        'outputs',{'hCost_12'},...
        'vars',{q,yy_12,d_12},...
        'file',sprintf('hessD2Mod_12_%i',m));
    
    
    yy_13 = sym('yy_13',[3 1]);
    yy_13 = sym(yy_13,'real'); % Model point on projection 13
    d_13 = sym('d_13',[2 1]);
    d_13 = sym(d_13,'real'); % Data point on projection 13
    
    err_13 = (nlp.pVw*R*yy_13-d_13)/sqrt(nD13);
    cost_13 = err_13'*err_13;
    cost_13 = a_13*cost_13;
    
    
    fprintf('==> Generate hessian, D2M_13 <==\n');
    
    hCost_13 = hessian(cost_13,[q;yy_13]);
    
    % Gen. Matlab function
    matlabFunction(hCost_13,...
        'outputs',{'hCost_13'},...
        'vars',{q,yy_13,d_13},...
        'file',sprintf('hessD2Mod_13_%i',m));
    
end

end


