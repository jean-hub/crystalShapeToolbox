function obj = solveNLP_IP(obj)
%
% Solve local NLP using IPOPT
%=============================
%

pDim = obj.estNLP.pDim;
nV = obj.estNLP.nV;
nD12 = obj.estNLP.n12;
nD13 = obj.estNLP.n13;
dimX = obj.estNLP.dimX;

%% Lower bounds on optimiser
opsIP.lb = [-1*ones(4,1);...
            -1e10*ones(dimX,1)];
opsIP.lb = (opsIP.lb)';

%% Upper bounds on optimizer
opsIP.ub = [ones(4,1);...
            1e10*ones(dimX,1)];
opsIP.ub = (opsIP.ub)';

%% Lower bound on constraints
opsIP.cl = [1;...
            -1e10*ones(size(obj.estNLP.bigH,1),1)];
opsIP.cl = (opsIP.cl)';

%% Upper bound on constraints
opsIP.cu = [1;...
            obj.estNLP.bigh];
opsIP.cu = (opsIP.cu)';

%% Wms on multiplier
opsIP.lambda = zeros(1,1+size(obj.estNLP.bigH,1));

%% IPOPT options
opsIP.ipopt = obj.optsIPOPT;

%% Callbacks
funcIP.objective = @(X)evalObj(X,obj.estNLP);
funcIP.gradient = @(X)evalGrdObj(X,obj.estNLP);
funcIP.constraints = @(X)evalCon(X,obj.estNLP);
funcIP.jacobian = @(X)evalJacCon(X,obj.estNLP);
funcIP.hessian = @(X,sig,lmb)hLagFun(X,sig,lmb,obj.estNLP);

%funcIP.iterfunc=@callback;

%% Set jacobian structure
jStruct = zeros(1+size(obj.estNLP.bigH,1),4+dimX);
jStruct(1,1:4) = ones(1,4);
jStruct(2:end,5:end) = obj.estNLP.spars.bigH;

funcIP.jacobianstructure = @()sparse(jStruct);

%% Set hessAL structure
hStruct = zeros(4+dimX);
% Quaternion & scaling
hStruct(1:pDim+4,1:pDim+4) = ones(pDim+4);
% yq, yt & yy
hStruct(1:pDim+4,pDim+5:pDim+4+4*nV) = ones(pDim+4,4*nV);
hStruct(pDim+5:pDim+4+4*nV,1:pDim+4) = ones(4*nV,pDim+4);
hStruct(pDim+5:pDim+4+4*nV,pDim+5:pDim+4+4*nV) = kron(eye(2*nV),eye(2));
% zq & zz
hStruct(4+pDim+4*nV+1:4+pDim+4*nV+3*(nD12+nD13),1:4) = ones(3*(nD12+nD13),4);
hStruct(1:4,4+pDim+4*nV+1:4+pDim+4*nV+3*(nD12+nD13)) = ones(4,3*(nD12+nD13));
hStruct(4+pDim+4*nV+1:4+pDim+4*nV+3*(nD12+nD13),...
        4+pDim+4*nV+1:4+pDim+4*nV+3*(nD12+nD13)) = kron(eye(nD12+nD13),ones(3));

funcIP.hessianstructure = @()sparse(tril(hStruct));

%
%% Solve using ipopt
%

[X,info] = ipopt(obj.estNLP.wms,funcIP,opsIP);

sol.x = X;
sol.q = X(1:4); X(1:4) = [];
sol.t = X(1:pDim); X(1:pDim) = [];
sol.Y12 = reshape(X(1:2*nV),2,nV); X(1:2*nV) = [];
sol.Y13 = reshape(X(1:2*nV),2,nV); X(1:2*nV) = [];
sol.Z12 = reshape(X(1:3*nD12),3,nD12); X(1:3*nD12) = [];
sol.Z13 = reshape(X(1:3*nD13),3,nD13); X(1:3*nD13) = [];

sol.lag = info.lambda;
sol.fval = funcIP.objective(sol.x);
sol.econ = funcIP.constraints(sol.x);
sol.nit = info.iter;


% Set estimated shape parameter & rotation
obj.tEst = sol.t;
obj.qEst = sol.q;
obj.repErr = sol.fval;

end


function v = evalObj(X,nlp)
%
% Get objective (re-projection error)
%====================================
%
Phi = evalPhi(X(1:4),nlp);

X = X(5:end);
v = Phi*X-nlp.bigb;
v = (v)'*v;
end


function g = evalGrdObj(X,nlp)
%
% Gradient of objective 
%======================
%
Phi = evalPhi(X(1:4),nlp);
dPhi = evalDevPhi(X(1:4),nlp);
g = zeros(4+nlp.dimX,1);

X = X(5:end);
dff = Phi*X-nlp.bigb;
g(1:4) = 2*([dPhi{1}*X,dPhi{2}*X,dPhi{3}*X,dPhi{4}*X])'*dff;
g(5:end) = 2*(Phi)'*(Phi*X-nlp.bigb);
end


function con = evalCon(X,nlp)
%
% Constraints
%=============
%
con = zeros(1+size(nlp.bigH,1),1);

con(1) = (X(1:4))'*X(1:4);
con(2:end) = nlp.bigH*X(5:end);
end


function jCon = evalJacCon(X,nlp)
%
% Jacobian of constraints
%========================
%
jCon = zeros(1+size(nlp.bigH,1),4+nlp.dimX);

jCon(1,1:4) = 2*(X(1:4))';
jCon(2:end,5:end) = nlp.bigH;

jCon=sparse(jCon);
end


function hessLag = hLagFun(X,sig,lmb,nlp)
%
% Hessian of Lagrangian ( in lower triangular form)
% ================================================
% Input:
%   - X, optimiser
%   - sig, tuning parameter sig_f (see Ipopt doc.)
%   - lmb, Lagrange multiplier 
%   - nlp, local NLP data
%


pDim = nlp.pDim;
nV = nlp.nV;
nD12 = nlp.n12;
nD13 = nlp.n13;

q = X(1:4);
X(1:4) = [];

t = X(1:pDim);
X(1:pDim) = [];

x12 = X(1:2*nV);
X(1:2*nV) = []; % Closest in data (View 12)
x12 = reshape(x12,2,nV);

x13 = X(1:2*nV);
X(1:2*nV) = []; % Closest in data (View 13)
x13 = reshape(x13,2,nV);

y12 = X(1:3*nD12);
X(1:3*nD12) = []; % Closest in Model (View 12)
y12 = reshape(y12,3,nD12);

y13 = X(1:3*nD13);
X(1:3*nD13) = []; % Closest in Model (View 13)
y13 = reshape(y13,3,nD13);

%% Hessians

hess_T_T =zeros(pDim);
hess_Q_Q = zeros(4);
hess_Q_T = zeros(pDim,4);

% View 1
hess_Z_Q_12 = zeros(4,3,nD12);
hess_Z_Z_12 = zeros(3,3,nD12);
% View 2
hess_Z_Q_13 = zeros(4,3,nD13);
hess_Z_Z_13 = zeros(3,3,nD13);

% View 1
hess_Y_T_12 = zeros(pDim,2,nV);
hess_Y_Q_12 = zeros(4,2,nV);
hess_Y_Y_12 = zeros(2,2,nV);

% View 2
hess_Y_T_13 = zeros(pDim,2,nV);
hess_Y_Q_13 = zeros(4,2,nV);
hess_Y_Y_13 = zeros(2,2,nV);


%% Model to data
for i = 1:nV
    
    % View 12
    hC12 = nlp.hssM2Dat12(q,t,x12(:,i),nlp.v{i});
    
    % View 13
    hC13 = nlp.hssM2Dat13(q,t,x13(:,i),nlp.v{i});
    
    % Hessian wrt [q;t;x]
    
    hess_Q_Q = hess_Q_Q+hC12(1:4,1:4)+hC13(1:4,1:4);
    
    hess_T_T = hess_T_T+hC12(5:4+pDim,5:4+pDim)+hC12(5:4+pDim,5:4+pDim);
    hess_Q_T = hess_Q_T+hC12(5:4+pDim,1:4)+hC13(5:4+pDim,1:4);
    
    % View 12
    hess_Y_T_12(:,:,i) = hC12(5:4+pDim,pDim+5:pDim+6);
    hess_Y_Q_12(:,:,i) = hC12(1:4,pDim+5:pDim+6);
    hess_Y_Y_12(:,:,i) = hC12(pDim+5:pDim+6,pDim+5:pDim+6);
    
    % View 13
    hess_Y_T_13(:,:,i) = hC13(5:4+pDim,pDim+5:pDim+6);
    hess_Y_Q_13(:,:,i) = hC13(1:4,pDim+5:pDim+6);
    hess_Y_Y_13(:,:,i) = hC13(pDim+5:pDim+6,pDim+5:pDim+6);
    
end

%% Distance to Poly

% View 12
for i = 1:nD12
    
    hC12 = nlp.hssD2Mod12(q,y12(:,i),nlp.d12(:,i));
    
    % Hessian wrt [q;y_12]
    hess_Q_Q = hess_Q_Q + hC12(1:4,1:4);
    
    hess_Z_Q_12(:,:,i) = hC12(1:4,5:7);
    hess_Z_Z_12(:,:,i) = hC12(5:7,5:7);
end

% View 13
for i = 1:nD13
    
    hC13 = nlp.hssD2Mod13(q,y13(:,i),nlp.d13(:,i));
    
    % Hessian wrt [q;y]
    hess_Q_Q = hess_Q_Q + hC13(1:4,1:4);
    
    hess_Z_Q_13(:,:,i) = hC13(1:4,5:7);
    hess_Z_Z_13(:,:,i) = hC13(5:7,5:7);
end

%% Hessian wrt [t;q;x12;x13;y12;y13]
hessLag = zeros(4+pDim+4*nV+3*(nD12+nD13));

hessLag(5:4+pDim,5:4+pDim) = hess_T_T;
hessLag(1:4,1:4) = hess_Q_Q;

hessLag(5:4+pDim,1:4) = hess_Q_T;
hessLag(1:4,5:4+pDim) = (hess_Q_T)';

% Distance to data

ind = pDim+4+1;

% View 12
for i = 1:nV
    
    hessLag(5:4+pDim,ind:ind+1) = hess_Y_T_12(:,:,i);
    hessLag(1:4,ind:ind+1) = hess_Y_Q_12(:,:,i);
    
    % Transpose
    hessLag(ind:ind+1,5:4+pDim) = (hess_Y_T_12(:,:,i))';
    hessLag(ind:ind+1,1:4) = (hess_Y_Q_12(:,:,i))';
    
    % Diagonal
    hessLag(ind:ind+1,ind:ind+1) = hess_Y_Y_12(:,:,i);
    
    ind = ind+2;
end


% View 13
for i = 1:nV
    
    hessLag(5:4+pDim,ind:ind+1) = hess_Y_T_13(:,:,i);
    hessLag(1:4,ind:ind+1) = hess_Y_Q_13(:,:,i);
    
    % Transpose
    hessLag(ind:ind+1,5:4+pDim) = (hess_Y_T_13(:,:,i))';
    hessLag(ind:ind+1,1:4) = (hess_Y_Q_13(:,:,i))';
    
    % Diagonal
    hessLag(ind:ind+1,ind:ind+1) = hess_Y_Y_13(:,:,i);
    
    ind = ind+2;
end

% Distance to poly

% View 12
for i = 1:nD12
    
    hessLag(1:4,ind:ind+2) = hess_Z_Q_12(:,:,i);
    
    % Transpose
    hessLag(ind:ind+2,1:4) = (hess_Z_Q_12(:,:,i))';
    
    % Diagonal
    hessLag(ind:ind+2,ind:ind+2) = hess_Z_Z_12(:,:,i);
    
    ind = ind+3;
end

% View 13
for i = 1:nD13
    
    hessLag(1:4,ind:ind+2) = hess_Z_Q_13(:,:,i);
    
    % Transpose
    hessLag(ind:ind+2,1:4) = (hess_Z_Q_13(:,:,i))';
    
    % Diagonal
    hessLag(ind:ind+2,ind:ind+2) = hess_Z_Z_13(:,:,i);
    
    ind = ind+3;
end

% Apply sigma coef & add norm-1 constraint
aux = zeros(size(hessLag)); aux(1:4,1:4) = 2*lmb(1)*eye(4);
hessLag = sig*hessLag + aux;

%% Convert to sparse lower triangular form
hessLag = tril(hessLag);
hessLag = sparse(hessLag);

end


function Phi = evalPhi(q,nlp)
%
% Get Phi(q) from data, quaternion and partition index
%========================================================
%


nV = nlp.nV;
nD12 = nlp.n12;
nD13 = nlp.n13;
pDim = nlp.pDim;

Phi = zeros(2*(nV+nV+nD12+nD13),nlp.dimX);

% Rotation
Rq = getRot(q);
Rq12 = Rq([1,2],:);
RR12_t = kron(eye(nV),Rq12)/sqrt(nV);
RR12_z = kron(eye(nD12),Rq12)/sqrt(nD12);
Rq13 = nlp.pVw*Rq;
RR13_t = kron(eye(nV),Rq13)/sqrt(nV);
RR13_z = kron(eye(nD13),Rq13)/sqrt(nD13);

V = nlp.V;

% Scaling prm columns
Phi(1:2*nV,1:pDim) = RR12_t*V;
Phi(2*nV+1:4*nV,1:pDim) = RR13_t*V;

% Fictive 2D projection points
Phi(1:4*nV,pDim+1:pDim+4*nV) = -1*eye(4*nV)/sqrt(nV);

% Data points
Phi(4*nV+1:4*nV+2*nD12,pDim+4*nV+1:pDim+4*nV+3*nD12) = RR12_z;
Phi(4*nV+2*nD12+1:4*nV+2*nD12+2*nD13,pDim+4*nV+3*nD12+1:nlp.dimX) = RR13_z;

end


function dPhi = evalDevPhi(q,nlp)
%
% Evaluate Frechet diff of Phi given q
%=====================================
%

[dR{1},dR{2},dR{3},dR{4}] = getDevRot(q);

nV = nlp.nV;
nD12 = nlp.n12;
nD13 = nlp.n13;
pDim = nlp.pDim;


for i=1:4
    
    dPhi{i} = zeros(2*(nV+nV+nD12+nD13),nlp.dimX);
    
    % Rotation
    dR12 = dR{i}([1,2],:);
    dR12_t = kron(eye(nV),dR12)/sqrt(nV);
    dR12_z = kron(eye(nD12),dR12)/sqrt(nD12);
    dR13 = nlp.pVw*dR{i};
    dR13_t = kron(eye(nV),dR13)/sqrt(nV);
    dR13_z = kron(eye(nD13),dR13)/sqrt(nD13);
    
    V = nlp.V;
    
    % Scaling prm columns
    dPhi{i}(1:2*nV,1:pDim) = dR12_t*V;
    dPhi{i}(2*nV+1:4*nV,1:pDim) = dR13_t*V;
    
    % Data points   
    dPhi{i}(4*nV+1:4*nV+2*nD12,pDim+4*nV+1:pDim+4*nV+3*nD12) = dR12_z;
    dPhi{i}(4*nV+2*nD12+1:4*nV+2*nD12+2*nD13,pDim+4*nV+3*nD12+1:nlp.dimX) = dR13_z;
    
end

end


function [d1R,d2R,d3R,d4R] = getDevRot(q)
%
% Partial devs of R
%

a=q(1);b=q(2);c=q(3);d=q(4);

% Dev wrt q_1
d1R=2*[a , -1*d , c;...
    d , a , -1*b;...
    -1*c , b , a];
% Dev wrt q_2
d2R=2*[b , c , d ;...
    c , -1*b , -1*a;...
    d , a , -1*b];
% Dev wrt q_3
d3R=2*[-1*c , b , a;...
    b , c , d;...
    -1*a , d , -1*c];
% Dev wrt q_4
d4R=2*[-1*d , -1*a , b;...
    a , -1*d , c;...
    b , c , d];
end