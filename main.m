clear;clc;

%% Load model
load crysModels/ace.mat

%% Generate data (from imaging set-up) specifying viewing angle
dat.npF = 10; 
dat.nMag = 0.00;
vWa = pi/2;pVw = genProjMat(vWa);

dat.q = randn(4,1);dat.q = dat.q/norm(dat.q,2);
P = Polyhedron(ace.part{1}.A,ace.part{1}.b);
dat.t = P.chebyCenter.x;

dat.P = Polyhedron(ace.A*(getRot(dat.q))',ace.B*dat.t);
dat.P12 = projection(dat.P,[1,2]);dat.P12.normalize();dat.P12.minHRep();
%dat.P13 = projection(dat.P,[1,3]);dat.P13.normalize();dat.P13.minHRep();
dat.P13 = getProjPlane(dat.P,pVw);dat.P13.normalize();dat.P13.minHRep();

% Sample projections boundaries
d12 = smpBndPolyBis(dat.P12,dat.npF,dat.nMag);
d13 = smpBndPolyBis(dat.P13,dat.npF,dat.nMag);


%% Run ShapeEstimator

% Instantiate ShapeEstimator
est = ShapeEstimator('aceta',ace,d12,d13,vWa);

% Compute partition from crystal model
est.computePartition();

% Extract relevant data points
est.processData();

% Set IPOPT options
optsIPOPT.mu_init = 1e-2; % Init. barrier parameter
optsIPOPT.mu_strategy = 'adaptive'; % Barrier update strategy
optsIPOPT.tol = 1e-7; % Tolerance on criticality
optsIPOPT.max_iter = 200; % Max. # iterations
optsIPOPT.hessian_approximation = 'exact';  %'exact' or 'limited-memory';
optsIPOPT.print_level = 0; 
optsIPOPT.print_user_options = 'no';
est.optsIPOPT = optsIPOPT;

% Set tolerance on re-projection error
est.tolRepErr = 1e-2;

% Set size of wms grid
est.nWms = 3;

% Run estimation loop
est.estimShape();


%% Compare estimated and true polytopes

Rest = getRot(est.qEst);
Pest = Polyhedron(ace.A*(Rest)',ace.B*est.tEst);

figure(9);clf;
plot(dat.P,'Color','b');hold on;plot(Pest,'Color','r');
