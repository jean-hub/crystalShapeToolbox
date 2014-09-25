clear;clc;

%
%% How to create a crystal model
%

% Matrices computed from Miller indices & dependent facets
A = [0.8084    0.5886         0
     0.8084   -0.5886         0
   -0.8084   -0.5886         0
   -0.8084    0.5886         0
   -0.4324         0    0.9017
    0.4324         0   -0.9017
   -0.8298         0    0.5581
    0.8298         0   -0.5581];
B = [1     0     0
     1     0     0
     1     0     0
     1     0     0
     0     1     0
     0     1     0
     0     0     1
     0     0     1];

% Matrices of parameter space polytope
H = [eye(size(B,2));...
     -1*eye(size(B,2))];
h = [10*ones(size(B,2),1);...
     zeros(size(B,2),1)];
 
% Generate model to be used by ShapeEstimator
testModel = genCrysModels('test',A,B,Polyhedron(H,h));

% Save model
% save('testModel.mat','testModel');
