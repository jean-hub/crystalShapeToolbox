function crysModel = genCrysModels(name,A,B,pSpace)
%
% Generate crystal model from {A,B} matrices & pSpace polytope
%==============================================================
% Input:
%   - {A,B}, matrices of parametric polytope P(A,Bt) = {x|Ax<=Bt}
%   - pSpace, polytope representing possible values for shape parameter t,
%             should be a box

crysModel.name = name;

crysModel.A = A;
crysModel.B = B;

pSpace.minHRep();

crysModel.H = pSpace.A;
crysModel.h = pSpace.b;

end