function Pprj = getProjPlane(Ply,pMat)
    %
    % Compute polytope projection from matrix pMat (2*3 matrix)
    %
       
    % Project vertices
    pVrts=pMat*(Ply.V)';
    % Minimal representation of the cv hull of the projected vertices
    Pprj=Polyhedron((pVrts)');
    Pprj=Pprj.minVRep;
    
    







