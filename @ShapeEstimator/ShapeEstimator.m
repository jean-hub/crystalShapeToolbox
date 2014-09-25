classdef ShapeEstimator < hgsetget
%
% ShapeEstimator class for polytopic shape estimation from projections
% =====================================================================
%
%   INPUTS
%   ------
%
%   - name, crystal name
%   - model, crystal parametric model {A,B}
%   - d12, data points from view 1 (array of size 2*n12)
%   - d13, data points from view 2 (array of size 2*n13)
%
%  AUTHOR(s)
%  ---------
%    
%   (c) Jean-Hubert Hours: EPF Lausanne
%   mailto:jean-hubert.hours@epfl.ch 
%  
%
%  LICENSE
%  -------
%    
%  This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License as published by the Free Software
%  Foundation; either version 2.1 of the License, or (at your option) any later
%  version.
%  This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%  You should have received a copy of the GNU General Public License along with
%  this library; if not, write to the  Free Software Foundation, Inc.,  59 Temple
%  Place, Suite 330,  Boston, MA 02111-1307 USA
%    

properties
    
    
    %
    %% Model
    %
    crysName  % Crystal name (existing: ace, alp, asc, asp, bet, cub, ibu)
    
    crysModel  % Crystal model (contains partition,...)
    
    pDim % Dim. of shape prm 
    
    prtComputed % 1 if partition has been computed, 0 otherwise
    
    
    %
    %% Data
    %
    
    vwAng % Viewing angle for view 2
    pVw % Proj. matrix for view 2
    
    d12  % Data points on view 1 (in R^{2*nD12})
    d13  % Data points on view 2 (in R^{2*nD13})
    
    e_d12  % 'Most relevant' points extracted from d12
    e_d13  % 'Most relevant' points extracted from d13
    
    Ad12  % A matrix of convex hull of e_d12
    bd12  % b matrix of convex hull of e_d12
    Ad13  % A matrix of convex hull of e_d13
    bd13  % b matrix of convex hull of e_d13
    
    nD12  % # data points in view 12
    e_nD12  % # 'most relevant' points in view 12
    nD13 % # data points in view 13
    e_nD13  % # 'most relevant' points in view 13
    
    a12 % Rel. coef. on view 12
    a13 % Rel. coef. on view 13
    
    
    %
    %% Estimator
    %
    nWms  % nWms*nWms*nWms grid of wms points in the estimation loop   
    
    tEst  % Current estimate of the shape parameter
    
    qEst  % Current estimate of the orientation (quaternion)
    
    repErr % Current re-projection error (objecitve to be optimised over)
    
    estNLP  % Local NLP for minimizing re-projection error
    
    optsIPOPT % IPOPT options
    
    tolRepErr % Tolerance on re-projection error (below tolerance => good fit)
    
end


methods
    
    %
    %% Constructor
    %
    function obj = ShapeEstimator(name,model,d12,d13,vwAng,varargin)
        % SHAPEESTIMATOR Creates a shape estimator object
        %
        %
        %
        %
        
        ip = inputParser;
        ip.addRequired('name');
        ip.addRequired('model');
        ip.addRequired('d12');
        ip.addRequired('d13');
        ip.addRequired('vwAng');
        ip.parse(name,model,d12,d13,vwAng,varargin{:});
        ipRes = ip.Results;
        
        obj.vwAng = ipRes.vwAng;
        obj.pVw = [cos(ipRes.vwAng),0,-sin(ipRes.vwAng);...
                   0,1,0];
        
        obj.crysName = ipRes.name;
        obj.crysModel = ipRes.model;
        obj.pDim = size(model.B,2);
        obj.prtComputed = 0;
        
        obj.d12 = ipRes.d12;
        obj.d13 = ipRes.d13;
        obj.nD12 = size(obj.d12,2);
        obj.nD13 = size(obj.d13,2);
        
        % Default value on a12 & a13
        obj.a12 = 1.;
        obj.a13 = 1.;
        
        obj.nWms = 0;
        
        obj.tEst = zeros(obj.pDim,1);
        obj.qEst = zeros(4,1);
    end
    
    %
    %% Set methods
    %
    function obj = set.d12(obj,d12)
        obj.d12 = d12;
    end
    
    function obj = set.d13(obj,d13)
        obj.d13 = d13;
    end
    
    function obj = set.vwAng(obj,vwAng) 
        obj.vwAng = vwAng;
    end
    
    function obj = set.a12(obj,a12) 
        obj.a12 = a12;
    end
    
    function obj = set.a13(obj,a13)
        obj.a13 = a13;
    end
    
    function obj = set.crysName(obj,name)
        obj.crysName = name;
    end
    
    function obj = set.crysModel(obj,model)
        obj.crysModel = model;
    end
    
    function obj = set.nWms(obj,nWms)
        obj.nWms = nWms;
    end
    
    function obj = set.optsIPOPT(obj,optsIPOPT)
        obj.optsIPOPT = optsIPOPT;
    end
    
    function obj = set.tolRepErr(obj,tolRepErr)
        obj.tolRepErr = tolRepErr;
    end
    
    %
    %% Set data points
    %
    obj = setData(obj,d12,d13);
    
    %
    %% Pre-processing of data points
    %
    obj = processData(obj);
    
    %
    %% Compute partition
    %
    obj = computePartition(obj);
    
    %
    %% Generate NLP (from partition polytope, data points and warm-start)
    %
    obj = genNLP(obj,n_prt);
    
    %
    %% 'Global optimization' loop
    %
    obj = estimShape(obj);
    
    %
    %% Solve local NLP using ipopt
    %
    obj = solveNLP_IP(obj);
    
end
    
 
end