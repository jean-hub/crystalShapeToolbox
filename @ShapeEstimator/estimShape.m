function obj = estimShape(obj)
%
% Global loop for shape estimation
%===================================
% Compute tEst & qEst minimizing the re-projection error 'globally'
%

if ~obj.prtComputed
    error('ShapeEstimator: Error, partition has not been computed.');
end

optSol.repErr = inf; 
for n = 1:length(obj.crysModel.part)        
 
    fprintf('Partition %i\n',n);
    
	%% Generate NLP matrices & callbacks from data
    obj.genNLP(n);
  
    %% Generate warm-starts & run local optimiser
    if obj.estNLP.nV>0
        
        wms = genWms(n,obj);
    
        %% Loop over warm-starts
        for k=1:size(wms,2)
        
            obj.estNLP.wms = wms(:,k);
            
            obj.solveNLP_IP();
            
            if obj.repErr < optSol.repErr
                optSol.repErr = obj.repErr;
                optSol.t = obj.tEst;
                optSol.q = obj.qEst;
            end
            
            fprintf('Best re-projection error so far:%1.10f\n',optSol.repErr);
            disp('Best tEst so far =');disp(optSol.t);
            
            if obj.repErr < obj.tolRepErr
                 optSol.repErr = obj.repErr;
                 optSol.t = obj.tEst;
                 optSol.q = obj.qEst;
                fprintf('STOP: current re-projection error below tolerance.\n'); 
                return;
            end
                
        end
    end
end

obj.repErr = optSol.repErr;
obj.tEst = optSol.t;
obj.qEst = optSol.q;

end


function wms = genWms(n_prt,obj)
%
% Generate wms on q and X
%=========================
%

ph = linspace(0,2*pi-1e-5,obj.nWms);
ps = linspace(0,pi,obj.nWms);
th = linspace(0,pi,obj.nWms);

wms = zeros(4+obj.estNLP.dimX,0);

for i=1:obj.nWms
    for j=1:obj.nWms
        for k=1:obj.nWms
            
            % Wms on scaling prm (Chebycenter of partition polytope)
            t = obj.crysModel.part{n_prt}.P.chebyCenter.x;
            
            % Wms on q
            n = [sin(ps(j))*cos(ph(i));...
                sin(ps(j))*sin(ph(i));...
                cos(ps(j))];
            q = [cos(th(k)/2);...
                n*sin(th(k)/2)];
            q = q/norm(q,2);
            
            R = getRot(q);
            
            % Fictive projected model vertices
            Y = zeros(4*obj.estNLP.nV,1);
            
            % Fictive model points
            Z = zeros(3*(obj.estNLP.n12+obj.estNLP.n13),1);
            
            wms(:,end+1)=[q;...
                          t;...
                          Y;...
                          Z];
        end
    end
end
    
end




