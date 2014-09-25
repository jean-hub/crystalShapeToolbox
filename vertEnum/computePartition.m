function [part] = computePartition(poly,varargin)

    p=inputParser;
    p.addRequired('poly');
    p.addOptional('tol',1e-4);
    p.parse(poly,varargin{:});
    opts=p.Results;
    
    %% Function to compute a partition of the parameters space

    % Generate possible vertices combinations
    ints=genCombinations(poly.m,poly.d);    
    nCk=ints.count-1;
    vMat=inf(nCk*poly.m,poly.prm.dim);

    % Among intersections, take as vertex if it is a basis
    for i=1:size(ints.mat,1)

        ineqs=ints.mat(i,:);
        
        if  abs(det(poly.A(ineqs',:)))>1e-3

            a=inv(poly.A(ineqs',:))*poly.B(ineqs',:);

            vMat((i-1)*poly.m+1:i*poly.m,:)=[poly.A*a,-poly.B]*[eye(poly.prm.dim);eye(poly.prm.dim)];  

        end        
    end
    
    % Partition set of parameters
    thetaFeas=Polyhedron(poly.H,poly.h);  
    diff=thetaFeas; 

    part=[]; 
    while diff.isFullDim
        
        if ~isempty(part)
            diff=thetaFeas\[part.P];
        else
            diff=thetaFeas;
        end

        
        th0=diff(1).chebyCenter.x;
        
        feas=inf(nCk*poly.m,poly.prm.dim);
        verts=zeros(nCk,1);

        % Test all vertices valid at theta_0
        for i=1:nCk
            mat=vMat((i-1)*size(poly.A,1)+1:i*size(poly.A,1),:);
            
            if isempty(th0)==0
                if (mat*th0<opts.tol) & isempty(find(mat==inf)) % If vertex feasible at theta0
                    feas((i-1)*size(poly.A,1)+1:i*size(poly.A,1),:)=mat;
                    verts(i)=1;
                end 
            end
        end
        indFeas=find(feas<inf); % Index of matrix of feasible thetas 
        l_indFeas=length(indFeas);
        % Feasible region around theta
        feasReg=Polyhedron([reshape(feas(indFeas),l_indFeas/poly.prm.dim,poly.prm.dim);poly.H],...
                [zeros(l_indFeas/poly.prm.dim,1);poly.h]);
           
        H=feasReg.H(:,1:end-1);
        K=feasReg.H(:,end);

        if feasReg.isFullDim
            part(end+1).P=Polyhedron(H,K);
            part(end).verts=find(verts);
            part(end).extreme=length(part(end).P.V);
            part(end).list=ints.mat;
            part(end).inqs=part(end).list(part(end).verts,:);
        else
            break;
        end
        
    end



%     p=inputParser;
%     p.addRequired('poly');
%     p.addOptional('tol',1e-4);
%     p.parse(poly,varargin{:});
%     opts=p.Results;
%     
%     %% Function to compute a partition of the parameters space
% 
%     % Generate possible vertices combinations
%     ints=genCombinations(poly.numInqs,poly.dim);
%     nCk=ints.count-1;
%     vMat=inf(nCk*poly.numInqs,poly.prm.dim);
% 
%     % Among intersections, take as vertex if it is a basis
%     for i=1:size(ints.mat,1)
% 
%         ineqs=ints.mat(i,:);
% 
%         if  abs(det(poly.A(ineqs',:)))>1e-3
% 
%             a=inv(poly.A(ineqs',:))*poly.B(ineqs',:);
% 
%             vMat((i-1)*poly.numInqs+1:i*poly.numInqs,:)=[poly.A*a,-poly.B]*[eye(poly.prm.dim);eye(poly.prm.dim)];  
% 
%         end        
%     end
% 
%     % Partition set of parameters
%     thetaFeas=polytope(poly.H,poly.h);  
%     diff=thetaFeas; 
% 
%     part=[]; 
%     while isfulldim(diff)
% 
%         if ~isempty(part)
%             diff=thetaFeas\[part.P];
%         else
%             diff=thetaFeas;
%         end
% 
%         
%         th0=chebyball(diff(1));
%         
%         feas=inf(nCk*poly.numInqs,poly.prm.dim);
%         verts=zeros(nCk,1);
% 
%         % Test all vertices valid at theta_0
%         for i=1:nCk
%             mat=vMat((i-1)*poly.numInqs+1:i*poly.numInqs,:);
%             if (mat*th0<opts.tol) & isempty(find(mat==inf)) % If vertex feasible at theta0
%                 feas((i-1)*poly.numInqs+1:i*poly.numInqs,:)=mat;
%                 verts(i)=1;
%             end 
%         end
%         indFeas=find(feas<inf); % Index of matrix of feasible thetas 
%         l_indFeas=length(indFeas);
%         % Feasible region around theta
%         feasReg=polytope([reshape(feas(indFeas),l_indFeas/poly.prm.dim,poly.prm.dim);poly.H],...
%                 [zeros(l_indFeas/poly.prm.dim,1);poly.h]);
%            
%         [H,K]=double(feasReg);
% 
%         if isfulldim(feasReg)
%             part(end+1).P=polytope(H,K);
%             part(end).verts=find(verts);
%             part(end).extreme=length(extreme(part(end).P));
%             part(end).list=ints.mat;
%             part(end).inqs=part(end).list(part(end).verts,:);
%         else
%             break;
%         end
%     
%     end

    





