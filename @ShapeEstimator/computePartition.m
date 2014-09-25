function obj = computePartition(obj)
%            
% Compute partition from crystal model
%======================================
%
%

nFacets = size(obj.crysModel.A,1);
dPoly = size(obj.crysModel.A,2);
dPrm = size(obj.crysModel.B,2);

% Generate possible vertices combinations
ints = genCombinations(nFacets,dPoly);
nCk = ints.count-1;
vMat = inf(nCk*nFacets,dPrm);

% Among intersections, take as vertex if it is a basis
for i = 1:size(ints.mat,1)
    
    ineqs = ints.mat(i,:);
    
    if  abs(det(obj.crysModel.A(ineqs',:)))>1e-3
        
        a = inv(obj.crysModel.A(ineqs',:))*obj.crysModel.B(ineqs',:);
        
        vMat((i-1)*nFacets+1:i*nFacets,:) = [obj.crysModel.A*a,-obj.crysModel.B]*[eye(dPrm);eye(dPrm)];
    end
end

% Partition set of parameters
thetaFeas = Polyhedron(obj.crysModel.H,obj.crysModel.h);
diff = thetaFeas;


part = [];
while diff.isFullDim
    
    if ~isempty(part)
        diff = thetaFeas\[part.P];
    else
        diff = thetaFeas;
    end
    
    
    th0 = diff(1).chebyCenter.x;
    
    feas = inf(nCk*nFacets,dPrm);
    verts = zeros(nCk,1);
    
    % Test all vertices valid at theta_0
    for i = 1:nCk
        mat = vMat((i-1)*nFacets+1:i*nFacets,:);
        
        if isempty(th0)==0
            if (sum(mat*th0<1e-4*ones(nFacets,1))==nFacets)&&isempty(find(mat==inf)) % If vertex feasible at theta0
                feas((i-1)*nFacets+1:i*nFacets,:) = mat;
                verts(i)=1;
            end
        end
    end
    indFeas = find(feas<inf); % Index of matrix of feasible thetas
    l_indFeas = length(indFeas);
    % Feasible region around theta
    feasReg = Polyhedron([reshape(feas(indFeas),l_indFeas/dPrm,dPrm);obj.crysModel.H],...
                         [zeros(l_indFeas/dPrm,1);obj.crysModel.h]);
    
    H = feasReg.H(:,1:end-1);
    K = feasReg.H(:,end);
    
    if feasReg.isFullDim
        
        part(end+1).P = Polyhedron(H,K);
        part(end).verts = find(verts);
        part(end).extreme = length(part(end).P.V);
        part(end).list = ints.mat;
        part(end).inqs = part(end).list(part(end).verts,:);
        
        % Compute vertices
        for i = 1:length(part)
            
            part(i).V = [];
            inqs = part(i).inqs;
            
            for j = 1:length(inqs)
                inq = inqs(j,:);
                part(i).v{j} = inv(obj.crysModel.A(inq',:))*obj.crysModel.B(inq',:);
                part(i).V = [part(i).V;...
                             part(i).v{j}];
            end
            
            part(i).nV = length(inqs);
            
        end
    else
        break;
    end
    
end

n = 1;
while 1
    
    if n>length(part)
        break;
    else
        if (part(n).P==thetaFeas) 
            n = n+1;
        else
            obj.crysModel.part{n}.A = part(n).P.A;
            obj.crysModel.part{n}.b = part(n).P.b;
            obj.crysModel.part{n}.P = part(n).P;
            obj.crysModel.part{n}.V = part(n).V;
            obj.crysModel.part{n}.nV = part(n).nV;
            for j = 1:part(n).nV
                obj.crysModel.part{n}.v{j} = part(n).v{j};
            end
            n = n+1;
        end
    end
end

obj.prtComputed = 1;
      
end


function [all] = genCombinations(n,p)

all.mat=zeros(nchoosek(n,p),p);
all.count=1;
if n==p
    all.mat(1,:)=1:n;
else
    list=zeros(1,p);
    index=1;
    all=computeList(index,n,p,list,all);
end

end


function [all] = computeList(index,n,p,list,all)

if index>=p+1
    all.mat(all.count,:)=list;
    all.count=all.count+1;
    return;
end
start=1;
if index>1
    start=list(index-1)+1;
end
for i=start:n
    list(index)=i;
    all=computeList(index+1,n,p,list,all);
end

end