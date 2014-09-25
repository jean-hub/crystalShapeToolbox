function [smp] = smpBndPolyBis(P,NpInt,noise)

    % P: Polytope
    % NpInt: Number of points per interval

    
    poly.dim = size(P.H,2)-1;
    poly.A = P.H(:,1:poly.dim);
    poly.B = P.H(:,end); 
    
    if poly.dim~=2
        error('Not a polygon, only works for 2-polytopes at the moment');
    end
    
    poly.V = P.V;poly.V = poly.V';
    poly.c = P.chebyCenter.x; % Cheby center
   
    
    % Compute angles
    poly.ang = zeros(1,size(poly.V,2));
    for i = 1:size(poly.V,2)
        
        n=norm(poly.V(:,i)-poly.c);
        cth=(poly.V(1,i)-poly.c(1))/n;
        sth=(poly.V(2,i)-poly.c(2))/n;
        
        if sth>=0
            poly.ang(i)=acos(cth);
        else
            poly.ang(i)=2*pi-acos(cth);
        end
        
    end
    poly.angSorted=sort(poly.ang);
    
    
    % Sample points
    
%     figure(23);clf;
%     plot(P,struct('color','w','edgecolor','b','linewidth',3));
%     hold on;
    
    smp=zeros(2,0);
    for i=1:size(poly.angSorted,2)-1
        
        angL=poly.angSorted(i);
        angU=poly.angSorted(i+1);
        vL=find(poly.ang==angL);
        vU=find(poly.ang==angU);
        
        cur.vec=linspace(0,1,NpInt);
         
        for j=1:size(cur.vec,2)
            % Compute the sample point as the intersection of an edge and a
            % ray from the chebychev center
            
            smp(:,end+1)=poly.V(:,vL)+cur.vec(j)*(poly.V(:,vU)-poly.V(:,vL));
        end
%         plot(smp(1,:),smp(2,:),'r.','MarkerSize',12);
%         hold on;
%         pause;
    end
    
    angL=poly.angSorted(i+1);
    angU=min(poly.angSorted);
    vL=find(poly.ang==angL);
    vU=find(poly.ang==angU);
    cur.vec=linspace(0,1,NpInt);
    
    for j=1:size(cur.vec,2)
        smp(:,end+1)=poly.V(:,vL)+cur.vec(j)*(poly.V(:,vU)-poly.V(:,vL));
    end
%      plot(smp(1,:),smp(2,:),'r.','MarkerSize',12);

    
    % Add noise
    
    smp=smp+noise*randn(size(smp));
    
    
    
    
    
    
    
    