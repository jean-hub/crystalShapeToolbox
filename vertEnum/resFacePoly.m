function [bool] = resFacePoly(R,S,poly)

    %%%%%%%%%%%%%%%
    % Restricted face of polyhedron:
    % Says if there is a face of poly that contains R and does not 
    % intersect S.
    %%%%%%%%%%%%%%%
    
    m=size(poly.A,1);
    
    
    AR=poly.A(find(R)',:);
    bR=poly.b(find(R)');
    
    AS=poly.A(find(S)',:);
    bS=poly.b(find(S)');
    
    un=mod(R+S,2);
    CMP=zeros(1,m);
    
    for i=1:m
        if un(i)==0
            CMP(i)=1;
        end
    end
    
    AT=poly.A(find(CMP)',:);
    bT=poly.b(find(CMP)');
    
    
    lp.H=zeros(poly.dim+1);
    
    lp.f=zeros(poly.dim+1,1);
    lp.f(end)=-1;
    
    if (isempty(AS)==0)&&(isempty(bS)==0)
    
        lp.A=[AR,zeros(size(AR,1),1)];
        lp.A=[lp.A;
              AS,ones(size(AS,1),1)];
        lp.A=[lp.A;
              AT,zeros(size(AT,1),1)];
        lp.A=[lp.A;
              zeros(1,poly.dim),1];
      
        lp.b=[bR;
              bS;
              bT;
              1];
    
        % Solve feasibility LP
        [XMIN,FMIN,SOLSTAT,DETAILS]=cplexint(lp.H,lp.f,lp.A,lp.b,1:size(AR,1),[],[],[],[],[],[]);
                      
        
        if XMIN(end)>0
            bool=true;
        else
            bool=false;
        end
    
    else
        
        bool=true;
        
    end
    
    
    
    