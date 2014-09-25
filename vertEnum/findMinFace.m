function [F] = findMinFace(R,poly)

    %%%%%%%%%%%%%
    % Find minimal face containing R
    %%%%%%%%%%%%%
    
    
    m=size(poly.A,1);
    F=R;    
    
    for i=1:m
        
        aux=zeros(1,m);
        
        
        if R(i)==0 
            
            aux(i)=1;
            
            res=resFacePoly(R,aux,poly);
            
            if res==0
                F(i)=1;
            end
        end
        
    end