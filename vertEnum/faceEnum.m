function [Fcs] = faceEnum(R,S,poly)

    %%%%%%%%%%%%%
    % Recursive function for face enumeration
    %%%%%%%%%%%%%
    
    m=size(poly.A,1);
    
    J=[];
    Fcs=[];
    
    if resFacePoly(R,S,poly)
        % Exists face
        
        % Find minimal face containing R
        F=findMinFace(R,poly);
        Fcs=[Fcs;F];
        
        % Build J
        if isempty(S)==0
        
            for i=1:m
                if (F(i)==0)&&(S(i)==0)
                    J=[J,i];
                end
            end
        
        else
           
            for i=1:m
                if (F(i)==0)
                    J=[J,i];
                end
            end   
        end    
    
        aux=[];
        for j=1:length(J)
            
            auxF=F;
            auxF(J(j))=1;
            
            auxS=S;
            auxS(aux)=1;
            
            Fcs=[Fcs;faceEnum(auxF,auxS,poly)];
            
            aux=[aux,J(j)];
        end
    
    end