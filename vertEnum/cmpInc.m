function [incM,satM] = cmpInc(poly)

    %%%%%%%%%%%%%%
    % Computes incidence matrix
    %%%%%%%%%%%%%%
    
    nr=size(poly.R,2);
    [m,d]=size(poly.A);
    
    incM=zeros(m,nr);
    satM=zeros(m,nr);
    
    for i=1:m
    
        a=poly.A(i,:);
        
        for j=1:nr
            
            ry=poly.R(:,j);
            
            satM(i,j)=abs(a*ry);
            
            if (abs(a*ry)<1e-9)
                incM(i,j)=1;
            end
        end
    end
   