function [same] = checkRays(sR1,sR2)

    %%%%%%%%%%%%%%
    % Function to check whether 2 sets of rays are equal or not
    %%%%%%%%%%%%%%%
    
    n1=size(sR1,1);
    n2=size(sR2,1);
    
    if n1~=n2
        error('Ray dimension not equal');
    end
    
    q1=size(sR1,2);
    q2=size(sR2,2);
    
    
    if q1~=q2
        
        same=false;
    else
        
        for j=1:q1
            
            isIn=false;
            
            for k=1:q2
                if norm(sR1(:,j)-sR2(:,k),2)<1e-4
                    isIn=true;
                    break;
                end
            end
            
            if isIn==false
                same=false;
                break;
            else
                same=true;
            end
            
        end
        
    end
        
    
    
    
    
    
    