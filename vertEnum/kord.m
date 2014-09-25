function [bool] = kord(key1,key2)

    %%%%%%%%%%%%%
    % Returns 1 if key1<key2, 0 if key1>key2
    % 2 keys cannot be equal
    %%%%%%%%%%%%%
    
    lk1=length(key1);
    lk2=length(key2);
    
    cnt=0;
    
    if lk1<lk2
        
        bool=true;
    elseif lk1>lk2
        
        bool=false;
    else
        
        % Check if 
        for i=1:lk1
            
            if (key1(i)==1)&&(key2(i)==0)
                bool=false;
                break;
            elseif (key1(i)==0)&&(key2(i)==1) 
                bool=true;
                break;
            elseif (key1(i)==key2(i))
                cnt=cnt+1;
            end
            
        end
        
        if cnt==lk1
            error('Keys are equal');
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    