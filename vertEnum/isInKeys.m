function [bool] = isInKeys(key,kList)

    %%%%%%%%%%%%
    % Check if a given key is in kList
    %%%%%%%%%%%%
    
    bool=false;
    
    for j=1:size(kList,2)
        
        if strcmp(key,kList(:,j)')
            bool=true;
            break;
        end
    end