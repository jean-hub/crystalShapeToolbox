function [bool,loc] = searchTree(ky,bst)

    %%%%%%%%%%%%%%%%
    % Find ky in bst and return index and occ
    %%%%%%%%%%%%%%%
    
    if strcmp(ky(1),'0')
        idx=2;
    elseif strcmp(ky(1),'1')
        idx=3;
    else
        error('Error in searching');
    end
    
    while 1
    
        if strcmp(bst{idx}.key,ky)
            bool=bst{idx}.occ;
            loc=idx;
            break;
        else
            
            if strcmp(ky(length(bst{idx}.key)+1),'0')
                
                idx=2^(floor(log2(idx))+1)+bin2dec(bst{idx}.lft);
                
            elseif strcmp(ky(length(bst{idx}.key)+1),'1')                
                
                idx=2^(floor(log2(idx))+1)+bin2dec(bst{idx}.rght);
                
            else
                error('Error in searching');
            end
      
        end
    end
   
    
    