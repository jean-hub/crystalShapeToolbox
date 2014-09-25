function [bst] = buildTree(vKeys)

    %%%%%%%%%%%%%%%
    % Function to build a binary tree from a set of ray keys 
    %%%%%%%%%%%%%%%
    
    nI=size(vKeys,1);
    nK=size(vKeys,2);
    
    fprintf('Number of keys:%i\n',nK);
    
    nNds=sum(2.^(0:nI));
    
    for k=1:nNds
    
        % If root
        if k==1
            bst{k}.key=[];
            bst{k}.rght=num2str(1);
            bst{k}.lft=num2str(0);
            bst{k}.occ=true;
        else
        
            ky=num2str(dec2bin(k-2^floor(log2(k))));
            lng=floor(log2(k));
        
            if length(ky)<lng
                str=[];
                for j=1:(lng-length(ky))
                    str=[str,'0'];
                end
                ky=[str,ky];
            end
           
            if isInKeys(ky,vKeys)  
           
                bst{k}.key=ky;
                bst{k}.rght=[];
                bst{k}.lft=[];
                bst{k}.occ=true;
                
            else
                
                bst{k}.key=ky;
                bst{k}.rght=[bst{k}.key,num2str(1)];
                bst{k}.lft=[bst{k}.key,num2str(0)];
                bst{k}.occ=false;
            end
           
        end
     
    end   
    
    occRate=0;
    for k=1:length(bst)
        if bst{k}.occ
            occRate=occRate+1;
        end
    end
    
    fprintf('Occupancy:%i\n',occRate);
    

    
    