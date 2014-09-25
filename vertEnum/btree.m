classdef btree < handle
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Implements binary search tree
    %%%%%%%%%%%%%%%%%%%%
    
    properties(SetAccess=private)
        keys
        endKs
        root
    end
    
    methods
        
        function tree = btree(n) 
            % Create binary tree with keys of max length n
            keys=zeros(n,sum(2.^(1:n)));
            endKs=1;
            root=btnode([-1*inf]);
            root.Occupied=false;
            root.kLeft=[];
            root.kRight=[];
        end
    
        function insert(key,occ,node)
            % Insert key in binary tree from node
            
            if (node.Key==-1*inf)
                
                % Either tree is empty
                if endKs==1
                    
                    
                    
                    
                    
                    
                else 
                    
                    
                    
                end
                
                
            
            else
            
                
                
                
            end
    
        end
        
        function delete
            % Delete tree
            keys=[];
            
        end
    
        function occ=search(key)
            % Function to return if node at key is occupied 
            
            
        end
    
    end
    
   
    
end