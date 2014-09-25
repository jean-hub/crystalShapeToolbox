classdef btnode < handle
    
    %%%%%%%%%%
    % Node in binary search tree
    %%%%%%%%
    
    properties(SetAccess = private)
        Key
        Occupied
        kLeft
        kRight
    end
    
    
    methods
    
        function node = btnode(Key)
            % Build node from key
            if nargin > 0 
                node.Key=Key;
                node.Occupied=false;
            end
        end
            
        function disp(node)
            % Display Nodes keys
            fprintf('Key:\n');
            disp(node.Key);
        end
        
    end
    
end