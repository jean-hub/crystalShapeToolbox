function [Latt] = cmpLatt(poly)

    %%%%%%%%%%%%%%%%%%
    % Function to compute face lattice 
    %   of any polyhedron
    %%%%%%%%%%%%%%%%%%%
    
    m=size(poly.A,1);
    
    %% Compute faces using backtrack search
       
    for k=1:m
        Latt{k}.num=0;
    end   
    
    aux=zeros(1,m);
    
    for k=1:m
        
        
        dk=zeros(1,m);
        dk(k)=1;
        
        %  Compute set of faces 
        bigF{k}=faceEnum(dk,aux,poly);
        
%         for i=1:size(bigF{k},1) 
%             if sum(bigF{k}(i,:))==m
%                 bigF{k}(i,:)=[];
%                 break;
%             end
%         end
        
        % Order them wrt their dimensions
        for i=1:size(bigF{k},1)
            
            Af=poly.A(find(bigF{k}(i,:))',:);
            
            if rank(Af)==poly.dim
                % Vertex
                
                nm=Latt{1}.num+1;
                
                Latt{1}.f(nm,:)=bigF{k}(i,:);
                
                Latt{1}.num=nm;
                
            elseif 1<=rank(Af)<=(poly.dim-1)
            
                % Edge,...
                
                r=rank(Af);
                nm=Latt{poly.dim-r+1}.num+1;
               
                Latt{poly.dim-r+1}.f(nm,:)=bigF{k}(i,:);
                
                Latt{poly.dim-r+1}.num=nm;
                
            else 
                error('Face has rank m');
            end
            
        end
        
        aux(k)=1;
    end
    
    
%     fprintf('Face lattice:\n');
%     fprintf('==============\n');
%     
%     for k=1:poly.dim
%         
%        fprintf('%i-faces:\n',k-1); 
%         
%        for i=1:size(Latt{k}.f,1)
%             disp(Latt{k}.f(i,:));
%        end
%     end
    
    
    
    
    
    
    
    
    
    
    
    
    