function [R] = cmpRays(A)

    %%%%%%%%%%%%%%%%%%%%%
    % Compute rays of polyhedral cone  
    % in form {x | Ax \leq 0}
    %%%%%%%%%%%%%%%%%%%%%%
    
    [m,d]=size(A);
    
    R=zeros(d+m,2*d+m);
    
    %% Compute rays of slack polytope
    
    for j=1:d
        
        ei=zeros(d,1);
        ei(j)=1;
        
        R(:,j)=[ei;...
                A*ei];
    end
    
    for j=d+1:2*d
    
        ei=zeros(d,1);
        ei(j-d)=1;
        
        R(:,j)=-1*[ei;...
                    A*ei];
    
    end
    
    for j=2*d+1:2*d+m
        
        ej=zeros(m,1);
        ej(j-2*d)=1;
        
        R(:,j)=[zeros(d,1);...
                ej];
        
    end
   
    
    %% Compute intersection with 0-subspace
    
    for k=1:m
    
        
        % Intersect with hyperplane H_k
        Raux=[];
        
        z_idx=find(R(d+k,:)==0);
        
        if ~isempty(z_idx)
            Raux=[Raux,R(:,z_idx)];
        end
        
       
        pos_idx=find(R(d+k,:)>0);
        neg_idx=find(R(d+k,:)<0);
       
        fprintf('Number of positive coefs=%i\n',length(pos_idx));
        fprintf('Number of neg coefs=%i\n',length(neg_idx));
        
        for i=1:length(pos_idx)
            for j=1:length(neg_idx)
                
                Raux=[Raux,R(d+k,pos_idx(i))*R(:,pos_idx(i))+...
                           (-1*R(d+k,neg_idx(j)))*R(:,neg_idx(j))];
                
            end
        end    
        
        R=Raux;
    end
    
    disp('Rays after F-M:');
    disp(R);
    
    
    
    
    
    
    
    
    
    
    
    
    