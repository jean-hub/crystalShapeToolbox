function [part] = computeParamVert(poly)


    %%%%%%%%%%%%%%%%%%
    % Compute parameterized vertices pairs [domain,vertex matrix]
    %%%%%%%%%%%%%%%%%
    
    %% Define polytope in combined space
    pc.A=[poly.A,-1.*poly.B];
    pc.b=zeros(size(poly.A,1),1);
    pc.dim=size(pc.A,2);
    pc.R=poly.R; % Rays of polytope in combined space
    pc.Inc=poly.Inc; % Incidence matrix
    
    %% Build binary tree (rays are leaves)
    fprintf('Start building tree...\n');
    bst=buildTree(pc.Inc);
    fprintf('...built\n');  
    
    %% Initialize face lattice
    
    Latt{1}=pc.Inc;
    
    disp('Rays:');
    disp(pc.Inc);
    
    %% Main loop to compute face lattice
    k=1;
    ff=zeros(1,poly.prm.dim);
    ff(1)=size(pc.Inc,2);
    
    while (ff(k)>=2)&&(k<=poly.prm.dim)
    
        ff(k+1)=0;
        Latt{k+1}=[];
        
        % Generate pairs
        prs=genCombinations(size(Latt{k},2),2);
        
        for i=1:size(prs.mat,1)
        
            p1=prs.mat(i,1);
            p2=prs.mat(i,2);
            
            f=Latt{k}(:,p1);
            fp=Latt{k}(:,p2);
            fpp=[];
            
            for i=1:length(f)
                fpp=[fpp,...
                     num2str(str2num(f(i))&str2num(fp(i)))];
            end
            
            [inTree,loc]=searchTree(fpp,bst);
            
            if inTree==0
            
%                 disp('Not in tree');
%                 disp('fpp=');disp(fpp);
%                 disp('fp=');disp(fp);
                
                if strcmp(fpp',f)||strcmp(fpp',fp)
                    fprintf('In delete\n');
                    aux_pp=strmatch(fpp',Latt{k});
                    Latt{k}(:,aux)=[];
                    ff(k)=ff(k)-1;
                end
                
                % Add fpp to bst
                bst{loc}.occ=true;
                Latt{k+1}=[Latt{k+1},fpp'];
                ff(k+1)=ff(k+1)+1;
            end
                
        end
        
        fprintf('%i-faces:\n',k);
        disp(Latt{k});
        fprintf('Number of %i-faces:\n',k);
        disp(size(Latt{k},2));
        
        disp('Number of faces:');
        disp(ff);
        pause;
        
        k=k+1;
    end
    
    
    
    
  
   
    
    
    