function [Latt]=cmpFaces(poly)

    %%%%%%%%%%%%%%
    % Compute face-lattice of a polytope
    %%%%%%%%%%%%%%%
    
    %% Build binary tree (rays are leaves)
    fprintf('Start building tree...\n');
    bst=buildTree(poly.Inc);
    fprintf('...built\n');  
    
%     disp('Occupancy:');
%     for k=1:length(bst)
%         disp(bst{k}.occ);
%     end
    
    %% Initialize face lattice
    
    Latt{1}=poly.Inc;
    
    disp('Incidence matrix:');
    disp(poly.Inc);
    
    %% Main loop to compute face lattice
    k=1;
    ff=zeros(1,poly.dim);
    ff(1)=size(poly.Inc,2);
    
    emptF=[];
    for i=1:size(poly.Inc,1)
        emptF=[emptF,'0'];
    end
    
    
    while (ff(k)>1)
    
        fprintf('k=%i\n',k);
        
        ff(k+1)=0;
        Latt{k+1}=[];
        
        % Generate pairs
        prs=genCombinations(size(Latt{k},2),2);
       
        for i=1:size(prs.mat,1)
            
            p1=prs.mat(i,1);
            p2=prs.mat(i,2);
            
            datPrs{i}=[Latt{k}(:,p1),Latt{k}(:,p2)];
        end
        
        fprintf('Number of combs:%i\n',size(prs.mat,1));
        
        for i=1:size(prs.mat,1)
            
            f=datPrs{i}(:,1);
            fp=datPrs{i}(:,2);
            fpp=[];
            
%             disp('f=');disp(f');
%             disp('fp=');disp(fp');
            
            for i=1:length(f)
                fpp=[fpp,...
                     num2str(str2num(f(i))&str2num(fp(i)))];
            end
            
%             disp('fpp=');disp(fpp);
            
%             disp('Comparison:');
%             disp(strcmp(fpp',f)|strcmp(fpp',fp));
            
            [inTree,loc]=searchTree(fpp,bst);
            
            %disp('inTree=');disp(inTree);
            
            if (inTree==0)&(strcmp(fpp,emptF)==0)

                if strcmp(fpp',f)|strcmp(fpp',fp)
                    fprintf('In delete\n');
                    aux_pp=strmatch(fpp',Latt{k});
                    Latt{k}(:,aux_pp)=[];
                    ff(k)=ff(k)-1;
                end

            %    disp('Not in tree');
%               disp('fpp=');disp(fpp);
%               disp('fp=');disp(fp);   
                
                % Add fpp to bst
                bst{loc}.occ=true;
                Latt{k+1}=[Latt{k+1},fpp'];
                ff(k+1)=ff(k+1)+1;
            
            else
                
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
    
    
   
    
    