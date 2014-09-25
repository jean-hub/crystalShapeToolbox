clear;

dir='/Users/jean/Desktop/Publications/PaperPolytopeRec/Figs/';

% poly.A=[-1*eye(poly.dim);...
%         eye(poly.dim)];
% poly.b=[zeros(poly.dim,1);...
%         ones(poly.dim,1)];

% poly.A=[ones(1,3);...
%         -1*eye(poly.dim)];
% poly.b=[1;...
%         zeros(poly.dim,1)];

% poly.A=-1*eye(3);
% poly.b=zeros(3,1);

chs=['cube','ibuprofen','aceta','aspirin','lglubeta','lglualpha','asc'];

% Crystal choice
crys_choice='aceta';

%% Create parameterized polytope

if strcmp(crys_choice,'random') % Random polytope

    A=randn(10,3);
    
    B=rand(10,3);
    
    % Rays 
    Cn=Polyhedron([A,-1*B],zeros(size(A,1),1));
    Cn.minHRep;
    Cn.normalize;
    R=Cn.R;
    R=R';

elseif strcmp(crys_choice,'cube') % cube
    
    fprintf('======== CUBE ========\n');
    
    A=[eye(3);-eye(3)];
    
    B=[eye(3);eye(3)];
    
elseif strcmp(crys_choice,'ibuprofen') % ibuprofen
    
    fprintf('======= IBUPROFEN =======\n');


    A=[1.0000         0.         0.
       -1.0000         0.         0.
        0.1627         0.    0.9867
       -0.1627         0.   -0.9867
        0.0972    0.8020    0.5894
       -0.0972    0.8020   -0.5894
        0.0972   -0.8020    0.5894
       -0.0972   -0.8020   -0.5894];

    B=[1.     0      0.
       1.     0      0.
       0.     1.      0.
       0.     1.      0.
       0.     0.      1.
       0.     0.      1.
       0.     0.      1.
       0.     0.      1.];
   
   % Rays of cone in combined space
   load ibuRays; 
   
   % Incidence matrix of cone in combined space
   load ibuInc;
   
    
elseif strcmp(crys_choice,'aceta') % aceta
    
    fprintf('======= ACETA ========\n');
    
    A=[0.8084, 0.5886, 0
    0.8084, -0.5886 , 0
   -0.8084, -0.5886, 0
   -0.8084, 0.5886, 0
   -0.4324 , 0 , 0.9017
    0.4324 , 0 , -0.9017
   -0.8298 , 0 , 0.5581
    0.8298 , 0 , -0.5581];
    
    B=[ 1  0  0
        1 0  0
        1  0  0
        1 0  0
        0  1  0
        0  1  0
        0  0   1
        0  0  1];

    % Rays of cone in combined space
   load acetaRays;
   
   % Incidence matrix of cone in combined space
   load aceInc; 


elseif strcmp(crys_choice,'aspirin') % aspirin
    
    fprintf('===== ASPIRIN ========\n');
   
    A=[ 1.0000         0         0
   -1.0000         0         0
   -0.0990         0    0.9951
    0.0990         0   -0.9951
   -0.0857    0.5007    0.8614
    0.0857    0.5007   -0.8614
    0.0857   -0.5007   -0.8614
   -0.0857   -0.5007    0.8614
    0.8663    0.4995         0
   -0.8663    0.4995         0
    0.8663   -0.4995         0
   -0.8663   -0.4995         0];

    B=[1     0     0     0
     1     0     0     0
     0     1     0     0
     0     1     0     0
     0     0     1     0
     0     0     1     0
     0     0     1     0
     0     0     1     0
     0     0     0     1
     0     0     0     1
     0     0     0     1
     0     0     0     1];
 
 
    % Rays of cone in combined space
    load aspiRays;
    
    % Incidence matrix of cone in combined space
    load aspInc;
    
    
elseif strcmp(crys_choice,'lglubeta') % lglubeta 
    
    fprintf('========= BETA =========\n');
    

 A=[0    1.0000         0
         0   -1.0000         0
    0.8023         0    0.5969
    0.8023         0   -0.5969
   -0.8023         0    0.5969
   -0.8023         0   -0.5969
         0    0.6255    0.7803
         0    0.6255   -0.7803
         0   -0.6255    0.7803
         0   -0.6255   -0.7803];
    
    B=[1     0     0
     1     0     0
     0     1     0
     0     1     0
     0     1     0
     0     1     0
     0     0     1
     0     0     1
     0     0     1
     0     0     1];
 
    % Rays of cone in combined space
    load betaRays;

    % Incidence matrix 
    load betaInc;
    
 
    
elseif strcmp(crys_choice,'lglualpha') % lglualpha
    
    fprintf('======== ALPHA ==========\n');
    
    
    A=[0         0   -1.0000
   -0.4630   -0.6755   -0.5738
    0.4630    0.6755    0.5739
    0.4630    0.6755   -0.5739
   -0.4630    0.6755   -0.5738
    0.4630   -0.6755   -0.5739
    0.4630   -0.6755    0.5739
         0         0    1.0000
   -0.4630   -0.6755    0.5739
   -0.4630    0.6755    0.5738];
    
    
    B=[0     1
     1     0
     1     0
     1     0
     1     0
     1     0
     1     0
     0     1
     1     0
     1     0];
 
 
    % Rays of cone in combined space
    load alphaRays;
    
    % Incidence matrix
    load alphaInc;
 

 
elseif strcmp(crys_choice,'asc') % Ascorbic 
 
    fprintf('========= ASCORBIC ========\n');
    
   
    A=[1.0000         0         0
   -1.0000         0         0
    0.2098         0    0.9777
   -0.2098         0   -0.9777
    0.1623         0   -0.9867
   -0.1623         0    0.9867
    0.3516   -0.9361         0
   -0.3516    0.9361         0
    0.3516    0.9361         0
   -0.3516   -0.9361         0];

    B=[1     0     0
     1     0     0
     0     1     0
     0     1     0
     0     1     0
     0     1     0
     0     0     1
     0     0     1
     0     0     1
     0     0     1]; 
 
    % Rays of cone in combined space
    load ascRays;
    
    % Incidence matrix
    load ascInc;
 
else
    error('Wrong choice');
end

% Parameter bounds
tmax=10;tmin=0;


%% Compute partiton using Bemporad's algo

pp.A=A;
pp.B=B;
pp.m=size(pp.A,1);
pp.d=size(pp.A,2);
pp.prm.dim=size(pp.B,2);

pp.H=[eye(pp.prm.dim);...
      -1*eye(pp.prm.dim)];
pp.h=[tmax*ones(pp.prm.dim,1);...
      -1*tmin*ones(pp.prm.dim,1)];

fprintf('Compute Bemporad partition...\n');
part=computePartition(pp);
fprintf('...done !\n');
  
bemp=[];
for i=1:length(part)
    bemp=[bemp,part(i).P];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute partition using faces projection algo

poly.A=[A,-1*B];
poly.b=zeros(size(poly.A,1),1);
poly.dim=size(poly.A,2);
poly.prm.dim=size(B,2);
% Rays
poly.R=R;

%% Polish rays

for i=1:poly.dim*size(poly.R,2)
    if abs(poly.R(i))<1e-4
        poly.R(i)=0;
    end
end

%% Parameter set
poly.prm.P=Polyhedron([eye(poly.prm.dim);-1*eye(poly.prm.dim)],...
            [tmax*ones(poly.prm.dim,1);-1*tmin*ones(poly.prm.dim,1)]);
poly.prm.H=poly.prm.P.H(:,1:end-1);
poly.prm.K=poly.prm.P.H(:,end);
    
%% Compute face lattice

fprintf('Compute face lattice...\n');
Latt=cmpLatt(poly);
fprintf('...done !\n');

%% Compute rays saturated by p-faces

for i=1:size(Latt{poly.prm.dim+1}.f,1) 
    
    pFcs{i}.f=Latt{poly.prm.dim+1}.f(i,:);
    pFcs{i}.R=[];
    pFcs{i}.pR=[];
    
    % Compute rays saturated by the p-face
    for j=1:size(poly.R,2)
        
        if abs(poly.A(find(pFcs{i}.f)',:)*poly.R(:,j))<1e-9
            
            pFcs{i}.R=[pFcs{i}.R,poly.R(:,j)];
            
            % Project rays onto parameter space
            pFcs{i}.pR=[pFcs{i}.pR,poly.R(end-poly.prm.dim+1:end,j)];
            
        end
    end
    %fprintf('Number of rays saturating face %i: %i\n',i,size(pFcs{i}.R,2));

end


%% Compute regions of the parameter space

cnt=1;

for i=1:size(Latt{poly.prm.dim+1}.f,1)
    
    if rank(pFcs{i}.pR)==poly.prm.dim
        Rgs{cnt}.pR=pFcs{i}.pR;
        cnt=cnt+1;
    end
end


%% Simplify rays sets

fprintf('Start simplifying ray sets...\n');

for i=1:length(Rgs)
    
    enc=[];
    for j=1:size(Rgs{i}.pR,2)
        
        isEnc=false;
        
        for k=1:size(enc,2)
            if norm(Rgs{i}.pR(:,j)-enc(:,k),2)<1e-4 
                isEnc=true;
                break;
            end
        end
        
        if isEnc==false
            enc=[enc,Rgs{i}.pR(:,j)];
        end
    end

    Rgs{i}.pR=enc;
end

fprintf('...simplification done !\n');


%% Transform conic regions into H-format

for i=1:length(Rgs)

    aux.R=Rgs{i}.pR';
    Rgs{i}.C=Polyhedron(aux);
    Rgs{i}.C.minHRep;
    Rgs{i}.C.normalize;
    % Compute H-representation
    Rgs{i}.C=Rgs{i}.C.computeHRep;
    % Intersect with the domain of parameters
    Rgs{i}.I=Polyhedron([poly.prm.H;Rgs{i}.C.H(:,1:poly.prm.dim)],...
                        [poly.prm.K;Rgs{i}.C.H(:,end)]);
    Rgs{i}.I.minHRep;
end


%% Count number of different regions

cnt=0;
enc=[];

for i=1:length(Rgs)

    isEnc=false;
    
    % Check whether current element is in the list of encountered elements
    for j=1:length(enc)
        %if Rgs{i}.I==enc{j}.I
        if checkRays(Rgs{i}.pR,enc{j}.pR)
            isEnc=true;
            break;
        end
    end
    
    if isEnc==false
        cnt=cnt+1; 
        enc{cnt}.I=Rgs{i}.I;
        enc{cnt}.pR=Rgs{i}.pR; 
    end
    
end

%fprintf('Number of different elements:%i\n',cnt);

% for i=1:length(enc)
%     
%       figure;clf;
%       plot(enc{i}.I);
%     
%       fprintf('Vertices of %i-th part:\n',i);
%       disp(enc{i}.I.V);
%       
%       disp('Rays:');
%       disp(enc{i}.pR);
%       
%       pause;
% end

aux=[];

for i=1:length(enc)
   aux=[aux,enc{i}.I];
end

U=PolyUnion(aux);
U=U.convexHull;

% Hypercube
Hcube=Polyhedron([eye(poly.prm.dim);-1*eye(poly.prm.dim)],[tmax*ones(poly.prm.dim,1);zeros(poly.prm.dim,1)]);

%disp('Union is equal to hypercube:');
%disp(U==Hcube);

%% Build partition

rCts=[];

% Initialize level-0 sets on enc
for i=1:length(enc)
    
    L{1}.S{i}.P=enc{i}.I;
    
    L{1}.S{i}.ints=zeros(1,0);
    
    % Check if set i has fulldim intersections with others 
    for j=1:length(enc)
    
        if j~=i
            
            iP=intersect(enc{i}.I,enc{j}.I);
            
            if iP.isFullDim
                L{1}.S{i}.ints(end+1)=j;
            end
        end
        
    end
    
%     disp('Intersections:');
%     disp(L{1}.S{i}.ints);
   
end


cnt=1;

while 1

    fprintf('Iteration %i\n',cnt);
    
    L{cnt+1}.S=[];
    np=0;
    
    for i=1:length(L{cnt}.S)
        
        ints=L{cnt}.S{i}.ints;
        
        % If does not intersect anybody add to critical regions
        if isempty(ints)
            rCts=[rCts,L{cnt}.S{i}.P];
        else
            
            for j=1:length(ints)
                
                jj=ints(j);
                
                inP=intersect(L{cnt}.S{i}.P,L{cnt}.S{jj}.P);
                cP1=L{cnt}.S{i}.P\inP;
                cP2=L{cnt}.S{jj}.P\inP;
               
                
                % Add inP if not already in L{cnt+1}.S
                if inP.isFullDim
                    isIn=false;
                    for k=1:length(L{cnt+1}.S)
                        if inP==L{cnt+1}.S{k}.P
                            isIn=true;
                            break;
                        end
                    end
                    
                    if isIn==false
                        % Add to L{cnt+1}.S
                        L{cnt+1}.S{np+1}.P=inP;
                        np=np+1;
%                         figure;
%                         plot(inP);
%                         title('inP');
%                         pause;
                    end
                end
                
                
                % Add cP1 if not already in L{cnt+1}.S 
                for r=1:length(cP1)
                    if cP1(r).isFullDim
                        isIn=false;
                        for k=1:length(L{cnt+1}.S)
                            if cP1(r)==L{cnt+1}.S{k}.P
                                isIn=true;
                                break;
                            end
                        end
                    
                        if isIn==false
                            % Add to L{cnt+1}.S
                            L{cnt+1}.S{np+1}.P=cP1(r);
                            np=np+1;
%                           figure;
%                           plot(cP1);
%                           title('cP1');
%                            pause;
                        end
                    end
                end
                
                % Add cP2 if not already in L{cnt+1}.S
                for r=1:length(cP2)
                    if cP2(r).isFullDim
                        isIn=false;
                        for k=1:length(L{cnt+1}.S)
                            if cP2(r)==L{cnt+1}.S{k}.P
                                isIn=true;
                                break;
                            end
                        end
                    
                        if  isIn==false
                            % Add to L{cnt+1}.S
                            L{cnt+1}.S{np+1}.P=cP2(r);
                            np=np+1;
%                           figure;
%                           plot(cP2);
%                           title('cP2');
%                           pause;
                        end
                    end
                end
                
            end
            
            
        end
      
    end
    
    % Compute intersections for sets at level cnt+1
    for i=1:length(L{cnt+1}.S)
        
        ints=zeros(1,0);
        
        for j=1:length(L{cnt+1}.S)
            if j~=i
                inP=intersect(L{cnt+1}.S{i}.P,L{cnt+1}.S{j}.P);
                
                if inP.isFullDim
                    ints(end+1)=j;
                end
            end
        end
        
        L{cnt+1}.S{i}.ints=ints;
        
%         figure;
%         plot(L{cnt+1}.S{i}.P);
%         pause;
    end
    
    
    
    
    if ~isempty(rCts)
        % Take union of critical regions
        U=PolyUnion(rCts);
        U=U.convexHull;
    
        % Check if union is equal to hybercube
        if U==Hcube
            fprintf('Union of sets at level %i equals hypercube\n',cnt+1);
            break;
        end
    end
    
    cnt=cnt+1;
   
end

%% Check intersections

for i=1:length(rCts)
        
    ints=zeros(1,0);
        
    for j=1:length(rCts)
        if j~=i
            inP=intersect(rCts(i),rCts(j));
                
            if inP.isFullDim
                ints(end+1)=j;
            end
        end
    end
        
    disp('Intersections:');
    disp(ints);
        
end

%% Compare Bemporad's output with ours

% figure(1);clf;
% plot(bemp);
% title('Bemporad');
% 
% figure(2);clf;
% plot(rCts);
% title('Face projections output');

%% Extract parameterization matrices from bemporad's output and ours

% Bemporad

fprintf('Compute Bemporad matrices...\n');

bmpMts=[];

for k=1:length(bemp)
    
    tc=bemp(k).chebyCenter.x;
    
    pl=Polyhedron(pp.A,pp.B*tc);
    
    
    for i=1:size(pl.V,1)
        vrt=pl.V(i,:)';
        
        inqs=find(abs(pp.A*vrt-pp.B*tc)<1e-6);
        % Take only the first d (degeneracy)
        inqs=inqs(1:pp.d);
        bmpMts{k}.M{i}=inv(pp.A(inqs,:))*pp.B(inqs,:);
    end
    
end

fprintf('...computed !\n');

% Ours

fprintf('Compute our matrices...\n');

fpMts=[];

for k=1:length(rCts)
    
    tc=rCts(k).chebyCenter.x; 
    
    pl=Polyhedron(pp.A,pp.B*tc);
   
    
    for i=1:size(pl.V,1)
    
        vrt=pl.V(i,:)';
        
        inqs=find(abs(pp.A*vrt-pp.B*tc)<1e-6);
        % Take only the first d (degeneracy)
        inqs=inqs(1:pp.d);
        fpMts{k}.M{i}=inv(pp.A(inqs,:))*pp.B(inqs,:);
    end
   
end


fprintf('...computed !\n');

%% Cluster critical regiosn that corrspond to same parameterization

enc=cell(1,length(rCts));
chk=1:length(rCts);
cnt=1;

while 1
    
    if cnt==length(rCts)+1
        break;
    end
    
    fprintf('Region %i\n',cnt);
    
    % Find critical region that has same parameterization as rCts(cnt)
    for k=1:length(fpMts)
        
        ok=false;
        
        if k~=cnt && length(fpMts{k}.M)==length(fpMts{cnt}.M)  
           
            allIn=true;
            
            % Check all matrices of cnt are in matrices list of k 
            for i=1:length(fpMts{cnt}.M)
                
                isIn=false;
                
                for j=1:length(fpMts{k}.M)
                    
                    if sum(sum(fpMts{cnt}.M{i}==fpMts{k}.M{j}))==pp.d*pp.prm.dim
                        isIn=true;
                        break;
                    end
                end
                
                if isIn==false
                    allIn=false;
                    break;
                end
            end
            
            ok=allIn;
        end
        
        if ok
            if ~isempty(enc{cnt})
                enc{cnt}.idx=[enc{cnt}.idx,k];
            else
               
                enc{cnt}.idx=[k];
            end
        end
    end
   
    
    % If non-empty, check if union equals whole set
    if ~isempty(enc)
        
        tst=[];
        for k=1:length(enc)
            if ~isempty(enc{k})
                tst=[tst,enc{k}.idx];
            end
        end
        
        allIn=true;
        
        for k=1:length(rCts)
            if isempty(find(tst==k))
                allIn=false;
                break;
            end
        end
        
        if allIn
            fprintf('First step of Clustering done\n');
            break;
        end
        
    end
    
    cnt=cnt+1;
end
 

% Check whether enc is empty or not
empt_enc=true;
for i=1:length(enc) 
    if ~isempty(enc{i})
        empt_enc=false;
        break;
    end
end

if ~empt_enc

    new_rCts=[];
    cnt=1;

    while 1

        fprintf('Iteration %i\n',cnt);
    
        if cnt==length(enc)+1
            break;
        end
    
        nP=PolyUnion([rCts(cnt),rCts(enc{cnt}.idx)]);
        nP=nP.convexHull;
    
        isIn=false;
    
        % Check if not already in
        if ~isempty(new_rCts)
        
            for i=1:length(new_rCts)
                if new_rCts(i)==nP
                    isIn=true;
                    break;
                end
            end
        
        end
    
        if ~isIn
            new_rCts=[new_rCts,nP]; 
        end
    
    
        % Evaluate convex hull of union
    
        U=PolyUnion(new_rCts);
        U=U.convexHull;
        
        if U==Hcube
            fprintf('Union equals Hcube: clustering finished\n');
            break;
        end
    
        cnt=cnt+1;
    end

else
    new_rCts=rCts;
end
    
%% Plot partitions in 3d
ft_size=50;

if pp.prm.dim<=3

    figure(1);clf;
    plot(bemp);
    title('Bemporad partition');

    pos.x=0;
    pos.y=0;
    pos.w=1000;
    pos.h=1000;
    
    
    fig1=figure('Color','w');clf;
    set(fig1,'Position',[pos.x,pos.y,pos.w,pos.h]);
    h=plot(new_rCts); 
    set(h,'LineWidth',1.5);
    axis equal;
    xlabel('\textbf{$t_1$}','Interpreter','LaTex','fontSize',ft_size);
    ylabel('\textbf{$t_2$}','Interpreter','LaTex','fontSize',ft_size);
    zlabel('\textbf{$t_3$}','Interpreter','LaTex','fontSize',ft_size);
    h=get(gca,'zlabel');
    set(h,'rotation',0);
    set(gca,'fontSize',ft_size,'fontname','Times New Roman');
    set(gcf, 'renderer','painters')
    grid off;
    fprintf('Arrange before printing to pdf...\n');
    pause;
    myaa;
    save2pdf([dir,'partAceta.pdf'],gcf,600);
      

else
    
    % First compute projections
    proj_bemp=[];
    for i=1:length(bemp)
        proj_bemp=[proj_bemp,projection(bemp(i),1:3,'fourier')];
    end
    
    proj_new_rCts=[];
    for i=1:length(new_rCts)
        proj_new_rCts=[proj_new_rCts,projection(new_rCts(i),1:3,'fourier')];
    end
    
    figure(1);clf;
    plot(proj_bemp);
    title('Bemporad partition');

    figure(2);clf;
    plot(proj_new_rCts); 
    title('Our partition');
    
end


%% Plot different polytopes for each partition polytope 

for i=1:length(new_rCts)
    
    tcur=new_rCts(i).chebyCenter.x;
    Pcur=Polyhedron(pp.A,pp.B*tcur);
    
    % Plot partition polytope
    fig1=figure('Color','w');clf;
    set(fig1,'Position',[pos.x,pos.y,pos.w,pos.h]);
    h=plot(new_rCts(i),struct('shade',0.5,'edgecolor',[.2 .2 .2],'linestyle','-','linewidth',5,'color','w','DisplayName','Fitted polytope'));
    set(h,'LineWidth',5);
    axis equal;axis off;
    set(gcf, 'renderer', 'opengl')
    grid off;
    fprintf('Arrange before printing to pdf...\n');
    pause;
    myaa;
    save2pdf([dir,'acePrt_',num2str(i),'.pdf'],gcf,600);
    
    % Plot sample polytope
    fig2=figure('Color','w');
    set(fig1,'Position',[pos.x,pos.y,pos.w,pos.h]);
    h=plot(Pcur,struct('shade',0.5,'edgecolor',[.2 .2 .2],'linestyle','-','linewidth',5,'color','w','DisplayName','Fitted polytope'));
    set(h,'LineWidth',5);
    axis equal;axis off;
    set(gcf, 'renderer', 'opengl')
    grid off;
    fprintf('Arrange before printing to pdf...\n');
    pause;
    myaa;
    save2pdf([dir,'aceSmp_',num2str(i),'.pdf'],gcf,600);
    
    
    
    pause;
end



















