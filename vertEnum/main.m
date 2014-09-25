clear;

%%%%%%
%% To run statictics on artificial images of your choice:
%% 'cube','ibuprofen','aceta','aspirin','lglubeta','lglualpha','asc'
%%%%%

clear;

chs=['cube','ibuprofen','aceta','aspirin','lglubeta','lglualpha','asc'];

% Crystal choice
crys_choice='ibuprofen';
dat.crys_choice=crys_choice;

disp('crys_choice=');disp(crys_choice);

% Number of stats
stats.N=1;

fprintf('Number of samples: %i\n',stats.N);

% Bounds on scaling parameter
dat.tmin=0;
dat.tmax=10;


if strcmp(crys_choice,'cube') % cube
    
    fprintf('======== CUBE ========\n');
    
    A=[eye(3);-eye(3)];
    
    B=[eye(3);eye(3)];
    
elseif strcmp(crys_choice,'ibuprofen') % ibuprofen
    
    fprintf('======= IBUPROFEN =======\n');

    %B=[eye(3);0,0,1]; B=[B;B];
    %A=[1,0,0;0,-2,0;0,1,1;0,1,-1]; A=[A;-A];
     
%     A=[ 0.1311   -0.5922   -0.7950
%     1.0000         0         0
%    -0.1627         0    0.9867
%    -0.1311   -0.5922    0.7950
%    -0.1311    0.5922    0.7950
%    -1.0000         0         0
%     0.1627         0   -0.9867
%     0.1311    0.5922   -0.7950]

%  B=[ 0     0     1
%         1     0     0
%         0     1     0
%         0     0     1
%         0     0     1
%         1     0     0
%         0     1     0
%         0     0     1];

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
    
    
%     A=[0    1.0000         0
%        0   -1.0000         0
%     0.5969         0    0.8023
%    -0.5969         0    0.8023
%    -0.5969         0   -0.8023
%     0.5969         0   -0.8023
%          0   -0.7803   -0.6254
%          0   -0.7803    0.6254
%          0    0.7803    0.6254
%          0    0.7803   -0.6254];
%     
% 
%     B=[1     0     0
%      1     0     0
%      0     1     0
%      0     1     0
%      0     1     0
%      0     1     0
%      0     0     1
%      0     0     1
%      0     0     1
%      0     0     1];

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

% elseif choice==chs(7) % Ascorbic
%  
%     A=[ 1.0000         0         0
%    -1.0000         0         0
%    -0.2098         0    0.9777
%     0.2098         0   -0.9777
%     0.9479         0   -0.3187
%    -0.9479         0    0.3187
%     0.9387   -0.3447         0
%    -0.9387    0.3447         0
%     0.9387    0.3447         0
%    -0.9387   -0.3447         0];
%     
%     
%     B=[1     0     0
%      1     0     0
%      0     1     0
%      0     1     0
%      0     1     0
%      0     1     0
%      0     0     1
%      0     0     1
%      0     0     1
%      0     0     1];
 
elseif strcmp(crys_choice,'asc') % Ascorbic 
 
    fprintf('========= ASCORBIC ========\n');
    
    
%     A=[1.0000         0         0
%    -1.0000         0         0
%    -0.2098         0    0.9777
%     0.2098         0   -0.9777
%     0.9479         0   -0.3187
%    -0.9479         0    0.3187
%     0.9387   -0.3447         0
%    -0.9387    0.3447         0
%     0.9387    0.3447         0
%    -0.9387   -0.3447         0];
%     
%     B=[1     0     0     0     0
%      1     0     0     0     0
%      0     1     0     0     0
%      0     1     0     0     0
%      0     0     1     0     0
%      0     0     1     0     0
%      0     0     0     1     0
%      0     0     0     1     0
%      0     0     0     0     1
%      0     0     0     0     1 ];

    
%     A=[1.0000         0         0
%    -1.0000         0         0
%    -0.2098         0    0.9777
%     0.2098         0   -0.9777
%     0.9479         0   -0.3187
%    -0.9479         0    0.3187
%     0.9387   -0.3447         0
%    -0.9387    0.3447         0
%     0.9387    0.3447         0
%    -0.9387   -0.3447         0];
% 
% 
%     B=[1.0000         0         0
%        1.0000         0         0
%        0    3.5000         0
%        0    3.5000         0
%        0    1.0000         0
%        0    1.0000         0
%        0         0    1.0000
%        0         0    1.0000
%        0         0    1.0000
%        0         0    1.0000];


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

%     % Simpler ascorbic
%     A=[1.0 0.0 0.0;
%        -1.0 0.0	0.0;
%        -0.20978921507880458	0.0	0.9777466365253419;
%        0.20978921507880458	0.0	-0.9777466365253419;
%        0.9387002229335936	-0.3447345231688028	0.0;
%       -0.9387002229335936	0.3447345231688028	0.0;
%       0.9387002229335936	0.3447345231688028	0.0;
%      -0.9387002229335936	-0.3447345231688028	0.0];
%     
%     B=[1.0	0.0	0.0
%        1.0	0.0	0.0
%        0.0	1.0	0.0
%        0.0	3.5	0.0
%        0.0	0.0	1.0
%        0.0	0.0	1.0
%        0.0	0.0	1.0
%        0.0	0.0	1.0];
 
else
    error('Wrong choice');
end
  
dat.d=3; % Polytope dimension
dat.n=size(A,1); % Number of facets
dat.m=size(B,2); % Dimension of scaling parameter

% Good t's for each crystal

% 'cube','ibuprofen','aceta','aspirin','lglubeta','lglualpha','asc','asc5'
    
t_cube=[1;
        1;
        1];

t_ibu=[2;
       2;
       3];

    
t_aceta=[1;
         2;
         1];

t_asp=[3;
       3;
       7;
       8];

t_beta=[1;
        5;
        2];

t_alp=[1;
       1.3];

%t_asc=[1;
%       7;
%       2;
%       4;
%       4];

t_asc=[1;
       2;
       4];

% % Simpler ascorbic
% t_asc=[1;
%        4;
%        4];

   
% Select t

if strcmp(crys_choice,'cube') % Cube
    t=t_cube;
elseif strcmp(crys_choice,'ibuprofen') % Ibuprofen  
    t=t_ibu;
elseif strcmp(crys_choice,'aceta') % Acetaminophen
    t=t_aceta;
elseif strcmp(crys_choice,'aspirin') % Aspirin
    t=t_asp;
elseif strcmp(crys_choice,'lglubeta') % Lglubeta
    t=t_beta;
elseif strcmp(crys_choice,'lglualpha') % Lglualpha
    t=t_alp;
elseif strcmp(crys_choice,'asc') % Ascorbic 
    t=t_asc;
else
    error('Wrong scaling parameter');    
end
    
disp('Right scaling=');disp(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate partitions from crystal model

stats.ts=zeros(dat.m,stats.N);
    
fprintf('Generate partition from crystal model...\n');

poly.Inc=Inc; % Incidence matrix

% Turn Inc into strings
auxInc=[];
for j=1:size(poly.Inc,2)
    str=[];
    for i=1:size(poly.Inc,1)
        str=[str;...
             num2str(poly.Inc(i,j))];
    end
    auxInc=[auxInc,str];
end
poly.Inc=auxInc;

poly.R=R; % Rays of cone in combined space
poly.A=A;
poly.B=B;
poly.numInqs=size(poly.A,1);
poly.prm.dim=dat.m;
poly.dim=size(poly.A,2);
poly.H=[eye(dat.m);-eye(dat.m)];
poly.h=[dat.tmax*ones(dat.m,1);-dat.tmin*ones(dat.m,1)];

part_1=computePartition(poly);

return;

% Compute vertices

for i=1:length(part_1)

    inqs=part_1(i).inqs;

    for j=1:length(inqs)
       
        inq=inqs(j,:);
        
        Prt.V{i}.v{j}=inv(poly.A(inq',:))*poly.B(inq',:);
    end
    
    Prt.V{i}.nV=length(inqs);
   
end
    
fprintf('Number of partition polytopes : %i\n',length(part_1));

% figure(1);clf;
% plot(part_1.P);
% title('Partition obtained by Bemporad algo');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use face projections for computing partition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

part_2=computeParamVert(poly);





























