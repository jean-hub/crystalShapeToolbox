function obj = processData(obj)
%
% Process data d12,d13 to get e_d12,e_d13 (Douglas-Peucker line simplification)
%=============================================================================
%
obj.e_d12 = splitAndMerge(obj.d12,0.15);
obj.e_d13 = splitAndMerge(obj.d13,0.15);

aux_12 = Polyhedron((obj.e_d12)');aux_12.minVRep();
obj.e_d12 = (aux_12.V)';
obj.Ad12 = aux_12.A;
obj.bd12 = aux_12.b;
    
aux_13 = Polyhedron((obj.e_d13)');aux_13.minVRep();
obj.e_d13 = (aux_13.V)';
obj.Ad13 = aux_13.A;
obj.bd13 = aux_13.b;

obj.e_nD12 = size(obj.e_d12,2);
obj.e_nD13 = size(obj.e_d13,2);
end


function DougLine = splitAndMerge(PntList,eps)
%
% Implementation of the Split-and-merge algo
%============================================
%

% Find the point furthest away from the segment [p_fir,p_las]
nPts=size(PntList,2);
seg=[PntList(:,1),PntList(:,nPts)];
d_max=0;
ind_max=1;
for i=2:size(PntList,2)-1
    d=cmpDstSeg(PntList(:,i),seg);
    if d>d_max
        d_max=d;
        ind_max=i;
    end
end

if d_max<eps
    DougLine=seg;
else
    % Split data
    DougLine_1=splitAndMerge(PntList(:,1:ind_max),eps);
    DougLine_2=splitAndMerge(PntList(:,ind_max:end),eps);
    DougLine=zeros(size(DougLine_1,1),size(DougLine_1,2)-1+size(DougLine_2,2));
    DougLine(:,1:size(DougLine_1,2)-1)=DougLine_1(:,1:end-1);
    DougLine(:,size(DougLine_1,2):end)=DougLine_2;
    %DougLine=[DougLine_1(:,1:end-1),DougLine_2];
end

end


function d = cmpDstSeg(pt,seg)
%
% Computes the distance of a point to a segment
%==============================================
%

if ((pt-seg(:,1))'*(seg(:,2)-seg(:,1))>0)&&((pt-seg(:,2))'*(seg(:,1)-seg(:,2))>0)
    d=abs((seg(2,1)-seg(2,2))*pt(1)+(seg(1,2)-seg(1,1))*pt(2)+det(seg))^2 ...
        /((seg(:,2)-seg(:,1))'*(seg(:,2)-seg(:,1)));
else
    d=min([(pt-seg(:,1))'*(pt-seg(:,1)),(pt-seg(:,2))'*(pt-seg(:,2))]);
end
end