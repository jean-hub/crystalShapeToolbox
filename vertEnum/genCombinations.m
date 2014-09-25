% Generate the list of all combinations of p elements among n

function [all] = genCombinations(n,p)
    all.mat=zeros(nchoosek(n,p),p);
    all.count=1;
    if n==p
        all.mat(1,:)=1:n;
    else
        list=zeros(1,p);
        index=1;
        all=computeList(index,n,p,list,all);
    end
end

function [all] = computeList(index,n,p,list,all)
    if index>=p+1
        all.mat(all.count,:)=list;
        all.count=all.count+1;
        return;
    end
    start=1;
    if index>1
        start=list(index-1)+1;
    end
    for i=start:n
        list(index)=i;
        all=computeList(index+1,n,p,list,all);
    end 
end