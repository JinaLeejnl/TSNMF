function indices=findlabels(per,labels,Class)
indices=[];
for i=1:Class
    pcd=find(labels==i);
    len=length(pcd);
    p=randperm(len);
    indices=[indices;pcd(p(1:per))];
end
end