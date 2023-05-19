function fea=NormalizeFea(fea)
[m,n]=size(fea);
for i=1:m
    fea(i,:)=fea(i,:)./max(1e-12,norm(fea(i,:)));
end
end