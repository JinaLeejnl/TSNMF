clc;
clear all;
load 'BinAlpha.mat'

fea=NormalizeFea(fea');
fulldata=fea';
Trans_fulldata=fea;

for k=12:3:36
[datasub{k},subN{k},labelsub{k},Trans_datasub{k}]=loadsub(fulldata,gnd,k);
%% W
nnparams=cell(1);
nnparams{1}='knn';
opts.K =floor(log2(subN{k})); 
opts.maxblk = 1e7;
opts.metric = 'eucdist';
nnparams{2}=opts;
T_G=slnngraph(datasub{k},[],nnparams);

W{k}=zeros(subN{k});
for i = 1:subN{k}
    tmp = find(T_G(:,i)~=0);
    W{k}(tmp,i) = exp(-0.25*T_G(tmp,i));  
end
W{k}=(W{k}+W{k}')/2;
GD=full(sum(W{k},2));
D_mhalf{k}=spdiags(GD.^-0.5,0,subN{k},subN{k});
graphL{k}=D_mhalf{k}*W{k}*D_mhalf{k};
imagesc(W{k})

fprintf('#of classes: %g\n',k);
for nter=1:10
    V_else{k}{nter}=rand(subN{k},k);
    Vexl{k}{nter}=rand(subN{k},k);
    V{k}{nter}=rand(subN{k},k,10);

    sumof=ceil(subN{k}*subN{k}*0.1);
    every=subN{k}/k;
    row=randi(subN{k},1,sumof);
    col=randi(subN{k},1,sumof);
    A{k}{nter}=zeros(subN{k});        
    for i=1:sumof
        if labelsub{k}(row(1,i),1)==labelsub{k}(col(1,i),1)  
            A{k}{nter}(row(1,i),col(1,i))=1;
        else
            A{k}{nter}(row(1,i),col(1,i))=-1;
        end
    end
    Z{k}{nter}=A{k}{nter};
    Z{k}{nter}(A{k}{nter}<0)=0;
end
end