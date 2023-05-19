function [V]=SNMF(W,para,V)

maxIter=para.maxIter;
k=para.k;
subN=size(W,1);
t=0;
while t<maxIter
    %更新V
    tp=(W*V+W'*V)./(2*V*V'*V+eps);
    V=V.*(tp.^.25);

    t=t+1;
    %求E
    E(t)=norm(W-V*V');
end