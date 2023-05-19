function [V_final,K,Qf] = PCPSNMF(W,mu,k,Z,alpha,Vb)

maxIter=500;
nIter=0;
tryNo=0;
obj=1e+7;
n=size(W,2);

graphD=zeros(size(W,1),size(W,1));
GD=sum(W,2);
graphD=spdiags(GD,0,n,n);
graphL=graphD-W;
D_mhalf = spdiags(GD.^-.5,0,n,n) ;
graphL = D_mhalf*graphL*D_mhalf;

Kb=rand(n);
[Vb] = NormalizeK(Vb);
[Kb] = NormalizeK(Kb);

for i=1:n
    tag=find(Z(i,:)==1);
    Kb(i,tag)=1;
end

W2=(Kb+Kb')/2;

FRb=norm(W2-Vb*Vb')+mu*trace(graphL*Kb)+alpha*norm(Kb-Z);

while nIter<maxIter
    % -----update K

    VV=Vb*Vb';

    T2=(2*alpha+1)*Kb+Kb'+mu*graphD;
    T1=2*VV+2*alpha*Z'+mu*W;

    K=Kb.*(T1./T2);    
    Kb=K;

    W2=(Kb+Kb')/2;

    % -----update V

    C=W2*Vb+W2'*Vb;
    D=2*Vb*Vb'*Vb;
    E=C./D;
    V=Vb.*(E.^(1/4));
    Vb=V;
    nIter=nIter+1;
end
FR=norm(W2-Vb*Vb')+mu*trace(graphL*Kb)+alpha*norm(Kb-Z);

if FR < obj
    [V_final] = NormalizeK(V);
    obj=FR;
    nIter=0;

end
end

function [K] = NormalizeK(K)
[m,n]=size(K);
for i=1:m
    K(i,:)=K(i,:)./max(1e-12,norm(K(i,:)));
end
end