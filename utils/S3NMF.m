function [V_final]=S3NMF(W,Z,para,V)

lamda1=para.lamda1;
lamda2=para.lamda2;
alpha=para.alpha;
mu=para.mu;
maxIter=para.maxIter;

n=size(W,2);
S=rand(n,n);
t=1;
GD=sum(W,2);
graphD=spdiags(GD,0,n,n);
graphL=graphD-W;

while t<maxIter
    %update S
    Snum=2*V*V'+alpha*W+2*mu*Z;
    Sden=S'+alpha*graphD+(1+2*mu+lamda2)*S+lamda1+eps;
    S=S.*(Snum./Sden);
    Sx=(S+S')/2;
    %update V
    Vnum=S*V+S'*V;
    Vden=2*V*V'*V+eps;
    M=(Vnum./Vden).^.25;
    V=V.*M;
    
    obj(t)=norm(Sx-V*V',"fro").^2+alpha*trace(graphL*S)+lamda1*norm(S,1)+lamda2/2*norm(S,"fro").^2;
    t=t+1;
end
V_final=V;