function [H,Vexl,a,h,obj]=Enhanced_SNMF(S,para,V,Vs)
%min{a,V,Vs}sum(ai*||S-V*V'||_{F}^2)+sum(ai*||V-Vs||_{F}^2)+beta*||a||_{F}^2
%——————————————————————————————
%input:
%     S   n*n ---------------------------------------affinity matrix
%     V   n*k*m -------------------------------------a set of initial embedding matrices
%     Vs  n*k ---------------------------------------initial embedding matrix
%Output:
%     H   n*k*m -------------------------------------a set of processed embedding matrices
%     Vexl  n*k -------------------------------------processed embedding matrix
%     a   m*1 ---------------------------------------weight
%     h   1*m ---------------------------------------residual

m=para.m;
beta=para.beta;
maxIter=para.maxIter;
a=rand(m,1);
for i=1:m
    H{i}=V(:,:,i);
end
Vexl=Vs;
y0=zeros(m,1);
h=zeros(1,m);
obj{1}=h*a+beta*norm(a).^2;

for iter=1:maxIter
    % update V
    for i=1:m
        NUM{i}=2*S*H{i}+Vexl;
        DEN{i}=2*H{i}*H{i}'*H{i}+H{i}+eps;
        H{i}=H{i}.*(NUM{i}./DEN{i});
        h(1,i)=norm(S-H{i}*H{i}').^2+norm(H{i}-Vexl).^2;
        y0(i,:)=(-0.5)*h(1,i)/beta;
    end
    % update a
    [a,res]=projection(y0,1);
    % update Vs
    V0=zeros(size(Vs));
    for i=1:m
        V0=V0+a(i,1).*H{i};
    end
    Vexl=V0;
    
    obj{iter+1}=h*a+beta*norm(a).^2;
    
    %disp(['the ', num2str(iter), ' obj is ', num2str(obj(iter))]);
    if (iter>300 && abs(obj{iter+1}-obj{iter})<10^-3)
        break;
    end
end
obj=cell2mat(obj);
plot(obj);
for i=1:m
    H{i}=NormalizeFea(H{i});
end
Vexl=NormalizeFea(Vexl);

end