function [ACCres,NMIres]=TSNMF(X,W,para,gnd,Z0,V,Vexl)

% initialize parameters
iter=para.iter;

HyperParameters.lamda=para.lamda;

SNMF_para.m=para.m;
SNMF_para.beta=para.beta;
SNMF_para.k=para.k;
SNMF_para.maxIter=para.maxIter;

Tensor_para.DEBUG = 1;
Tensor_para.tol = 1e-3;

Omega = find(Z0~=0);

S=W;
it=0;
Z_tensor=Z0;
while it<iter
    it=it+1;
    % adjust S with Z
    S(Z_tensor>0)=1-(1-Z_tensor(Z_tensor>0)).*(1-S(Z_tensor>0));
    S(Z_tensor<0)=(1+Z_tensor(Z_tensor<0)).*S(Z_tensor<0);
    S=(abs(S)+abs(S'))/2;
    
    % update Vs
    [H,Vs,a,h]=Enhanced_SNMF(S,SNMF_para,V,Vexl);
    [ACCres,NMIres]=cal_ACC_NMF_symNMF_v3(Vs,gnd);
    
    % update A and Z
    [A_tensor,Z_tensor,E]=tensor_low_nuclear(X,Vs*Vs',Z0,Omega,Tensor_para,HyperParameters);
    
    % update S
    S=A_tensor;
end
end
