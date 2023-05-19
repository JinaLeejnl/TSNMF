function [V,S,D,obj]=SANMF(W,P,Z1,Z2,para,V)
%%%%%%% function readme %%%%%%%%
% solves the DP and SP guided SymNMF model
% i.e., \min |S-VV^T|+alpha|D\odot VV^T|+beta Tr(DLD^T) + eta Tr(SLS^T)
% s.t., V\geq 0, D \geq 0, S \geq 0, D,S meet the constraints.
% which is reformulate as
% \min |S-VV^T|+alpha|D\odot VV^T|+beta Tr(DLD^T) + eta Tr(SLS^T) + gamma(norm(P\odot(S-Z1))+norm(P\odot(D-Z2)))
% s.t., V\geq 0, D\geq 0 and S \geq 0.

% W is the local similarity matrix for 
% A is the diagonal matrix for W
% L is the Laplacian matrix for W 
% P is the position of supervisory information position matric
% Z1 is must-link matrix
% Z2 is the cannot-link matrix 

%%%%%%% parameter setting
A=diag(sum(W,2));
L=A-W;
alpha=para.alpha;
beta=para.beta;
gamma=para.gamma;
eta=para.eta;
c=para.c; % number of clusters

maxiter=para.maxiter;

S=rand(size(W))+W;
D=rand(size(W));
% V=rand(length(W),c);

obj(1)=norm(S-V*V').^2+alpha*trace(D*L*D')+beta*trace(S*L*S')...
    +gamma*norm(P.*(D-Z2))+gamma*norm(P.*(S-Z1))+eta*sum(sum(D.*S));

for iter=1:maxiter
    % update S 
    S=S.*((2*gamma*P.*Z1+2*beta*S*W+2*V*V')./(2*S+eta*D+2*beta*S*A+2*gamma*P.*S+eps));
    % update D
    D=D.*((2*gamma*(P.*Z2)+2*alpha*D*W)./(eta*S+2*alpha*D*A+2*gamma*(P.*D)+eps));
    % update V 
    V=V.*((S*V+S'*V)./(4*V*V'*V+eps)).^0.25;
    
    % objective function value 
    obj(iter+1)=norm(S-V*V').^2+alpha*trace(D*L*D')+beta*trace(S*L*S')...
    +gamma*norm(P.*(D-Z2))+gamma*norm(P.*(S-Z1))+eta*sum(sum(D.*S));

    %disp(['the ',num2str(iter),' iteration. obj value: ',num2str(obj(iter+1))]);
    
    % Check convergence condition
    if max(max(abs(obj(iter)-obj(iter+1))))<10^-3 
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S=(S+S')/2;
    D=(D+D')/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



