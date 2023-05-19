function [revS, Dv, L] = MHPCP(X,Z,opts)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

p = opts.p;
mu = opts.mu;
view = 1;

revS = cell(1,view);
Dv = cell(1,view);
L = cell(1,view);


for i=1:view

    [H, W] = constructHW(X,p);

    invDe = diag(1./sum(H,2));

    S = H'*W*invDe*H;
    S=(S+S')/2;

    revS{i} = reviseHyperW(S,Z,mu);
    
    diag_revS = diag(revS{i});
    diag_revS_half = diag(diag_revS.^(-0.5));
    tmp_resS = diag_revS_half*revS{i}*diag_revS_half;
    revS{i} = (tmp_resS + tmp_resS') / 2;    
    
    Dv{i} = diag(sum(revS{i},1));
    L{i} = Dv{i} - revS{i};
    
end



function [H, W] = constructHW(fea, p)

% Each row of H consists of a hyper edge.
% Each column of H means the sample is used to consist of the corresponding
% hyper edge.


[nSmp, mFeat] = size(fea);
G = zeros(nSmp*(p+1),3);

dist = EuDist2(fea);
D = dist.^2;
sigma = mean(mean(dist));
%sigma=sqrt(2);
%sigma=2;
A = exp(D/(-sigma.^2));

dump = zeros(nSmp, p+1);
idx = dump;

for j = 1:p+1
    [dump(:,j),idx(:,j)] = min(dist,[],2);
    temp = (idx(:,j)-1)*nSmp+[1:nSmp]';
    dist(temp) = 1e100;
end

G(1:nSmp*(p+1),1) = repmat([1:nSmp]',[p+1,1]);
G(1:nSmp*(p+1),2) = idx(:);
G(1:nSmp*(p+1),3) = 1;

H = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
W = diag(sum(A.*H'));


function reW = reviseHyperW(W, Z, alpha)

cp_num = 30;

n = size(W,1);
gd = full(sum(W,2));
d_mhalf = spdiags(gd.^-.5,0,n,n);
L = d_mhalf*W*d_mhalf;


Fv = 2*rand(size(W))-1;
Fh = 2*rand(size(W))-1;

for i=1:cp_num
    Fv = alpha*L*Fv + (1-alpha)*Z;  
end

for i=1:cp_num
    Fh = alpha*Fh*L + (1-alpha)*Fv;  
end

F = Fh;

idx = (F>=0);
diagW = diag(W);
W_ii = repmat(diagW,1,n);

reW = ((W_ii-(1-F).*(W_ii-W)) .* idx) +  (((1+F).*W) .* (~idx));



function D = EuDist2(fea_a,fea_b,bSqrt)
%EUDIST2 Efficiently Compute the Euclidean Distance Matrix by Exploring the
%Matlab matrix operations.
%
%   D = EuDist(fea_a,fea_b)
%   fea_a:    nSample_a * nFeature
%   fea_b:    nSample_b * nFeature
%   D:      nSample_a * nSample_a
%       or  nSample_a * nSample_b
%
%    Examples:
%
%       a = rand(500,10);
%       b = rand(1000,10);
%
%       A = EuDist2(a); % A: 500*500
%       D = EuDist2(a,b); % D: 500*1000
%
%   version 2.1 --November/2011
%   version 2.0 --May/2009
%   version 1.0 --November/2005
%
%   Written by Deng Cai (dengcai AT gmail.com)


if ~exist('bSqrt','var')
    bSqrt = 1;
end

if (~exist('fea_b','var')) || isempty(fea_b)
    aa = sum(fea_a.*fea_a,2);
    ab = fea_a*fea_a';
    
    if issparse(aa)
        aa = full(aa);
    end
    
    D = bsxfun(@plus,aa,aa') - 2*ab;
    D(D<0) = 0;
    if bSqrt
        D = sqrt(D);
    end
    D = max(D,D');
else
    aa = sum(fea_a.*fea_a,2);
    bb = sum(fea_b.*fea_b,2);
    ab = fea_a*fea_b';

    if issparse(aa)
        aa = full(aa);
        bb = full(bb);
    end

    D = bsxfun(@plus,aa,bb') - 2*ab;
    D(D<0) = 0;
    if bSqrt
        D = sqrt(D);
    end
end






