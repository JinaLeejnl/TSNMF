function [A,Z,E]=tensor_low_nuclear(X,A0,Z0,Omega,opts,HyperParameters)
%min_{C,Z,A,D,E} ||C||_{*}+lamda*||E||_{F}^2
%s.t.A>=0, Z=C(:,:,1), A=C(:,:,2), Z=D
%
%--------------------------------------------------
%Input:
%       X       -    d*n matrix
%       A0       -    n*n matrix
%       Z0       -    n*n matrix
%       W       -    n*n matrix affinity
%       Omega   -    index of the observed entries
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%Output:
%       A       -    n*n matrix
%       Z       -    n*n matrix
tol = 1e-3; 
max_iter = 500;
rho = 1.1;
mu = 1e-3;
max_mu = 1e10;
DEBUG = 1;
if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end

lamda=HyperParameters.lamda;
% weight=HyperParameters.w;
% p=HyperParameters.p;

X_size = size(X);   
n2 = X_size(2); % number of samples

Ck = zeros(n2,n2,2);
Y2 = Ck;
A =rand(n2,n2);
Z = A;
D = A;
Y3 = D;
sZ = Z0(Omega);
D(Omega) = sZ;  

I = eye(n2);
E=zeros(n2);
Y1 = E;

for iter = 1 : max_iter
    Ak = A;
    Zk = Z;
    Ek = E;
    Dk = D;
    Ck(:,:,2) = A;  
    Ck(:,:,1) = Z;
    % update C
    [C,tnnC,~] = prox_tnn(Ck + Y2/mu,1/mu);
    
    %update Z
    Z=(mu*C(:,:,1)-Y2(:,:,1)+mu*D-Y3)/(2*mu*I);
    
    %update A
    Cpo2=(abs(C(:,:,2))+C(:,:,2))/2;
    Cne2=(abs(C(:,:,2))-C(:,:,2))/2;
    Y2po2=(abs(Y2(:,:,2))+Y2(:,:,2))/2;
    Y2ne2=(abs(Y2(:,:,2))-Y2(:,:,2))/2;
    Y1po=(abs(Y1)+Y1)/2;
    Y1ne=(abs(Y1)-Y1)/2;
    Epo=(abs(E)+E)/2;
    Ene=(abs(E)-E)/2;
    Anumerator=A0+Y1po/mu+Cpo2+Ene+Y2ne2/mu;
    Adenominator=2*I*A+Epo+Y2po2/mu+Y1ne/mu+Cne2;
    A=A.*(Anumerator./(Adenominator+eps));
    
    %update E
    E=(mu*A0+Y1-mu*A)/((2*lamda+mu)*I);
    
    %update D
    D=Z+Y3/mu;
    D(Omega) = sZ;
    
    dY1 = A0 - A - E;
    dY2 = Ck - C;
    dY3 = Z - D;
    chgC = max(abs(Ck(:)-C(:)));
    chgA = max(abs(Ak(:)-A(:)));
    chgZ = max(abs(Zk(:)-Z(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chgD = max(abs(Dk(:)-D(:)));
    chg = max([ chgC chgA chgZ chgE chgD max(abs(dY1(:))) max(abs(dY2(:))) max(abs(dY3(:)))]);
    
    %ERR(iter) = chg;
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnC + lamda*norm(E,"fro"); 
            err = norm(dY1(:)); 
        end
    end
    if chg < tol
        break;
    end 
    
    Y1 = Y1 + mu*dY1;
    Y2 = Y2 + mu*dY2;
    Y3 = Y3 + mu*dY3;
    mu = min(rho*mu,max_mu);   
end
obj = tnnC+lamda*norm(E,"fro");  
err = norm(dY1(:));
end