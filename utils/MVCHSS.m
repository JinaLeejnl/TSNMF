function [V_final, obj_record] = MVCHSS(X,Z,opts,Vb,V_star)

k = opts.k;  % class number
subN=opts.subN ;  % sample number

view = 1;   % view number 

alpha=opts.alpha; % hyperparameter
beta=opts.beta;  % hyperparameter

% per=opts.per; % percentage of pairwiese constraints (0-1)

Times=opts.randtimes; % randtimes = selected from 1 to 5
maxIter=opts.maxIter; 

%%%% Pairwise Constraints Propagation on Hypergaph %%%%
[S, Dv, L] = MHPCP(X,Z,opts);


%%%% Learning Procedures %%%%
tryNo=0;
obj=1e+15;

while tryNo < Times 
    tryNo = tryNo+1;
    nIter=0;
    
    obj_temp_record = [];
    obj_temp = 0;
    for i=1:view
        obj_temp = obj_temp + sum(sum((S{i}-Vb*(Vb')).^2)) + ...
                   alpha*sum(sum((L{i}*Vb).*Vb)) + beta*sum(sum((V_star-Vb).^2));
    end
    obj_temp_record = [obj_temp_record, obj_temp];

    while nIter<maxIter
        
        % update V
        for i=1:view
            
            SV = (2+alpha)*S{i}*Vb + beta*V_star;
            VVV = 2*Vb*(Vb')*Vb + alpha*Dv{i}*Vb + beta*Vb;

            V = (SV./max(VVV,1e-10)) .* Vb;
            V = NormalizeK(V);
            Vb = V;

        end
        
        % update V_star
        sum_Vb=zeros(subN,k);
        for i=1:view
            sum_Vb = sum_Vb + Vb;
        end
        V_star = sum_Vb / view;
        
        nIter=nIter+1;      
    
        obj_temp = 0;
        for i=1:view
            obj_temp = obj_temp + sum(sum((S{i}-Vb*(Vb')).^2)) + ...
                       alpha*sum(sum((L{i}*Vb).*Vb)) + beta*sum(sum((V_star-Vb).^2));
        end
        obj_temp_record = [obj_temp_record, obj_temp];

    end
    
    FR = obj_temp_record(end);

    if FR < obj
        V_final = V_star;
        obj = FR;
        obj_record = obj_temp_record;
    end

end

function [NSelLoc,count] = SelLabSam_Semi_2(gnd, Per)
%  randomly	Select Per% samples from each class for semi-supervised NMF methods 
% but do not move the selected sample in the front
% where
% NSelLoc ... the location of all selected samples
% gnd ... the label information N*1
% Per ... the percent that select from each class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[N,M] = size(fea);
Nclass = length(unique(gnd));
NLabel = zeros(1,Nclass);
 
for i=1:Nclass
    LocCla = find(gnd==i);           % obtain the location for class i.
    NLabel(i) = round(Per*length(LocCla));
end
NSelLoc = zeros(1, sum(NLabel));    % obtain the location of all labeled sample
count = 0;
for i=1:Nclass
    LocCla = find(gnd==i);          % obtain the location for class i.
    SelLab = randperm(length(LocCla),NLabel(i));
    NSelLoc(count+1:count+NLabel(i)) = LocCla(SelLab);
    count = count + NLabel(i); 
end
if sum(NLabel) ~= length(NSelLoc)
    error('Error!');
end

function [K] = NormalizeK(K)
    % the dim of K is nSamp by mFeat
    % the function is used to normalize the row of a matrix
    K = K';
    n = size(K,2);
    norms = max(1e-15,sqrt(sum(K.^2,1)))';
    K = K*spdiags(norms.^-1,0,n,n);
    K = K';