clc;
clear all;
load 'Yale32_0.1_7.15_10.mat'
load 'Yale32_V_else_7.15.mat'
load 'Yale32_V_Vexl_7.15.mat'

for k=7:1:15
    num=1;
    fprintf('#of classes: %g\n',k);
    for nter=1:10
        fprintf('#of iterations: %g\n',nter);
        %% new proposed 10
        TSNMF_para.maxIter=500;
        TSNMF_para.iter=2;
        TSNMF_para.k=k;
        TSNMF_para.m=10;
        TSNMF_para.lamda=1;
        TSNMF_para.beta=0.1;
        [ACC,NMI]=TSNMF(datasub{k},W{k},TSNMF_para,labelsub{k},A{k}{nter},V{k}{nter},Vexl{k}{nter});
        AC_weight(num,k-6) = ACC;
        MIhat_weight(num,k-6) = NMI;
        
        %% PCPSNMF
        PCPS_mu=1;
        PCPS_alpha=10;%optimal
        
        [Vn_PCPS,K_PCPS] = PCPSNMF(W{k},PCPS_mu,k,Z{k}{nter},PCPS_alpha,V_else{k}{nter});
        [tmp_PCPS label_PCPS] = max(Vn_PCPS, [], 2); 
        label_PCPS = bestMap(labelsub{k},label_PCPS);
            
        AC_PCPS(num,k-6) = length(find(labelsub{k} == label_PCPS))/length(labelsub{k});
        MIhat_PCPS(num,k-6) = MutualInfo(labelsub{k},label_PCPS); 
        
        %% S3NMF
        S3NMF_para.lamda1=1;
        S3NMF_para.lamda2=1;
        S3NMF_para.mu=100;
        S3NMF_para.alpha=10;
        S3NMF_para.k=k;
        S3NMF_para.maxIter=500;
        [V_S3NMF] = S3NMF(W{k},Z{k}{nter},S3NMF_para,V_else{k}{nter});
        [tmp_S3NMF,label_S3NMF] = max(V_S3NMF, [], 2); 
        label_S3NMF = bestMap(labelsub{k},label_S3NMF);
            
        AC_S3NMF(num,k-6) = length(find(labelsub{k} == label_S3NMF))/length(labelsub{k});
        MIhat_S3NMF(num,k-6) = MutualInfo(labelsub{k},label_S3NMF); 


        %% MVCHSS
        MVCHSS_para.p=5;
        MVCHSS_para.k=k;
        MVCHSS_para.subN=subN{k};
        MVCHSS_para.mu=0.1;
        MVCHSS_para.alpha=10;
        MVCHSS_para.beta=10;
        MVCHSS_para.randtimes=1;  % randtimes = selected from 1 to 5
        MVCHSS_para.maxIter=100;
        [V_MVCHSS, obj_record] = MVCHSS(Trans_datasub{k},A{k}{nter},MVCHSS_para,V_else{k}{nter},Vexl{k}{nter});
        [tmp_MVCHSS,label_MVCHSS] = max(V_MVCHSS, [], 2); 
        label_MVCHSS = bestMap(labelsub{k},label_MVCHSS);
            
        AC_MVCHSS(num,k-6) = length(find(labelsub{k} == label_MVCHSS))/length(labelsub{k});
        MIhat_MVCHSS(num,k-6) = MutualInfo(labelsub{k},label_MVCHSS); 

        
        %% SANMF
        P_SANMF=A{k}{nter};
        P_SANMF(A{k}{nter}<0)=1;
        Z2=P_SANMF-Z{k}{nter};
        SANMF_para.gamma=10;
        SANMF_para.c=k;
        SANMF_para.maxiter=500;
        SANMF_para.alpha=0.1;
        SANMF_para.beta=1;
        SANMF_para.eta=10;
        [V_SANMF]=SANMF(W{k},P_SANMF,Z{k}{nter},Z2,SANMF_para,V_else{k}{nter});
        [tmp_SANMF,label_SANMF]=max(V_SANMF,[],2);
        label_SANMF = bestMap(labelsub{k},label_SANMF);
        AC_SANMF(num,k-6) = length(find(labelsub{k} == label_SANMF))/length(labelsub{k});
        MIhat_SANMF(num,k-6) = MutualInfo(labelsub{k},label_SANMF);


        %% SNMF
        SNMF_para.maxIter=500;
        SNMF_para.k=k;
        [V_SNMF]=SNMF(W{k},SNMF_para,V_else{k}{nter});
        [tmp1_SNMF,label_SNMF]=max(V_SNMF,[],2);
        label_SNMF=bestMap(labelsub{k},label_SNMF);
        AC_SNMF(num,k-6) = length(find(labelsub{k} == label_SNMF))/length(labelsub{k});
        MIhat_SNMF(num,k-6) = MutualInfo(labelsub{k},label_SNMF);

        %% NMF
        [~,V_NMF] = nmf(datasub{k},k);
        V_NMF=V_NMF';
        label_NMF = kmeans(V_NMF, k);
        label_NMF=bestMap(labelsub{k},label_NMF);
        AC_NMF(num,k-6) = length(find(labelsub{k} == label_NMF))/length(labelsub{k});
        MIhat_NMF(num,k-6) = MutualInfo(labelsub{k},label_NMF);
        
        %% GNMF
        GNMF_para = [];
        GNMF_para.maxIter = 500;
        GNMF_para.alpha = 100;
        GNMF_para.nRepeat = 5;
        [~,gH] = GNMF(datasub{k},k,W{k},GNMF_para,V_else{k}{nter});

        label_GNMF = litekmeans(gH,k,'Replicates',20);
        label_GNMF = bestMap(labelsub{k},label_GNMF);        

        AC_GNMF(num,k-6) = length(find(labelsub{k} == label_GNMF))/length(labelsub{k});
        MIhat_GNMF(num,k-6) = MutualInfo(labelsub{k},label_GNMF);
        
        %% GSNMF
        GSNMF_para = [];
        GSNMF_para.maxIter = 500;
        GSNMF_para.alpha = 100;
        GSNMF_para.nRepeat = 5;
        GSNMF_para.triF = 0; %bi factorization
        
        [~,V_GSNMF] = GNMF_S(graphL{k},k,Z{k}{nter},GSNMF_para,V_else{k}{nter});
        [~, label_GSNMF] = max(V_GSNMF, [], 2); 
        label_GSNMF = bestMap(labelsub{k},label_GSNMF);
        AC_GSNMF(num,k-6) = length(find(labelsub{k} == label_GSNMF))/length(labelsub{k});
        MIhat_GSNMF(num,k-6) = MutualInfo(labelsub{k},label_GSNMF);
        
        %% SNMFCC
        Q=A{k}{nter};
        Q(A{k}{nter}==1)=0.1;
        Q(A{k}{nter}==-1)=3;

        Q=Q-diag(diag(Q));
        
        V_SNMFCC = SNMFCC(graphL{k},Q,Z{k}{nter},k,V_else{k}{nter});
        [~, label_SNMFCC] = max(V_SNMFCC, [], 2); 
        label_SNMFCC = bestMap(labelsub{k},label_SNMFCC);
        AC_SNMFCC(num,k-6) = length(find(labelsub{k} == label_SNMFCC))/length(labelsub{k});
        MIhat_SNMFCC(num,k-6) = MutualInfo(labelsub{k},label_SNMFCC); 


        
        num=num+1;
    end
end
%% compute result
%ACC
ans = mean(AC_SNMF,1);
bns = mean(AC_PCPS,1);
cns = mean(AC_SANMF,1);
dns = mean(AC_weight,1);
ens = mean(AC_S3NMF,1);
sns = mean(AC_MVCHSS,1);
kns = mean(AC_NMF,1);
lns = mean(AC_GNMF,1);
mns = mean(AC_GSNMF,1);
nns = mean(AC_SNMFCC,1);
%NMI
fns = mean(MIhat_SNMF,1);
gns = mean(MIhat_PCPS,1);
hns = mean(MIhat_SANMF,1);
ins = mean(MIhat_weight,1);
jns = mean(MIhat_S3NMF,1);
tns = mean(MIhat_MVCHSS,1);
ons = mean(MIhat_NMF,1);
pns = mean(MIhat_GNMF,1);
qns = mean(MIhat_GSNMF,1);
rns = mean(MIhat_SNMFCC,1);