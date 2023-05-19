function [ACCssnmf,NMIssnmf,Cres]=cal_ACC_NMF_symNMF_v3(H,gnd)
[~,res]=max(H');
labelnew = res;
            %                 gndnew=gnd;
NMIssnmf= MutualInfo(gnd,labelnew);

labelnew = bestMap(gnd,labelnew);
ACCssnmf= length(find(gnd == labelnew))/length(gnd);
            
Pur=purity(max(gnd),labelnew,gnd);
ARI=RandIndex(labelnew,gnd);

[f,p,r] = compute_f(gnd,labelnew);
Cres.ACC=ACCssnmf;
Cres.NMI=NMIssnmf;
Cres.Pur=Pur;
Cres.ARI=ARI;
Cres.F1=f;
Cres.Pre=p;
Cres.Rec=r;
%% add the metrics for F1 score, Precision and Recall