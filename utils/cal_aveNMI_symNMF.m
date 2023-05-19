function [aveNMI,PairwiseNMI,DifNMI]=cal_aveNMI_symNMF(H)
for i=1:length(H)
    [~,res{i}]=max(H{i}');
end
% PairwiseNMI=zeros(length(H))*(zeros(length(H))-1)/2;
aveNMI=0;
PairwiseNMI=0;

n=length(H);
DifNMI=zeros(n,1);

for i=1:length(H)
    for j=1:length(H)
        try 
            %求两两之间的nmi
            temp= MutualInfo(res{i},res{j});
        catch
            temp=mean(PairwiseNMI);
        end
        %把两两之间的nmi存起来
        PairwiseNMI((i-1)*n+j)=temp;
        %两两之间nmi的加和
        aveNMI=temp+aveNMI;
        %存某个H和其他H之间的nmi和
        DifNMI(i)=DifNMI(i)+temp;
    end
end
aveNMI=aveNMI/(n*n);
