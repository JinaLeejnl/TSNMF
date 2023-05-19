function [datasub,subN,labelsub,Trans_datasub]=loadsub(fulldata,gnd,k)

tureLabel=gnd;
temp=[];
labelsub=[];
for i=1:k
    t=find(tureLabel==i);
    temp=[temp;t];
    labelsub=[labelsub;tureLabel(t)];
end
subN=size(labelsub,1);
datasub=fulldata(:,temp);
Trans_datasub=datasub';

end