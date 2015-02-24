function [dat,DX]=getvismotdata(ncond,eres);

mypath='c:\matlab6p5\work\indscal\vismotor_data\';
submat={'10lb' '10rh' '10td'  '28kh' '28ms' '28ns' '28rm' '11dc' '11pd' '11ps' };
names ={ 'Rvis' 'Lvis' 'Rmot' 'Lmot'};

disp('creating design matrix');
load([mypath,'vnl_condf']);
events=makebinary(vnl_condf);
for n=1:ncond;sf{n}=squeeze(events(n,:));end;
[DX sf1]=tor_make_deconv_mtx3(sf,eres,1);

for s=1:length(submat);
    for r=1:length(names);
        disp([mypath,'roi_0105',char(submat(s)),'_',char(names(r)),'.mat']);
        name=['roi_0105',char(submat(s)),'_',char(names(r))];
        eval(['load ' [mypath,name]]);
        eval(['ROI= ' name]);  
        data=trimts(ROI.adjustedy,3,[],1);
        dat(:,r,s)=data;                
    end
end


save vismotdat.mat dat  DX names