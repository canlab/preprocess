clear; close all
thepath='\\delta\chris\fhmem\';
%roimat={'vox_lerc3','vox_rerc3','vox_lpara2','vox_rpara2','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa','vox_lppa','vox_rerc2','vox_hippo','vox_rdlpfc','vox_lpar'}; %31
%roimat={'vox_lerc3','vox_rerc3','vox_lpara2','vox_rpara2','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa','vox_lppa','vox_rerc2','vox_hippo'}; %32
%roimat={'vox_lerc3','vox_rerc3','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa','vox_lppa','vox_hippo','vox_rerc2'}; %34
%roimat={'voxrffa','vox_lffa','vox_rppa','vox_lppa'}; %99
%roimat={'vox_lerc3','vox_rerc3','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa','vox_lppa','vox_hippo','vox_rdlpfc','vox_lpar'}; %35
%roimat={'vox_lerc3','vox_rerc3','vox_lvis','','vox_rloc','voxrffa','vox_lffa','vox_rppa','vox_lppa','vox_hippo'}; %33
% forties
%roimat={'vox_lerc4','vox_rerc4','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa4','vox_lppa4','vox_hippo'}; %43
%roimat={'vox_lerc3','vox_rerc3','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa4','vox_lppa4','vox_hippo'}; %44
%roimat={'vox_lerc3','vox_rerc3','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa4','vox_lppa4','vox_hippo','vox_rdlpfc','vox_lpar'}; %45
roimat={'vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa4','vox_lppa4','vox_lerc3','vox_rerc3','vox_hippo','vox_lpar','vox_rdlpfc','vox_lvlpfc'}; %46

submat=[2 3 4 6 8 9 11:20];

% get onsets
load('\\delta\chris\fhmem\onsetmat');
start=onsetmat([1:6 8:17],:,1);  %leave out the 7th subject (mem10)
%%%%%%%%%%

%%% load all roi data into a 1*subs cell with each entry:
%%% regions*blocks*dp

for ss=1:length(submat);
     s=submat(ss);
    for rr=1:length(roimat);
        r=roimat{rr};
        clear ts;clear sdat;clear in;
        disp([thepath,'mem',num2str(s),'\sntaskTD\ts_',r])
        load([thepath,'mem',num2str(s),'\sntaskTD\ts_',r]);
        substart=start(ss,:);
        for k=1:sum(substart~=0);
        in=ts(round(substart(k)+1:substart(k)+40));
        sdat(k,:)=trimts(in,4);        
         end
        dat{ss,rr}=sdat;
      end
    trigs{ss}=squeeze(onsetmat(ss,1:sum(substart~=0),:));    
end

for n=1:length(trigs);t=trigs{n};t(1:3,:)=[];trigs{n}=t;end
for s=1:size(dat,1);for r=1:size(dat,2);d=dat{s,r};d(1:3,:)=[];dat{s,r}=d;end;end
save('C:\Program Files\MATLAB704\work\indscal\new\mem17\mem_dat47','dat');
save('C:\Program Files\MATLAB704\work\indscal\new\mem17\mem_trigs','trigs');
names=roimat;
save('C:\Program Files\MATLAB704\work\indscal\new\mem17\mem_names47','names');
