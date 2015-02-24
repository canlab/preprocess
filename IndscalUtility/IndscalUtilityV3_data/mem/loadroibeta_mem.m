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
roimat={'vox_lerc3','vox_rerc3','vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa4','vox_lppa4','vox_hippo','vox_rdlpfc','vox_lpar','vox_lvlpfc'}; %46

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
        disp([thepath,'mem',num2str(s),'\sntaskTD\beta',r])
        load([thepath,'mem',num2str(s),'\sntaskTD\beta',r]);
        cb=reshape(cbeta,2,3,length(cbeta)/6);
        cb1=squeeze(cb(1,1:2,:));       %take only first beta (canonical HRF)
        cb2=mean(cb1(1,:)-cb1(2,:));    %for hits-misses
        meanbeta(ss,rr)=cb2;
    end
end

