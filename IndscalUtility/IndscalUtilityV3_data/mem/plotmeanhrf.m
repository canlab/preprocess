clear;%
close all;
submat = [2 3 4 6 8 9 11 12 13:20];

path = '\\delta\chris\fhmem\';

%Screen('CloseAll')
%clear screen

%regmat={'L_VLPFC','R_VLPFC','L_DLPFC','R_DLPFC','L_PUTA', 'R_PUTA',  'MF_antinf','MF_ant'};
%regmat = { 'L_FFA','R_FFA','L_postHIPP','R_postHIPP','L_antHIPP','R_antHIPP','L_HIPP'};
%regmat = {'lvis','rvis', 'lpara','rpara','lhippo','rhippo'};
%regmat = {'lvis5','rvis5','lppa5','rppa5','lerc5','rerc5','lffa5','lpar5','rdlpfc5','rloc5','hippo5','rpara25','lpara25','rerc25'};
%regmat = {'vox_lerc3', 'vox_rerc3','vox_lppa','vox_rppa','vox_lvis','vox_rvis'};
%regmat={'voxrffa','vox_lffa','vox_rppa4','vox_lppa4','vox_rerc4','vox_lerc4','vox_rvis','vox_lvis'};
regmat={'vox_lvis','vox_rvis','vox_rloc','voxrffa','vox_lffa','vox_rppa4','vox_lppa4','vox_lerc3','vox_rerc3','vox_hippo','vox_lpar','vox_rdlpfc','vox_lvlpfc'}; %46
%regmat = {'lhippo', 'rhippo'};
%regmat = {'lvis' 'lvis'};
%regmat = {'lvis' ,'lvis','lpara', 'rpara'};
sessavg = zeros(length(submat),length(regmat),3,16);

for n=1:  length(submat);
    disp(['subject....',num2str(submat(n))]);
numsess=4;if submat(n)==3 | submat(n)==20,numsess=3;,end;
    for m = 1:length(regmat);
        s = submat(n);
        r = regmat{m};
        a = load([path, 'mem', num2str(s), '\sntaskTD\hrf', r]);
      % load([path, 'mem', num2str(s), '\ntask\beta', r]);
       % betaavg(n,m,:)=squeeze(mean(reshape(cbeta,3,numsess)'));
        sessavg(n,m,:,:) = squeeze(mean(a.hrf.psth,1));
     end
end

%figure;
for r=1:length(regmat);
    %subplot(1,length(regmat),r);
    figure('color','w');
    steplot(squeeze(sessavg(:,r,:,:)),'l');
    title(regmat{r});
    legend('hit','miss','retrieval');legend boxoff;
    set(gca,'box','off')
    set(gca,'Xticklabel',0:4:32)
    Xlabel('time (s)','FontSize',16)
    Ylabel('fMRI signal','FontSize',16)
    %saveas(gcf,['c:\chris\fhmem\sntaskTD\ttest\',regmat{r},'fig.tif'],'tif');
end

%figure;
%for r=1:length(regmat);
%  %  subplot(1,length(regmat),r);
%  figure;
%  steplot(squeeze(betaavg(:,r,:)),'b');
%end

        