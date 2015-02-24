clear region1;
regmat={'L VIS','R VIS','R LOC','R FFA','L FFA','R PPA','L PPA','L ERC','R ERC','L HIPPO','L PAR','R DLPFC','L VLPFC'}; %45 regions
load sessavg

for r=1:length(regmat);
    sessavg(:,r,:,:)=scaler(sessavg(:,r,:,:));
end

region=squeeze(sessavg(:,:,1,:)-sessavg(:,:,2,:));  %hits -misses

for reg=1:size(region,2);
    for bet=1:size(sig,2);
        for n=1:size(sig,1);
            dw=smoothn(squeeze(region(:,reg,n)),1);
            [r p]=corrcoef(dw,sig(:,bet));
            region1(reg,bet,n)=r(1,2);
            regionp(reg,bet,n)=p(1,2);
        end;
         [r p]=corrcoef(squeeze(mean(region(:,reg,8:14),3)),sig(:,bet));        
         reggie1(reg,bet)=r(1,2);
         reggiep(reg,bet)=p(1,2);
    end;
end
figure;
for r=1:size(region1,1);
    subplot(4,4,r);
    r1=squeeze(region1(r,:,:))';
%    p1=squeeze(regionp(r,:,:))'>0.05;
%    r2=r1.*p1;
    plot(r1,'linewidth',3);
%    plot(r2,'linewidth',3,'color','w');
    
    set(gca,'Ylim',[-1 1]);
    title(regmat{r});
end

%figure;
%for r=1:size(region,2);subplot(4,4,r);steplot(squeeze(region(:,r,:)));set(gca,'Ylim',[-0.3 0.3]);title(regmat{r});end

regs=[1 2 8 9];
rg=reggie1(regs,:);

[d d zrg]=r2z(rg(:),16);
zrg=reshape(zrg,size(rg,1),size(rg,2));

rp=reggiep(regs,:);
figure;bar(rg);
title('14-28s poststimulus');
legend('v-ffa','v-ppa','v-dlpfc');legend boxoff
set(gca,'Xticklabel',['lvis','rvis','lerc','rerc']);

for n=1:size(zrg,1);
    sd=sqrt(2/13);
    nzrg=zrg(n,:);
    a1vs3=nzrg(1)-nzrg(3)/sd
    a1vs3Z = a1vs3 / sd; 
    a1vs3p  = (1 - normcdf_t(abs(a1vs3Z))) .* 2
    a2vs3=(nzrg(2)-nzrg(3))/sd
    a2vs3Z = a2vs3 / sd; 
    a2vs3p  = (1 - normcdf_t(abs(a2vs3Z))) .* 2
end

keyboard

doplot=0;
if dopplot
rr1=region1(regs,:,1:7);
rr2=region1(regs,:,8:14);
[d rsig1 zregions1]=r2z(rr1(:),16);
zrr1=reshape(zregions1,size(rr1,1),size(rr1,2),size(rr1,3));    %regions * connections * dps
rrsig1=reshape(rsig1,size(rr1,1),size(rr1,2),size(rr1,3));
zrr11=mean(zrr1,3);
thefigure;subplot(2,1,1);
bar(zrr11);
title('0-14s poststimulus');
legend('v-ffa','v-ppa','v-dlpfc');legend boxoff
set(gca,'Xticklabel',['lvis','rvis','lerc','rerc']);
[d rsig2 zregions2]=r2z(rr2(:),16);
zrr2=reshape(zregions2,size(rr2,1),size(rr2,2),size(rr2,3));
rrsig2=reshape(rsig2,size(rr2,1),size(rr2,2),size(rr2,3));
zrr22=mean(zrr2,3);
subplot(2,1,2);
bar(zrr22);
title('14-28s poststimulus');
legend('v-ffa','v-ppa','v-dlpfc');legend boxoff
set(gca,'Xticklabel',['lvis','rvis','lerc','rerc']);
end

figure;
for n=1:length(regs);
subplot(1,length(regs),n);
plot(squeeze(region1(regs(n),:,:))','linewidth',3);
end


