function timeinfo=dotime(xl,names,clustermat,numperm)

useflip=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%       BETWEEN CLUSTER            %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['BETWEEN CLUSTER STATS...']);
mlat = mean(xl,3);          %between points anova for time values
tmp = abs(mlat);
[tmp2,e] = cmdscale(tmp);   % e gives eigenvalues for y y' - should be positive and zeros for euclidean
%thefigure;plot(tmp2(:,1),tmp2(:,2),'o','MarkerFaceColor','k')

tor_fig;
subplot(2,1,1);
colors = {'r' 'g' 'b' 'y' 'c' 'm' 'k'};
while max(clustermat) > length(colors), colors = [colors colors];,end
for j=1:size(tmp2,1);
hold on;plot(tmp2(j,1),0,'ok','MarkerFaceColor',char(colors{clustermat(j)}))
end

for i = 1:length(names),text(tmp2(i),.1,names{i}),end
set(gca,'FontSize',16,'YTickLabel',[]); title(' CLUSTER Timeline','FontSize',18),xlabel('Time'),
%saveas(gcf,'cluster_timeline','fig')

tmp3=tmp2(:,1);
cell{1}=clustermat;
[p,t]=anovan(tmp3,cell(1),[],[],[],'off'); 
timeinfo.meanp=p;
timeinfo.meant=t{2,6};

subplot(2,1,2);
for c=1:max(clustermat);
        meanc=mean(tmp3(find(clustermat==c)));
        hold on;plot(meanc,c,'ok','MarkerFaceColor',char(colors{c}),'MarkerSize',10);         %plot unflipped
        line([max(tmp3(find(clustermat==c))) min(tmp3(find(clustermat==c)))],[c c],'color',char(colors{c}),'linewidth',3);
end 
YLim([0 max(clustermat)+1]);

%%FIND MIN AND MAX AS GUIDE TO FLIP INDIVIDUALS
for c=1:max(clustermat);
meanc(c)=mean(tmp3(find(clustermat==c)));
end
if useflip==1;                      %choose left-right across all clusters
clustmin=find(meanc==min(meanc));
clustmax=find(meanc==max(meanc));
elseif useflip==2;                  %assume min and max from 2D mds.
    if meanc(1)<meanc(max(clustermat));
        clustmin=1;
        clustmax=max(clustermat);
    else
        clustmax=1;
        clustmin=max(clustermat);
    end
end
        
disp(['omnibus p-value for between cluster stats...',num2str(p)]);

%DO POST-HOC TESTING

numclust=max(clustermat);
pairs=possibly(numclust);
for test=1:size(pairs,1);
    newcm=clustermat;
    newtmp3=tmp3;
    newcm(find(clustermat~=pairs(test,1) & clustermat~=pairs(test,2)))=[];
    newtmp3(find(clustermat~=pairs(test,1) & clustermat~=pairs(test,2)))=[];
    ncm{1}=newcm;
    [p,t]=anovan(newtmp3,ncm(1),[],[],[],'off');
    close;
    disp(['p-value for comparison of cluster ',num2str(pairs(test,1)),' and cluster ',num2str(pairs(test,2)),'...',num2str(p)]);
    clear ncm;clear newcm;clear newtmp3;        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%       BETWEEN SUBJECT            %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['BETWEEN SUBJECT STATS...']);
disp(['flipping on the basis of clusters ',num2str(clustmin),' and ',num2str(clustmax)]);

for n=1:size(xl,3)                  %calculate between subjects anova for time values 
    smlat=squeeze(xl(:,:,n));  
    tmp = abs(smlat);
    [tmp2,e] = cmdscale(tmp); 
    timeline=squeeze(tmp2(:,1))';
    for c=1:max(clustermat);           %find mean cluster values from unflipped nmds
        meanc(c)=mean(timeline(find(clustermat==c)));
    end

    if meanc(clustmax)<meanc(clustmin);         %if needs flipped
        timeline=timeline*-1;              %flip timeline        
        meanc=meanc*-1;                    %flip mean values for stats
    end
    meanlat(n,:)=meanc;                         %generate meanlat
    tmp5(n,:)=timeline;    
end

for n=1:max(clustermat);grup(n,:)=ones(1,size(xl,3))*n;end
grup=grup';
grup=grup(:);
meanlat=meanlat(:);

[p,t]=kruskalwallis(meanlat,grup,'off');                      %do stats on flipped values
%[p,t]=anova2(meanlat,1,'off');                      %do stats on flipped values
timeinfo.subp=p(1);
timeinfo.subt=t{2,5};

for np=1:numperm                    %%%PERMUTE    
    clear meanlat;
    clear grup;
    
    if np/(numperm/10)==fix(np/(numperm/10));
        disp(['permutation...',num2str(np)]);
    end
    
    for n=1:size(xl,3)                  %calculate between subjects anova for time values 
    smlat=squeeze(xl(:,:,n));  
    tmp = abs(smlat);
    [tmp2,e] = cmdscale(tmp); 
    timeline=shuffle(squeeze(tmp2(:,1))');
    for c=1:max(clustermat);           %find mean cluster values from unflipped nmds
        meanc(c)=mean(timeline(find(clustermat==c)));
    end

    if meanc(clustmax)<meanc(clustmin);         %if needs flipped
        timeline=timeline*-1;              %flip timeline        
        meanc=meanc*-1;                    %flip mean values for stats
    end
    meanlat(n,:)=meanc;                         %generate meanlat
    tmp5(n,:)=timeline;    
end

meanlat=meanlat(:);
for n=1:max(clustermat);grup(n,:)=ones(1,size(xl,3))*n;end
grup=grup';grup=grup(:);

[p,t]=kruskalwallis(meanlat,grup,'off');    


%[p,t]=anova2(meanlat,1,'off');                      %do stats on flipped values
timeinfo.permsubp(np)=p(1);
timeinfo.permsubt(np)=t{2,5};

end             %ENDING PERMUTE

%sorted=sort(timeinfo.permsubt);
pval=1-(sum(timeinfo.subt>timeinfo.permsubt)./numperm);
disp(['p-value for between subject stats...',num2str(pval)]);


thefigure;                          %PLOT
for n=1:size(xl,3)                  %calculate between subjects anova for time values 
    timeline=squeeze(tmp5(n,:));
    for c=1:max(clustermat);
        meanc=mean(timeline(find(clustermat==c)));
        hold on;plot(meanc,(n*2)+(c/3),'ok','MarkerFaceColor',char(colors{c}),'MarkerSize',10);         %plot unflipped
        line([max(timeline(find(clustermat==c))) min(timeline(find(clustermat==c)))],[(n*2)+(c/3) (n*2)+(c/3)],'color',char(colors{c}),'linewidth',3);
    end 
end
   

        figure; 
        histvals=hist(timeinfo.permsubt);
        hist(timeinfo.permsubt);line([timeinfo.subt timeinfo.subt],[0 max(histvals)],'color','r');
        title(['permuted latency t-values']);


   
keyboard

