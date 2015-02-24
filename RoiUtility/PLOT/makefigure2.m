function f = makefigure2(dat,varargin)
% f = makefigure2(subclusters(j).BARPLOT.dat,str,cnames,splitat,docenter,covs);
% This one does a line plot
warning off

    if length(varargin) > 2,splitat=varargin{3};,else,splitat = 0;,end
    if length(varargin) > 3,docenter=varargin{4};,else,docenter = 0;,end
    if length(varargin) > 4,covs = varargin{5}; covs(:,end+1) = 1;,else,covs = [];,end
    if length(varargin) > 5,covs2 = varargin{6}; ,end
    
    %if ~isempty(covs), 
        %disp('Removing covariates')
        %disp('Before'), mean(dat),std(dat)
        %covs(:,1:end-1) = scale(covs(:,1:end-1));
        
     %   for i = 1:size(dat,2),
     %       b = pinv(covs) * dat(:,i);
     %       r = dat(:,i) - covs * b + b(end);
     %       dat(:,i) = r;
     %   end
        %disp('After')
        %mean(dat),std(dat)
    end
    
    if docenter,
        disp('Centering rows of data for plotting.')
        dat = dat - repmat(nanmean(dat')',1,size(dat,2));
        
    end
    
    doprintt = 1;
    
    if doprintt
        % print t values
        for i = 1:size(dat,2)
            [h,p,ci,stats] = ttest(dat(:,i));
            nums{i} = sprintf('%3.2f',stats.tstat);
        end
    end

    ste = nanstd(dat) ./ sqrt(size(dat,1));
    dat1 = dat;
    dat = nanmean(dat);
    mysum = sign(dat) .* (ste + abs(dat) + .15 .* abs(dat));
    
    if splitat,
        f = figure('Color','w'); 
        
        h(1) = subplot(1,2,1); set(gca,'FontSize',16); hold on; grid on;
        bar(dat(1:splitat)); tor_bar_steplot(dat(1:splitat),ste(1:splitat),{'b'});
        set(gca,'XTick',1:splitat)
        xlabel('Conditions','FontSize',18),ylabel('fMRI Signal','FontSize',18)
        if length(varargin) > 1, set(gca,'XTickLabel',varargin{2}(1:splitat)),end
        if length(varargin) > 0, title(varargin{1},'FontSize',20),end
        
        for i = 1:splitat, text(i-.5,mysum(i),nums{i},'Color','k','FontWeight','b','FontSize',14);,end
        
        
        h(2) = subplot(1,2,2); set(gca,'FontSize',16); hold on; grid on;
        bar(dat(splitat+1:end)); tor_bar_steplot(dat(splitat+1:end),ste(splitat+1:end),{'b'});
        set(gca,'XTick',1:length(dat)-splitat)
        xlabel('Conditions','FontSize',18)
        if length(varargin) > 1, set(gca,'XTickLabel',varargin{2}(splitat+1:end)),end
        
        for i = 1:length(dat)-splitat, text(i,mysum(splitat+i),nums{splitat+i},'Color','k','FontWeight','b','FontSize',14);,end
        
        equalize_axes(h)
        set(gcf,'Position',[195         308        1259         674]), drawnow
        
    else
        
        f = figure('Color','w'); set(gca,'FontSize',16); hold on; % grid on;
        %bar(dat);
        %tor_bar_steplot(dat,ste,{'b'});
    
        % special insert for intext
        % within s.e. {'DE' 'OE' 'AE' 'NE' 'DI' 'OI' 'AI' 'NI'}
        %for n = 1:size(dat1,1)  % subtract out no switch as baseline
        %    dat1(n,1:4) =dat1(n,1:4) - dat1(n,4);
        %    dat1(n,5:8) =dat1(n,5:8) - dat1(n,8);
        %end
        %c = [1 1 -1 -1 1 1 -1 -1; 1 -1 1 -1 1 -1 1 -1; 1 -1 -1 1 1 -1 -1 1];
        %
        
c = [-1 -1 -1 -1 1 1 1 1; ...
        1 1 -1 -1 1 1 -1 -1; ...
        1 -1 1 -1 1 -1 1 -1; ...
        1 -1 -1 1 1 -1 -1 1; ...
        0 1 -1 0 0 1 -1 0; ...
        0 1 0 -1 0 1 0 -1; ...
        0 0 1 -1 0 0 1 -1; ...
        1 -1 -1 1 0 0 0 0; ...
        0 0 0 0 1 -1 -1 1; ...
        0 1 0 -1 0 0 0 0; ...
        0 0 0 0 0 1 0 -1; ...
        0 0 1 -1 0 0 0 0; ...
        0 0 0 0 0 0 1 -1; ...
        0 1 1 -2 0 0 0 0; ...
        0 0 0 0 0 1 1 -2
];

c(1:4,:) = c(1:4,:) ./ 4; % sum to 1 so scale is the same as means on plot
c(5:9,:) = c(5:9,:) ./ 2; % sum to 1 so scale is the same as means on plot
 

cov2 = covs(:,1:end-1); cov2 = cov2 - repmat(mean(cov2),size(cov2,1),1);
[str,res]=repeated_anovan(dat1,cov2,c,1);
%res.flag = ones(1,size(dat1,1));

str2 = [sprintf(str)]; disp(str2)
       %'IEme\tOSme\tASme\tINT\tOS>AS\tOSonly\tASonly\tIntE\tIntI\tOSonlE\tOSonlI\tASonlyE\tASonlyI:\t ')
       
        tmp = (c * dat1(find(res.flag),:)')'; 
        ste = mean(std(tmp) ./ sqrt(size(tmp,1)));
        tmps = repmat(ste,2,4);
        
        dat = nanmean(dat1(find(res.flag),:));  % re-compute mean w/o outliers
 
        % r-square from behavioral regs
        %covs2(:,end+1) = 1;
        %for i = 10:13
        %    [B,BINT,R,RINT,STATS] = regress(tmp(:,i),covs2(find(res.flag),:));
        %    z3(i-9) = STATS(2);
        %end
        
        % subtract no switch condition as baseline
        dat = dat - mean(dat([4 8]));
        %dat(1:4) = dat(1:4) - dat(4); dat(5:8) = dat(5:8) - dat(8);
        
        tmp = reshape(dat,2,4); tmp = flipud(tmp);  % flip to get no switch first

        mk = {'bo-' 'rd--' 'bs-' 'r^--'};
        hh(1)= subplot(1,3,1); set(gca,'FontSize',14)
        for lin = 1:2    %size(tmp,2), 
            h(lin) = plot(tmp(:,lin),mk{lin},'LineWidth',2);,
            tor_bar_steplot(tmp(:,lin)',tmps(:,lin)',mk(lin));
        end

        xlabel('Attribute Switch'),set(gca,'XLim',[.5 2.5],'XTick',[1 2],'XTickLabel',{'No Sw.' 'Switch'})
        ylabel('BOLD Regression weights')
        %tmp3 = get(gca,'YLim'); set(gca,'YLim',[-.05 tmp3(2)]);
        title('External')
        
        hh(2) = subplot(1,3,2); set(gca,'FontSize',14)
        for lin = 3:4    %size(tmp,2), 
            h(lin) = plot(tmp(:,lin),mk{lin},'LineWidth',2);,
            tor_bar_steplot(tmp(:,lin)',tmps(:,lin)',mk(lin));
        end
        %legend(h,{'Obj. SW (Ext.)' 'Obj. NS (Ext.)' 'Obj. SW (Int.)' 'Obj. NS (Int.)'},2);
        xlabel('Attribute Switch'),set(gca,'XLim',[.5 2.5],'XTick',[1 2],'XTickLabel',{'No Sw.' 'Switch'})
        equalize_axes(hh);
        title('Internal')
        
       subplot(1,3,1);legend(h,{'Object Switch' 'No Object Switch'},0);
                
        set(gcf,'Position',[85         532        1482         524])
        
        %title(str,'FontSize',10)
        drawnow
    
        % second plot!
        
        % EO ob IO intern IA att EA extern
        ord = [10 6 11 15 13 7 12 14];
        
        subplot(1,3,3)
        %z = res.t(1,10:13); z(z < 0) = .2;
        z = res.t(1,ord); z(z < 0) = .2;
        % obj correl with obj cue, att correl with att cue
        %z2 = [res.t(3,[10 11]) res.t(2,[12 13])]; z2(z2 < 0) = .2;
        z2 = [res.t(3,[10 6 11]) res.t(2,[15 13 7 12 14])]; z3=-z2; z2(z2 < 0) = .2;
        z3(z3 < 0) = .2;
        
        %x = [0 1 0 -1]; y = [1 0 -1 0];
        x = [0 .5 1 .5 0 -.5 -1 -.5]; y = [1 .5 0 -.5 -1 -.5 0 .5];
        
        fill(x.*z,y.*z,'r'); hold on; fill(x.*z2,y.*z2,'b');
        
        % extra F value
        %fill(x.*z3,y.*z3,'g');
        
        alpha(.5)
        hold on; zm = max([4.5 z z2]); plot([-zm zm],[0 0],'k','LineWidth',2)
        hold on; plot([0 0],[-zm zm],'k','LineWidth',2)

        set(gca,'FontSize',14);xlabel('t-value');ylabel('t-value');
        set(gca,'XTick',-4:4,'YTick',-4:4,'XTickLabel',[4 3 2 1 0 1 2 3 4],'YTickLabel',[4 3 2 1 0 1 2 3 4])
        text(.1,zm,'Ext. Object')
        text(.5*zm,.2,'Int. Object')
        text(.1,-zm,'Int. Attribute')
        text(-zm,.2,'Ext. Attribute')
        set(gca,'XLim',[-4.5 4.5],'YLim',[-4.5 4.5])
        lim = tinv(.95,res.dfe);
        %plot([-lim lim],[lim lim],'k');plot([-lim lim],[-lim -lim],'k')
        %plot([-lim -lim],[-lim lim],'k');plot([lim lim],[-lim lim],'k')
        
        [x,y]=ellipse((randn(1000,1)),lim,lim,1); 
        hh = plot(x(2:end-1),y(2:end-1),'k-');
        
        drawnow
    end  
    
    return