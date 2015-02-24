function f = makefigure(dat,varargin)
% f = makefigure(subclusters(j).BARPLOT.dat,str,cnames,splitat,docenter,covs);

    if length(varargin) > 2,splitat=varargin{3};,else,splitat = 0;,end
    if length(varargin) > 3,docenter=varargin{4};,else,docenter = 0;,end
    if length(varargin) > 4,covs = varargin{5}; covs(:,end+1) = 1;,else,covs = [];,end
    
    % ----------------------------------------------------  
    % > REMOVE covariates before getting means
    % ---------------------------------------------------- 

    if ~isempty(covs), 
        %disp('Removing covariates')
        %disp('Before'), mean(dat),std(dat)
        %covs(:,1:end-1) = scale(covs(:,1:end-1));
        for i = 1:size(dat,2),
            b = pinv(covs) * dat(:,i);
            r = dat(:,i) - covs * b + b(end);
            dat(:,i) = r;
        end
        %disp('After')
        %mean(dat),std(dat)
    end
    
    if docenter,
        disp('Centering rows of data for plotting.')
        dat = dat - repmat(nanmean(dat')',1,size(dat,2));
        
    end
    
    
    % ----------------------------------------------------  
    % > print t-values, get averages and ste
    % ---------------------------------------------------- 

    doprintt = 1;
    
    if doprintt
        % print t values
        for i = 1:size(dat,2)
            [h,p,ci,stats] = ttest(dat(:,i));
            nums{i} = sprintf('%3.2f',stats.tstat);
        end
    end

    ste = nanstd(dat) ./ sqrt(size(dat,1));
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
        
        
        % ----------------------------------------------------  
        % > Special INTEXT stuff - 2 x 3 grouped?
        % ---------------------------------------------------- 

        f = figure('Color','w'); set(gca,'FontSize',16); %hold on; grid on;
        %bar(dat);
        %tor_bar_steplot(dat,ste,{'b'});
    
        % special insert for intext
        cla
        xp = dat; h = bar([xp(1:3);xp(4:6)],'grouped');,cm = [0 1 0;1 0 0; 0 0 1];colormap(cm)
        xe = ste;
        tor_bar_steplot([xp(1:3)], xe(1:3),{'k'},.55,.225)
        tor_bar_steplot([xp(4:6)], xe(4:6),{'k'},1.55,.225)
        set(gca,'FontSize',16,'XTickLabel',{'External' 'Internal'});legend(h,{'Object switching' 'Attribute Switching' 'Interaction'})
        ylabel('BOLD Contrast')
        
        %set(gca,'XTick',1:size(dat,2))
        %xlabel('Conditions','FontSize',18),ylabel('fMRI Signal','FontSize',18)
        %if length(varargin) > 1, set(gca,'XTickLabel',varargin{2}),end
        if length(varargin) > 0, title(varargin{1},'FontSize',20),end

        set(gcf,'Position',[464   283   930   827]), drawnow
    
        
            
    end  
    
    return