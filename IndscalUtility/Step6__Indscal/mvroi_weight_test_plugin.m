% only really finished for 2 task states
% could add statistics also

comp = DATA.SPEC.comps(1,:)';
complong = repmat(comp,size(W,1)./length(comp),1);

W = DATA.INDSCAL.W;

% get a vector of task states (integers) corresponding to W
% store in taskvector
% this assumes that task states are nested within subjects in W
% (same format as DATA.DATA.xc)

tv = zeros(size(W,1),1);

for i = 1:length(comp)
    wh = find(complong == comp(i));
    tv(wh) = i;
end
    
DATA.INDSCAL.taskvector = tv;

numstates = max(tv);

figure;
plot(W(tv==1,1),W(tv==1,2),'ro','MarkerFaceColor','r');hold on;
if numstates > 1
    plot(W(tv==2,1),W(tv==2,2),'bo','MarkerFaceColor','b');
end

try,legend(DATA.SPEC.tasknames);,catch, disp('Cannot find DATA.tasknames'), end

if numstates > 1
    disp('Step 6: Drawing arrows between successive states (arrows assume 2 states.)');
    for i = 1:numstates:size(W,1)-1
        arrow([W(i,1) W(i,2)],[W(i+1,1) W(i+1,2)],'Length',10,'Tipangle',20)
    end
end

xlabel('Dimension 1'); ylabel('Dimension 2');
title('Subject weights in 1st 2 dimensions by task state');

drawnow

try, saveas(gcf,'step6_indscal_weights','fig');,saveas(gcf,'step6_indscal_weights','tif'); close,     catch, disp('Error saving figure.'), end




if numstates > 1
    
    
% difference between weights for task states, for MANOVA
wdiff = W(tv==1,:) - W(tv==2,:);

% now stats

Wdif = contrast2d(DATA.INDSCAL.W,DATA.SPEC.comps(1,:)');
Wavg = contrast2d(DATA.INDSCAL.W,ones(size(DATA.SPEC.comps(1,:)))'./length(DATA.SPEC.comps(1,:)));


[h,p,ci,stat] = ttest(Wdif);

fprintf(1,['Differences between conditions in indscal weights, ' DATA.SPEC.comptitle '\n']);
fprintf(1,'Mean diff\t');
fprintf(1,'%3.2f\t',mean(Wdif));
fprintf(1,'\n');

fprintf(1,'t\t');
fprintf(1,'%3.2f\t',stat.tstat);
fprintf(1,'\n');

fprintf(1,'p-value\t');
fprintf(1,'%3.4f\t',p);
fprintf(1,'\n');

if isfield(DATA.SPEC,'beh')
    if ~isempty(DATA.SPEC.beh)
        fprintf(1,'\n');
        
        x = [Wavg Wdif];
        for i = 1:size(Wavg,2), nms{i} = ['Dim ' num2str(i) ' Average'];,end
        for i = 1:size(Wdif,2), nms{end+1} = ['Dim ' num2str(i) ' Contrast'];,end
        
        if length(DATA.SPEC.beh) ~= size(x,1), error('Behavior vector does not match num. of subjects.'),end
    
        fprintf(1,'\n');
        fprintf(1,['Behavioral prediction of indscal weights\n']);
        fprintf(1,['Avg weights then contrast across weights\n']);

        DATA.INDSCAL.STEPWISE = stepwise_tor(x,DATA.SPEC.beh,nms); 

        %[b,se,pval,inmodel,stats] = stepwisefit(x,DATA.SPEC.beh,'penter',.10,'display','on');
        
    end
end

DATA.INDSCAL.Wdif = Wdif;
DATA.INDSCAL.Wavg = Wavg;

else
    % only 1 state
    
    
    if isfield(DATA.SPEC,'beh')
    if ~isempty(DATA.SPEC.beh)
        fprintf(1,'\n');

        for i = 1:size(W,2), nms{i} = ['Dim ' num2str(i)];,end
        
        if length(DATA.SPEC.beh) ~= size(W,1), error('Behavior vector does not match num. of subjects.'),end
    
        fprintf(1,'\n');
        fprintf(1,['Behavioral prediction of indscal weights\n']);

        DATA.INDSCAL.STEPWISE = stepwise_tor(W,DATA.SPEC.beh,nms); 

        %[b,se,pval,inmodel,stats] = stepwisefit(x,DATA.SPEC.beh,'penter',.10,'display','on');
        
    end
end
    
    
    
    
    
    
    
    
end


