function Collapse_Trial(trialtype, extent, cleardat)
    plotdata = get(gcf, 'UserData');
    
    scrsz = get(0,'ScreenSize');
    
    
    if nargin<3
        cleardat = 0;
    end;
    if nargin<2
        extent = 'global';
    end;
    
    tfig = findobj('Name', 'Trial Data Display');
    if isempty(tfig)
        tfig = figure('Position',[scrsz(3)/3 1 scrsz(3)*2/3 scrsz(4)/2], 'UserData', trialtype, 'Name', 'Trial Data Display');
        ax = axes;
        set(ax, 'UserData', 'tplot')
    else
        if ~cleardat
            trialtype = [trialtype; get(tfig, 'UserData')];
        end;
        set(tfig, 'UserData', trialtype);
    end;
    
    x=[];
    y=[];
    
    for tcnt = 1:length(trialtype)
    
    switch extent
       case 'session'
        subj = plotdata.currentSub;
        sess = plotdata.currentSession;
        fs = plotdata.fs{subj}{sess};
        [newx, newy] = collectPlotdata(trialtype(tcnt), subj, sess, fs, plotdata);
       case 'subj'
        subj = plotdata.currentSub;
        for i = 1:length(plotdata.fs{subj});
            fs = plotdata.fs{subj}{sess};
            [newx, newy] = collectPlotdata(trialtype(tcnt), subj, i, fs, plotdata);
            [x,y] = resolvexy(x,y,newx,newy);
        end;
       case 'global'
        for j = 1:length(plotdata.fs)
            for i = 1:length(plotdata.fs{j});
                fs = plotdata.fs{j}{i};
                [newx, newy] = collectPlotdata(trialtype(tcnt), j, i, fs, plotdata);
                [x,y] = resolvexy(x,y,newx,newy);
            end;
        end;
    end;
    
    
    
    end;
    
    ax = findobj(tfig, 'UserData', 'tplot');
    
    plot(ax,x,y);
    set(ax, 'UserData', 'tplot')
    
end



function [x,y] = collectPlotdata(tt, subj, sess, fs, plotdata) 
    times = plotdata.trials{subj}{sess}(find(plotdata.trials{subj}{sess}(:,3) == tt),1:2);
    n = length(times(:,1));
    if ~isempty(times)
        maxx = max(times(:,2)-times(:,1))+1;
        x = zeros(maxx,n);
        y = zeros(maxx,n);
        
        for i = 1:n
            temp = plotdata.signal{subj}{sess}(times(i,1):times(i,2));
            if length(temp) < maxx
                temp = [temp; zeros(maxx-length(temp),1)];
            end;
            y(:,i) = temp;
            x(:,i) = (1:maxx)/fs;
            
        end;
    else
        x=[];
        y=[];
    end;
end
    
    
function [x,y] = resolvexy(oldx,oldy,newx,newy)

    if ~isempty(oldx) && ~isempty(newx)
        olds = size(oldx);
        news = size(newx);
        
        if olds(1) < news(1)
            oldx = [zeros(news(1)-olds(1),olds(2)); oldx];
            oldy = [zeros(news(1)-olds(1),olds(2)); oldy];
        elseif news(1) < olds(1)
            
            newx = [zeros(olds(1)-news(1),news(2)); newx];
            newy = [zeros(olds(1)-news(1),news(2)); newy];
        end;
        
 
    end;
    
    x = [oldx newx];
    y = [oldy newy];
end    






