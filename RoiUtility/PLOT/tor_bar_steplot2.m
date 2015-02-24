function tor_bar_steplot(avg,ste,plotcolor)
% function tor_bar_steplot(avg,ste,plotcolor)
%
% for grouped axes, with 2 bars per group
%

hold on
if isempty(plotcolor), plotcolor{1} = 'b';,end
wid = get(gca,'XLim');
wid = (wid(2) - wid(1))/50;
            
   for g = 1:size(ste,1)
       
            %ste = STE{g};
            
            % set x axis values
            % ---------------------------------------------------------------
            myx = 1:size(ste,1);
            
            %if isfield(Op,'window'),myx = Op.window(1):Op.window(2);,end
            
            if length(myx) ~= size(avg,1),
                warning('Window is wrong length, using default')
                myx = 1:size(avg,1);
            end
            
      			%for i = 1:size(ste,2)
                    
                    if length(plotcolor) < i, plotcolor{i} = plotcolor{1};,end
                    

                    % first one in group (left)
                    if avg(g,1) < 0, mymult = -1; else, mymult = 1;, end
                    j = myx(g) - .15;
         			try
                        plot([j j],[avg(g,1) avg(g,1) + mymult*ste(g,1)],plotcolor{1}(1))
         			    plot([j-wid j+wid],[avg(g,1) + mymult*ste(g,1) avg(g,1) + mymult*ste(g,1)],plotcolor{1}(1))
                    catch
                        disp('Can''t plot error bar.')
                        1
                        avg
                        ste
                        plotcolor
                    end
                    
                     % second one in group (right)
                     if avg(g,2) < 0, mymult = -1; else, mymult = 1;, end
                     j = myx(g) + .15;
                    try
                        plot([j j],[avg(g,2) avg(g,2) + mymult*ste(g,2)],plotcolor{1}(1))
         			    plot([j-wid j+wid],[avg(g,2) + mymult*ste(g,2) avg(g,2) + mymult*ste(g,2)],plotcolor{1}(1))
                    catch
                        disp('Can''t plot error bar.')
                        2
                        avg
                        ste
                        plotcolor
                    end
                    
                    %end
    end
    
return
    