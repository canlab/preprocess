function [time,resp,avg,ste,bincenter] = continuousavg(y,window,events,groups,binres,varargin)
% function [time,resp,avg,ste,bincenter] = continuousavg(y,window,events,groups,binres,varargin)
%
% y:		 time series of voxel or region data, column vector
%
% window: scans prior to and succeeding event time to collect averages over.  
%			 0 is the closest scan to event onset
%
% events: a cell array. In each cell,1 column vector containing onsets in scans for each condition
%			 event times start with 0.
%
% groups: a row vector of which event columns go in which group.
%			 like - [1 1 2 2 2 3 3] for 7 conditions in 3 groups.
%
% variable args:
%			 'voifig', followed by voifigh (handles for voifig figure)
%
% example:
% [time,resp,avg,ste,bincenter] = continuousavg(y,[-1 6],events,[1 1 2 2],.5,'voifig',voifigH)
%
% by Tor Wager, 4/16/01

windowlen = window(2) - window(1);
colors = {'bs','rd','go','c^','kv','mo'};
plotpoints = 0;
wid = .1;
lwid = 2;

index = 1;
while index <= size(varargin,2)
   switch varargin{index}
   case 'voifig', voifigH = varargin{index+1};,index = index+1;
   case 'nopoints', plotpoints = 0;
   end
   index = index + 1;
end

if ~exist('groups'), groups = 1:size(events,2);,end

disp('Start concatenation');drawnow
% concatenate events in groups
% ============================================================
for i = 1:max(groups)
   newevents{i} = [];
   for j = 1:size(groups,2)
		if groups(j) == i, newevents{i} = [newevents{i};events{j}]; end
   end
   newevents{i} = newevents{i}(newevents{i}>0,:);
end


disp('Start event listing');drawnow
% make long list of times from event onset and responses for each group
% ============================================================
for i = 1:max(groups)
   list = newevents{i};
   index = 1;
   for j = 1:size(list,1)
      if list(j,1) > 0
      		first = round(list(j,1)) + window(1);
      		last = round(list(j,1)) + window(2);
      		for k = first:last
             if k > 0 & k <= size(y,1)								% make sure not to go off the end.
            		time{i}(index) = k - list(j,1);					% time relative to stim onset
            		resp{i}(index) = y(k);								% response at time k.
            		index = index + 1;
         		end
          end
      end
   end
end

disp('Start averaging');drawnow
% average responses in specified bins
% ============================================================
for i = 1:max(groups)
	index = 1;
   for j = window(1):binres:window(2)
      bincenter{i}(index) = j + .5 * binres;
      avg{i}(index) = mean(resp{i}(time{i}>j & time{i}<j + binres));
      ste{i}(index) = std(resp{i}(time{i}>j & time{i}<j + binres)) ./ sqrt(sum(time{i}>j & time{i}<j + binres));
      index = index +1;
   end
end

% plotting stuff
% ============================================================
if exist('voifigH'), axes(voifigH(4)), else figure, end
   
for i = 1:max(groups)
   hold on
   if plotpoints,plot(time{i},resp{i},colors{i}),end
   plot(bincenter{i},avg{i},colors{i}(1),'LineWidth',lwid)
   index = 1;
   for j = bincenter{i}
      plot([j j],[avg{i}(index)-ste{i}(index) avg{i}(index)+ste{i}(index)],colors{i}(1),'LineWidth',lwid)
      plot([j-wid j+wid],[avg{i}(index)-ste{i}(index) avg{i}(index)-ste{i}(index)],colors{i}(1),'LineWidth',lwid)
      plot([j-wid j+wid],[avg{i}(index)+ste{i}(index) avg{i}(index)+ste{i}(index)],colors{i}(1),'LineWidth',lwid)
      index = index + 1;
   end
end

return
