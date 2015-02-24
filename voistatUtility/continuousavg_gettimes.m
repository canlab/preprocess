function [time,flindex,bincenter] = continuousavg_gettimes(y,window,events,groups,binres)
% function [time,flindex,bincenter] = continuousavg_gettimes(y,window,events,groups,binres)
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
tslen = length(y);

if ~exist('groups'), groups = 1:size(events,2);,end

%disp('Start concatenation');drawnow
% concatenate events in groups
% ============================================================
for i = 1:max(groups)
   grpevents = [];
   for j = 1:size(groups,2)
		if groups(j) == i, grpevents = [grpevents;events{j}]; end
   end
   newevents{i} = grpevents(grpevents>0,:);
end


% make long list of times from event onset and responses for each group
% ============================================================
for i = 1:max(groups)
   list = newevents{i};
   index = 1;
   for j = 1:size(list,1)
      if list(j,1) > 0
      		first = max(round(list(j,1)) + window(1),1);
          last = min(round(list(j,1)) + window(2),length(y));
          flindex{i}(j,1) = first;
          flindex{i}(j,2) = last;
      		for k = first:last
             %if k > 0 & k <= size(y,1)								% make sure not to go off the end.
            		time{i}(index) = k - list(j,1);					% time relative to stim onset
                  %resp{i}(index) = y(k);								% response at time k.

            		index = index + 1;
         		%end
          end
      end
   end
end

% define bin centers to get number of bins.
% ============================================================
for i = 1:max(groups)
	index = 1;
   for j = window(1):binres:window(2)
      bincenter{i}(index) = j + .5 * binres;
      index = index +1;
   end
end


return
