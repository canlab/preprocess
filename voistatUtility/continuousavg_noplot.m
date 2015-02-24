function [avg,ste,bincenter] = continuousavg_noplot(y,window,events,groups,binres)
% function [avg,ste,bincenter] = continuousavg_noplot(y,window,events,groups,binres)
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


%disp('Start event listing');drawnow
% make long list of times from event onset and responses for each group
% ============================================================
%t1 = cputime;
for i = 1:max(groups)
   list = newevents{i} > 0;
   index = 1;
   for j = 1:size(list,1)
      		first = max(round(list(j,1)) + window(1),0);
      		last = min(round(list(j,1)) + window(2),tslen);
      		for k = first:last
             if k > 0 & k <= size(y,1)								% make sure not to go off the end.
            		time{i}(end+1:) = k - list(j,1);					% time relative to stim onset
                resp{i}(index) = y(k);								% response at time k.
                flindex{i}(j,1) = first; 
                flindex{i}(j,2) = last; 
            		index = index + 1;
         		end
          end
   end
end

% 'slave' voxels: pass in time and flindex; compute response.
for i = 1:max(groups)
	for j = 1:size(list,1)
		resp{i} = [resp{i} y(flindex{i}(j,1):flindex{i}(j,2))];				% response at time first:last.
	end
end

   
%disp(['Event listing took ' num2str(cputime - t1) 's']);

%disp('Start averaging');drawnow
% average responses in specified bins
% ============================================================
%t1 = cputime;
for i = 1:max(groups)
	index = 1;
   for j = window(1):binres:window(2)
      bincenter{i}(index) = j + .5 * binres;
      avg{i}(index) = mean(resp{i}(time{i}>j & time{i}<j + binres));
      ste{i}(index) = std(resp{i}(time{i}>j & time{i}<j + binres)) ./ sqrt(sum(time{i}>j & time{i}<j + binres));
      index = index +1;
   end
end
%disp(['Averaging took ' num2str(cputime - t1) 's']);


return
