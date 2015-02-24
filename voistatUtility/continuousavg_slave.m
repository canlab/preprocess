function [avg,ste,bincenter] = continuousavg_slave(y,window,groups,binres,time,flindex)
% [avg,ste,bincenter] = continuousavg_slave(y,window,groups,binres,time,flindex)


windowlen = window(2) - window(1);
myresp = [];

if ~exist('groups'), groups = 1:size(events,2);,end

% 'slave' voxels: pass in time and flindex; compute response.
for i = 1:max(groups)
   for j = 1:size(flindex{i},1)
      myresp = [myresp; y(flindex{i}(j,1):flindex{i}(j,2))];				% response at time first:last. 
   end
   resp{i} = myresp;
   myresp = [];
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
