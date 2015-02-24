function [sub,c] = getscans3(stimes,outname,ltimes)
% function [sub,c] = getscans3(input times matrix,output name,ltimes)
% fixed right now for 8 scans, for the intext study
% irf = filename for your irf scan, if any, in single quotes - if using spm gamma, leave out.
% if no irf is specified, convolution will not be done.
% ltimes = lat matrix, with 1's and 2's for left and right
%
% sub = getscans2('subject6.txt','s6model4','s6unsmooth_irf_interp.ref');
% sub = getscans2('subject6.txt','s6model4','spm');
%
% note: the program automatically adds .txt to the output name when saving the regressors, 
% so you don't have to include it.
% function:
% 			ordinary: 1 regressor for each condition
%			group: 1 regressor for each group of conditions
%					4th - nth arguments: define groups by condition $
%					e.g., [1 2] [3 4]
%           eachevt: 1 reg for each evt, random effects style...

%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up inputs
%%%%%%%%%%%%%%%%%%%%%%%%%
numscans = 8;					% number of scans
stimoffset = 0;				    % time from start of trial to stimulus onset
behoffset = 0;					% time from stim onset to cog process of interest (now not used)
scanlength = 150;				% length of scans in imgs (for cutoff)
tr = 2000;						% tr in ms
numconds = 12;                  % number of conditions per scan


%(['stimes = load(''' fname ''');'])

if stimes(1,1) < stimes(1,2)
   disp('switching columns for each regressor because they''re backwards!')
   clear a,clear b;
   for i = 1:2:15
      a = stimes(:,i);
      b = stimes(:,i+1);
      stimes(:,i) = b;
      stimes(:,i+1) = a;
   end
end

if ltimes(1,1) < ltimes(1,2)
   disp('switching lat columns for each regressor because they''re backwards!')
   clear a,clear b;
   for i = 1:2:15
      a = ltimes(:,i);
      b = ltimes(:,i+1);
      ltimes(:,i) = b;
      ltimes(:,i+1) = a;
   end
end

   disp('your result: a sample slice of the times matrix.')
   stimes(1:10,1:6)
   ltimes(1:10,1:6)
   


sub.times = stimes;
sub.ltimes = ltimes;


%%%%%%%%%%%%%%%%%%%%%%%%%
% Break into scans
% And compute mode
%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = stimes(:,1:2);
s1 = s1(s1(:,1) > 0,:);
mode(1) = max(hist(s1(:,1),s1(end,1)));
for i = 2:numscans
	eval(['s' num2str(i) ' = stimes(:,' num2str(2*i-1) ':' num2str(2*i) ');'])
   eval(['s' num2str(i) ' = s' num2str(i) '((s' num2str(i) '(:,1) > 0),:);'])
   eval(['mode(i,1) = max(hist(s' num2str(i) '(:,1),s' num2str(i) '(end,1)));'])
end

disp(['maximum mode for conditions - time responses; should be 2 for intext'])
max(mode)

disp(['Num conditions: ' num2str(numconds)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make adjustments to scans 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Adjusting Condition 6 to add 4000 ms'])
disp(['Changing trials after probes (6) to 11 (1st trial)'])
disp(['Adjusting times so that 1st event in start of each scan occurs at 0 ms and times are in TRs'])
disp(['Adding ' num2str(stimoffset) ' ms to conditions 1-4'])

% make first trial start at time 0;
for i = 1:2:size(stimes,2)
   firsttimes(i) = stimes(1,i);
end
firstevent = min(firsttimes(1,:));

cindex = 1;
for i = 1:numscans
   eval(['scan = s' num2str(i) ';'])
   lscan = ltimes(1:size(scan,1),i * 2);
   
   % add 4000 ms to condition 6, to offset it from the probe onset time.
   scan(scan(:,2) == 6,1) = scan(scan(:,2) == 6,1) + 4000;
   mode(i,1) = max(hist(scan(:,1),scan(end,1)));
     
   % change 1st trial to be 11 (1st trial)
   scan(1,2) = 11;
   
   % change trials after probes (6) to condition 7
   for j = 2:size(scan,1)
      if scan(j-1,2) == 6, scan(j,2) = 11;,end
   end
      
   % change numbers for conditions 5-6 to 9-10
   scan(scan == 5) = 9;
   scan(scan == 6) = 10;
   
   test = scan(1:10,:);
   
   % change 1-4, where lat = 2, to 5-8, and 11 where lat = 2 to 12 (1st trial L/R)
   mynums = find(scan(:,2) < 5 & lscan == 2);
   scan(mynums,2) = scan(mynums,2) + 4;
   
   mynums = find(scan(:,2) == 11 & lscan == 2);
   scan(mynums,2) = 12;
   
   disp(['scan ' num2str(i) ': orig lat output'])
   [test lscan(1:10) scan(1:10,2)] 
   
   % make first trial start at time 0;
   scan(:,1) = scan(:,1) - scan(1,1);
   
   % add stimoffset to appropriate conditions
   scan(scan(:,2) == 1,1) = scan(scan(:,2) == 1,1) + stimoffset;
   scan(scan(:,2) == 2,1) = scan(scan(:,2) == 2,1) + stimoffset;
	scan(scan(:,2) == 3,1) = scan(scan(:,2) == 3,1) + stimoffset;
	scan(scan(:,2) == 4,1) = scan(scan(:,2) == 4,1) + stimoffset;
   
   scan(:,1) = scan(:,1) ./ tr;
   scan(scan(:,1) > scanlength,:) = [];
   
   sub.scan{i} = scan;
   
   % make c variable with all times in each condition in cell array
   for j = 1:numconds
       c{cindex} = scan(scan(:,2) == j,1);
       cindex = cindex + 1;
   end
       
   
end


sub.c = c;

% save workspace
eval(['save ' outname ' sub c'])

return
