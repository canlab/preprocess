function events = cond2evt(file,numscans,boldreps)
% events = cond2evt(file,numscans,boldreps)
% boldreps is images acquired per scan

% load file and get the numbers
% ====================================
f = textread(file,'%s');
startf = max(find(strcmp(f,';'))) + 1;
endf = length(f);

index = 1;
for i = startf:endf
   cond(1,index) = str2num(f{i});
   index = index + 1;
end

events.cond = cond;


% check for divisibility and stuff
% ====================================
leng = length(cond);
imgsperscan = leng/numscans;

if ~(round(imgsperscan) == imgsperscan)
     warning('Cond function is not equally divisible into scans.')
     leng
     numscans
     imgsperscan
     error('exiting now...')
end

totalbr = boldreps * numscans;
elsperscan = leng / totalbr;

events.boldreps = totalbr;
events.elsperscan = elsperscan;


if ~(elsperscan == 1)
     warning([elsperscan ' elements per scan detected.'])
end

times = 0:1/elsperscan:totalbr-1;

if ~(length(times) == leng)
     warning('times length does not equal elements length.')
     length(times)
     leng
     cond(1:10)
     times(1:10)
     cond(end-10:end)
     times(end-10:end)
end

numconds = max(cond);


for i = 1:numscans
   % ---- define the scan stuff ----
   scanstart = (i-1) * boldreps * elsperscan + 1; 
   scanend = i * boldreps * elsperscan;

   disp(['Scan ' num2str(i) ': elements ' num2str(scanstart) ' - ' num2str(scanend) ' of ' num2str(leng)])

   scanc = cond(scanstart:scanend);
   scant = 0:1/elsperscan:boldreps-1;
   
   if ~(length(scanc) == length(scant))
      warning(['Scan ' num2str(i) ': times and cond list are different lengths.'])
   end
   
   % ---- determine max event time length ----
   for j = 1:numconds
      maxlen(j) = max(sum(scanc == j));
   end
   maxlen = max(maxlen);

   for j = 1:numconds
      % ----- get event times -----
      evttimes = scant(scanc == j);
      addnum = maxlen - length(evttimes);
      
      % ----- add buffers -----
      if addnum > 0, evttimes = [evttimes -1*ones(1,addnum)];,end
      
      % ----- store in output -----
      events.offset{i}(j,:) = evttimes;
   end % loop through conds

end % loop through scans

