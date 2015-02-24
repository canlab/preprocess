for ss = [1 2 4 6 8 9 10 11 12]
   
   disp(['Subject ' num2str(ss)])
	eval(['cd /data4/intext2/sub' num2str(ss) '/task'])
	!rm nra*
	!gunzip ravol*
   
end

try
   taskpreprocessintext
catch
   error('error in moving scans')
end

tor_script_stat_model15

