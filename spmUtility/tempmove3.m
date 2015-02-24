for ss = [6]
   
   disp(['Subject ' num2str(ss)])
	eval(['cd /data4/intext2/sub' num2str(ss) '/task'])
%	!rm nra*
	!gunzip ravol*
   
end

try
   taskpreprocessintext
catch
   error('error in moving scans')
end

tor_script_stat_model5


for ss = [1 2 4 6 8 9 10 11 12]
   
   disp(['Subject ' num2str(ss)])
	eval(['cd /data4/intext2/sub' num2str(ss)])
   !mkdir Anatomy
   cd Anatomy
   try,eval(['!cp /data1/intext/sub' num2str(ss) '/Anatomy/* .']),catch,end
   try,eval(['!cp /data2/intext/sub' num2str(ss) '/Anatomy/* .']),catch,end
	try,eval(['!cp /data4/intext/sub' num2str(ss) '/Anatomy/* .']),catch,end

end

model = 5;
subs = ss;
normsmooth_modelX

model = 15;
normsmooth_modelX