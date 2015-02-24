goOK = 1;
subs = {'sub12'}; %'sub1' 'sub2' 'sub4' 'sub6' 'sub8' 'sub9' 'sub10' 'sub11'};
cd C:\Tor_Documents\CurrentExperiments\intext2
load EXPT

for snum = subs
	
	SubjCode = snum{1};

	cd C:\Tor_Documents\CurrentExperiments\intext2
	clear c
	intext_model1

end