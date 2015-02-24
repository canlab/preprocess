
calcacross = 0;
randthresh = 3;

%ss = [1 2 4 6 8 9 11 12];
%drive = [1 1 1 1 2 2 4 4];

printfile = 'spm99.ps';

global PRINTSTR
PRINTSTR = [spm_figure('DefPrintCmd'),printfile];

for i= 1:size(ss,2)
   
   if i > 1,calcacross = 0;,end
 
   datapath = ['/data4/intext2/RESULTS/model' model '/sub' num2str(ss(i))];
   checkss
   spm_figure('Print')
end