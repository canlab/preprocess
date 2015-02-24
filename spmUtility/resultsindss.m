% resultsindss:
% script by Tor Wager, 6/06/01
% -----------------------------------------------------------
% This script loops through individual subjects and prints a
% page for each that includes the ResMS image, some betas,
% and some t-maps overlaid on the betas.  
%
% for betas, pick your own beta images based on your model -
% they should be some of the intercept betas from the model.
%
% for t-maps, likewise, specify the ones from your model.
%
% the script will look in the datapath you specify.  I use:
% /data4/intext2/RESULTS/model#/sub# 
% to store the individual subjects' results for each model.

model = input('Enter model number: ');


for snum = [1 2 4 6 8 9 10 11 12]
   
   do0003 = 1;
   	do0004 = 1;

datapath = ['/data4/intext2/RESULTS/model' model '/sub' num2str(snum)];

   P = str2mat(...
   [datapath '/ResMS.img'],...
   [datapath '/beta_0038.img'],...
   [datapath '/beta_0035.img'],...
   [datapath '/beta_0035.img'] ...
   );

warning off
spm_imcalc_ui([datapath '/spmT_0002.img'], ...
   [datapath '/thresh_spmT_0002.img'],'i1 + 0./(i1>3)');

if do0003
try
	spm_imcalc_ui([datapath '/spmT_0003.img'], ...
   [datapath '/thresh_spmT_0003.img'],'i1 + 0./(i1>3)');
catch
end
end

if do0004
try
	spm_imcalc_ui([datapath '/spmT_0004.img'], ...
   [datapath '/thresh_spmT_0004.img'],'i1 + 0./(i1>3)');
catch
end
end


warning on

spm_check_registration(P);
spm_orthviews('Interp',0);


spm_orthviews('window',1,[-30 30])


spm_orthviews('addimage',2,[datapath '/thresh_spmT_0002.img']);

if do0003
try
	spm_orthviews('addimage',3,[datapath '/thresh_spmT_0003.img']);
catch   
end
end

if do0004
try
	spm_orthviews('addimage',4,[datapath '/thresh_spmT_0004.img']);
catch   
end
end


figure(1)

% text(-200,250,['Model ' model ' sub ' num2str(snum) ' T0002 and 0003'],'FontSize',18,'Color','b')

gtext(['Model ' model ' sub ' num2str(snum) ' T0002, 0003, 0004'],'FontSize',18,'Color','b')

input('Position crosshairs and press return.')

spm_figure('Print')

end
