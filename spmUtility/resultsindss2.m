% resultsindss2:
% script by Tor Wager, 6/17/01
% -----------------------------------------------------------
% This script loops through individual subjects and prints a
% page for the whole group that has the contrast of interest
% for each subject, overlaid on the beta image of your choice.  
%
% for betas, pick your own beta images based on your model -
% they should be some of the intercept betas from the model.
%
% for t-maps, likewise, specify the ones from your model.
%
% the script will look in the datapath you specify.  I use:
% /data4/intext2/RESULTS/model#/sub# 
% to store the individual subjects' results for each model.
%
% initial values for P are templates, etc. to compare against

disp('resultsindss2')
disp('---------------------------')
disp('Exampless:')
disp('Enter model number: 15')
disp('Enter contrast number: 0003')
disp('Enter img file stem: nspmT')
disp('Enter t-threshold: 3')
disp('     *      *     *')

model = input('Enter model number: ','s');
connum = input('Enter contrast number: ','s');
imgbase = input('Enter img file stem: ','s');
thresh = input('Enter t-threshold: ','s');

subjects = [1 2 4 6 8 9 10 11 12];
studypath = '/data4/intext2';

% ---------------------------------
% * build list of background images
% ---------------------------------

P = str2mat(['/data3/biman/spm99/spm99/templates/T1.img']);
for snum = subjects   

	datapath = [studypath filesep 'sub' num2str(snum) filesep 'Anatomy'];
	P = strvcat(P,str2mat([datapath filesep 'nrt1.img']));
      
end
   
% -----------------------
% * calculate thresholded
% -----------------------

warning off

for snum = subjects   

	datapath = [studypath filesep 'RESULTS/model' model '/sub' num2str(snum)];
   fname = [datapath filesep 'thresh_' thresh '_' imgbase '_' connum '.img'];
   
   % eval(['cd ' datapath])
   timg = getfiles(fname(1:end-4));
   
   if isempty(timg)
   
		spm_imcalc_ui([datapath filesep imgbase '_' connum '.img'], ...
   		fname,['i1 + 0./(i1>' thresh ')']);

   end
   
end


% -----------------------
% * calculate group mean
% -----------------------

meanname = [studypath filesep 'RESULTS/model' model filesep 'mean_' imgbase '_' connum '.img'];
meanthrname = [studypath filesep 'RESULTS/model' model filesep 'mean_thr_' thresh '_' imgbase '_' connum '.img'];

timg = getfiles(meanname(1:end-4));
   
if isempty(timg)
   
	 
 	% -------------------------------------------
 	% * build list of component parametric maps
 	% -------------------------------------------

    Q = [];
    calcstr = [];
	index = 1;
	for snum = subjects   

		datapath = [studypath filesep 'RESULTS/model' model '/sub' num2str(snum)];
   		fname = [datapath filesep 'thresh_' thresh '_' imgbase '_' connum '.img'];
   
   		Q = strvcat(Q,str2mat(fname));
   		calcstr = [calcstr ' i' num2str(index) ' + '];
         
      index = index + 1;
         
	end

	calcstr = calcstr(1:end-3);		% get rid of final +
	calcstr = [calcstr ' ./ ' num2str(length(subjects))];
   
   disp(['resultsindss2: calculating mean image: ' meanname])
   disp(['calculation string is: ' calcstr])
   
   spm_imcalc_ui(Q, meanname, calcstr);
   
end


% -----------------------
% * threshold group mean
% -----------------------


timg = getfiles(meanthrname(1:end-4));
   
if isempty(timg)

	spm_imcalc_ui(meanname, meanthrname, ['i1 + 0./(i1>' thresh ')']);

end

warning on


% -----------------------
% * display background
% -----------------------

spm_check_registration(P);
spm_orthviews('Interp',0);



% -----------------------
% * add thresholded
% -----------------------

spm_orthviews('addimage',1,meanthrname);

index = 2;
for snum = subjects
   
   datapath = [studypath filesep 'RESULTS/model' model '/sub' num2str(snum)];
   fname = [datapath filesep 'thresh_' thresh '_' imgbase '_' connum '.img'];

   spm_orthviews('addimage',index,fname);
   index = index + 1;
   
end




% -----------------------
% * text labels and print
% -----------------------

figure(1)

% text(-200,250,['Model ' model ' sub ' num2str(snum) ' T0002 and 0003'],'FontSize',18,'Color','b')

gtext(['Model ' model ' ' num2str(length(subjects)) ' subj. ' imgbase '_' connum],'FontSize',14,'Color','b')

input('Position crosshairs and press return.')

spm_figure('Print')


