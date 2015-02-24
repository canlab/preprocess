% user-prompted input
% generic script for checking normalization
% across multiple subjects
%
% Tor Wager, 12/5/01



substyle = 2;		% name of ind sub directory - choose 1 for numbers, 2 for letter/number string.

disp('1	registration')
disp('2	normalization')
disp('3	reg + norm')
cchoice = input('Enter your choice of which to check: ');

if ~(exist('subs') == 1), error('Define subs cell array or vector before running.'), end

if ~(exist('canonicalT1') == 1),
	canonicalT1 = spm_get(1,'img','Choose canonical T1',[],0);
end

if ~(exist('fmriDIR') == 1),
	fmriDIR = '/data/biman5';
end

%if ~(exist('segtemplate') == 1),
%	segtemplate = spm_get(1,'img','Choose segmented canonical T1',[],0);
%end

eval(['cd ' fmriDIR])

myt1 = spm_get(1,'img','Choose reg, UNnorm. T1 from one subject',[],0);

if cchoice == 1 | cchoice == 3
	mytarget = spm_get(1,'img','Choose target image for reg for one subject',[],0);
end

clear L, clear M
L{1} = canonicalT1;
M{1} = canonicalT1;

% get image and path
% ----------------------------------------------
[dummypath,rt1name,t1ext] = fileparts(myt1);
rt1name = [rt1name t1ext];

if cchoice == 2 | cchoice == 3
	nrt1name = ['n' rt1name];
end

if cchoice == 1 | cchoice == 3
	[dummypath,targetname,targetext] = fileparts(mytarget);
	targetname = [targetname targetext];
end



	

index = 1;
for i = subs

	if substyle == 2
   		% if subs get letter/number code   
   		snum = i{1};
	elseif substyle == 1
   		% if subs are numbered
   		snum = ['sub' num2str(i)];
   	else
		error('Unknown substyle - choose 1 or 2')
	end

	t1path = [snum filesep 'anatomy'];
	targetpath = [snum filesep 'anatomy'];

	if cchoice == 2 | cchoice == 3
		%segt1 = spm_get(1,'img',['Choose funct image for ' snum],[],0);
		a = getfiles(['/data/biman5/' snum '/task/scan1/ra*0001.img']);
		segt1 = a{1};
		[segt1path,segt1name,segext] = fileparts(segt1);
		segt1name = [segt1name segext];

		L{index+1} = [t1path filesep nrt1name];
		M{index+1} = [segt1path filesep segt1name]; 
		index = index + 1;
	end

	if cchoice == 1 | cchoice == 3
		clear N
		N{1} = canonicalT1;
		N{2} = [targetpath filesep targetname];
		N{3} = [t1path filesep rt1name];
		N = str2mat(N);
		N
		str = ['Registration (canonical, target, T1 (object)) for ' snum];
		spm_check_registration(N);
		%spm_orthviews('Interp',0);
		gtext(str,'FontSize',14,'Color','b')
		str = input('Press return to continue - will print where xhairs are...');
		spm_print
	end

end

L = str2mat(L);
M = str2mat(M);
L
spm_check_registration(L);
%spm_orthviews('Interp',0);
str = input('Position crosshairs and then input text title for this page: ','s');
gtext(str,'FontSize',14,'Color','b')

spm_print

M
spm_check_registration(M);
spm_orthviews('Interp',0);
str = input('Position crosshairs and then input text title for this page: ','s');
gtext(str,'FontSize',14,'Color','b')

spm_print