% user-prompted input
% generic script for checking normalization
% across multiple subjects
%
% Tor Wager, 12/5/01

substyle = 2;		% name of ind sub directory - choose 1 for numbers, 2 for letter/number string.


if ~(exist('subs') == 1), error('Define subs cell array or vector before running.'), end

if ~(exist('canonicalT1') == 1),
	canonicalT1 = spm_get(1,'img','Choose canonical T1',[],0);
end

if ~(exist('segtemplate') == 1),
	segtemplate = spm_get(1,'img','Choose segmented canonical T1',[],0);
end

if ~(exist('myt1') == 1),
	myt1 = spm_get(1,'img','Choose registered (un-normalized) T1 from one subject',[],0);
end

if ~(exist('segt1') == 1),
	segt1 = spm_get(1,'img','Choose other r (unnorm) image',[],0);
end

clear L, clear M
L{1} = canonicalT1;
M{1} = segtemplate;

% get name of segmented t1 image.
% ----------------------------------------------
[dummypath,rt1name,t1ext] = fileparts(myt1);
rt1name = [rt1name t1ext];
[dummypath,segt1name,segext] = fileparts(segt1);
segt1name = [segt1name segext];


% add n to t1 and segt1
% ----------------------------------------------
rt1name = ['n' rt1name];
segt1name = ['n' segt1name];	

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

	L{index+1} = [t1path filesep rt1name];
	M{index+1} = [t1path filesep segt1name]; 
	index = index + 1;
end

L = str2mat(L);
M = str2mat(M);

spm_check_registration(L);
spm_orthviews('Interp',0);
str = input('Position crosshairs and then input text title for this page: ','s');
gtext(str,'FontSize',14,'Color','b')

spm_print


spm_check_registration(M);
spm_orthviews('Interp',0);
str = input('Position crosshairs and then input text title for this page: ','s');
gtext(str,'FontSize',14,'Color','b')

spm_print