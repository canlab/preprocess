function batch_group_outline(EXPT,gt)
% function batch_group_outline(EXPT,gt)
% Tor Wager
% EXPT is structure, as in get_expt_info
% gt is thresholds for how many subjects must activate
% to show up in a given color, e.g., [10 15 20]
% for 3 thresholds.  Colors are 'b' 'r' 'y', hard-coded
% in this script.
%
% Start in directory above individual results directories
% and run this.
%
% Thresholds are .05 = height and 5 = extent, also
% hard-coded here.

Po = which('single_subj_T1.img'); 
if isempty(Po), Po = spm_get(1,'*img','Choose anatomical overlay'),end

str = num2str(gt); str(str == ' ') = '_';

d = dir;
P = [];
mypwd = pwd;
di = 1;
while isempty(P)
	if d(di).isdir
		eval(['cd ' d(di).name])
		P = dir('con*img');
	end
	di = di + 1;
	eval(['cd ' mypwd])

	if di > length(d), error('No con*img files found in subdirectories!'), end
end

P = str2mat(P.name);
nums = str2num(P(:,end-7:end-4))';
disp('Testing these contrasts:')
nums

nnums = input('Enter vector [e.g., 2 4 5] of which to test, or return for all: ');
if ~isempty(nnums)
	nums = nnums;
	disp('Testing these contrasts:')
	nums
	pause(4)
end

for i = 1:length(nums)

	[gvol,clusters] = group_outline2(EXPT.subjects,.05,5,nums(i),gt,Po,{'b' 'r' 'y'});
	if ~isempty(clusters{1}.XYZmm)
		saveas(gcf,['group_con' num2str(nums(i)) str],'fig')
		saveas(gcf,['group_con' num2str(nums(i)) str],'jpg')
	end

end

return
