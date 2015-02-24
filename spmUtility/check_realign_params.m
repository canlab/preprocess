% for SPM realignment
%
% requires subject number codes
%
% start in the main study directory.
% run this.

for i = 1:length(subs)

	mysub = num2str(subs(i));
	taskpath = ['sub' mysub filesep 'task'];

	[file,dirs] = spm_list_files(taskpath,'real*');
	file = file(1,:);
	realfile = [taskpath filesep file];
	eval(['load ' realfile])
	
	len = min(31,length(file));

	realvar = file(1:len);
	figure(findobj('Tag', 'Graphics'))
	clf
	eval(['plot(' realvar ')'])
	title(taskpath,'FontSize',18)
	xlabel(taskpath,'FontSize',18)
	legend({'x' 'y' 'z' 'roll?' 'pitch?' 'yaw?'})
	spm_print

end

!gv spm99.ps