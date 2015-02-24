goOK = 1;

EXPT.mask = '/usr/private/spm99/templates/scalped_avg152T1_graymatter_smoothed.img'
cd /data/placebo
			%load EXPT

% load EXPT and check image files 
% to make sure they match subdirs
% ===============================
for i = 1:length(EXPT.subjects), 
	a(i) = any(findstr(EXPT.im_files{i}(1,:),EXPT.subjects{i}));,
end
if any(~a), error('Check image files in EXPT!'), end
 
subs = EXPT.subjects

for snum = 1:length(EXPT.subjects)
	
	SubjCode = EXPT.subjects{snum};

	cd /data/placebo
	clear c

	if strcmp(EXPT.order{snum},'CACA')
		placebo_model3_manip_CACA
    elseif strcmp(EXPT.order{snum},'ACAC')
		placebo_model3_manip_ACAC
	end

	cd /data/placebo
	clear c
	if strcmp(EXPT.order{snum},'CACA')
		placebo_model3_test_CACA
    elseif strcmp(EXPT.order{snum},'ACAC')
		placebo_model3_test_ACAC
	end
end