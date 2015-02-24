function EXPT = get_correl_zimgs_for_robfit(EXPT)
% EXPT = get_correl_zimgs_for_robfit(EXPT)
%
% start in dir above individual subject dirs
% each subdir contains z-maps for within subject correlations
% saved from whole_brain_correl  / wb_multisubject_correl.m
%
% saves z-images in EXPT.SNPM.P for random effects analysis
% with robfit.m or robseed.m
%
% tor wager


EXPT.SNPM = [];

d = dir([EXPT.subjects{1} filesep 'correl_z_0*.img']); 

for i = 1:length(d)
    EXPT.tmp = [];
    EXPT = getfunctnames2(EXPT,d(i).name,'tmp');
    EXPT.tmp = str2mat(EXPT.tmp);
    EXPT.SNPM.P{i} = EXPT.tmp;
end

EXPT.SNPM.connames = str2mat(d.name);
EXPT.SNPM.connums = 1:length(d);

EXPT = rmfield(EXPT,'tmp');

if isfield(EXPT,'seednames'), EXPT.SNPM.conimgs = EXPT.SNPM.connames; EXPT.SNPM.connames = EXPT.seednames;,end

EXPT.SNPM.connames = str2mat(EXPT.SNPM.connames{:});

return
