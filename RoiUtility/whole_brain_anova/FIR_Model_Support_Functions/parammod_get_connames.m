function EXPT = parammod_get_connames(EXPT)
%
% 
% uses contrast names entered in EXPT.RT.connames
% assumes those names are used as basis for running
% wb_multisubject_parammod
% gets contrast images with same names
% and stores them in EXPT.SNPM.P for random effects analysis with robfit.
%
% tor wager


EXPT.SNPM.connames = str2mat(EXPT.RT.connames{:});
EXPT.SNPM.connums = 1:length(EXPT.RT.connames);
EXPT.SNPM.P = {};

for i = 1:length(EXPT.RT.connames)
    
    EXPT = getfunctnames2(EXPT,[EXPT.RT.connames{i} '.img'],'tmp');
    tmp = str2mat(EXPT.tmp{:});
    
    if ~isempty(tmp)
        EXPT.SNPM.P{i} = tmp;
    else
        warning('Contrast names empty!! No contrast images?');
    end
    
end


return
