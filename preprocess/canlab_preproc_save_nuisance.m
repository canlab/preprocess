function canlab_preproc_save_nuisance(PREPROC)
% canlab_preproc_save_nuisance(PREPROC)
% Build a set of nuisance regressors and save in matrix R in .mat file
% for input into SPM batch model jobs.

N = fieldnames(PREPROC.Nuisance);

if isempty(N) || ~iscell(N)
    disp('Did not save nuisance regressor .mat file: No fields in PREPROC.nuisance.');
    return
end

for i = 1:length(N)
    R{1, i} = PREPROC.Nuisance.(N{i});
    v(1, i) = size(R{1, i}, 1);
end

if any(diff(v))
    error('PREPROC.Nuisance fields are not the same size. Should all equal num. volumes in study.');
end

outname = fullfile(PREPROC.basedir, 'Functional', 'Preprocessed', 'Nuisance_covariates_R.mat');
save(outname, 'R');

fprintf('Saved nuisance regressors in R variable in this file:\n%s\n', outname);

end % function
