% regressor = combined_regressor(design_matrix, wh_cond, num_runs)

function regressor = combined_regressor(design_matrix, wh_cond, num_runs)
    [num_vols num_regressors] = size(design_matrix);
    num_regressors_per_run = (num_regressors - num_runs) / num_runs;
    num_vols_per_run = num_vols / num_runs;

    wh_regressors = wh_cond:num_regressors_per_run:(num_regressors - num_runs);

    regressor = zeros(num_vols,1);
    for i=1:length(wh_regressors)
        wh_vols = (1:num_vols_per_run) + ((i-1)*num_vols_per_run);
        regressor(wh_vols) = design_matrix(wh_vols, wh_regressors(i));
    end
end