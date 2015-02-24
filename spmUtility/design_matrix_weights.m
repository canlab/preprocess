% weights = design_matrix_weights(conditions, condition_weights, design_matrix_columns)
%
% E.g.:
% contrasts = {'A > B' 'A + B > C + D + E'};
% conditions = {'A' 'B' 'C' 'D' 'E'};
% condition_weights = parse_contrast(contrasts, conditions);
% ...
% SPM = spm_fmri_spm_ui(SPM); % This must have been run already or else the design matrix won't be filled out
% ...
% weights = design_matrix_weights(conditions, condition_weights, SPM.xX.name);

function weights = design_matrix_weights(conditions, condition_weights, design_matrix_columns)
    num_contrasts = size(condition_weights, 1);
    num_conditions = length(conditions);
    num_design_matrix_columns = length(design_matrix_columns);

    weights = zeros(num_contrasts, num_design_matrix_columns);

    for i=1:num_contrasts
        for j=1:num_conditions
            if(condition_weights(i,j) ~= 0)
                wh = ~cellfun(@isempty, regexp(design_matrix_columns, [conditions{j} '\>'], 'match', 'once'));
                weights(i, wh) = condition_weights(i, j);
            end
        end
    end
end