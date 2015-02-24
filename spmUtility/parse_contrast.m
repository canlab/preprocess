% condition_weights = parse_contrast(contrasts, regressor_names)
%
% Returns the weights of the contrasts for each condition, in the order the
% regressor_names are passed in. For use with SPM, if all regressor_names are present
% in each run in the same order, you can just tile the results across the runs. If
% not, you will have to figure out where to place the weights manually. Or,
% if you just duplicate the regressor_names across the whole experiment
% properly, the weights will be calculated correctly.
%
%
% contrasts: cellstr array along the lines of 'A > B'
%   or 'A + B > C + D + E'
%   or '(A > B) > (C > D)'
% regressor_names: cellstr array of condition names
%
% condition_weights: matrix of weights for the regressor_names that represents the
% contrast, e.g. [ 0  1  0 -1  0]

function condition_weights = parse_contrast(contrasts, regressor_names)
    contrasts = cellstr(contrasts);
    regressor_names = cellstr(regressor_names);

    if(~iscellstr(contrasts) || ~iscellstr(regressor_names))
        error('At least one of the inputs to %s is not of type cellstr or char.', mfilename());
    end

    idx = strfind(contrasts, '>');

    for i=1:length(idx)
        if(isempty(idx{i}))
            error('Contrast string "%s" is ill-formed. (It needs a ">".)', contrasts{i});
        elseif(length(idx{i}) > 1)
            parenidx = strfind(contrasts{i}, '(');
            if(isempty(parenidx))
                error('Contrast string "%s" is ill-formed. (When contrasting contrasts, nesting is required. Uses parentheses to group. E.g. "(A > B) > (C > D)".)', contrasts{i});
            end
        end
    end

    condition_weights = zeros(length(contrasts), length(regressor_names));

    for i=1:length(contrasts)
        condition_weights(i,:) = parse_contrast_line(contrasts{i}, regressor_names);
    end
end


function condition_weights = parse_contrast_line(current_contrast, regressor_names)
    subcontrasts = regexp(current_contrast, '\(([^\)]*)\)', 'tokens');
    if(isempty(subcontrasts))
        condition_weights = parse_single_contrast(current_contrast, regressor_names);
    else
        subweights(1,:) = parse_single_contrast(subcontrasts{1}{1}, regressor_names);
        subweights(2,:) = -1 * parse_single_contrast(subcontrasts{2}{1}, regressor_names);
        condition_weights = sum(subweights, 1) / 2;
    end
end


function condition_weights = parse_single_contrast(current_contrast, regressor_names)
    condition_weights = zeros(1, length(regressor_names));
    split = regexp(strtrim(current_contrast), '(.*)>(.*)', 'tokens', 'once');
    for j=1:2
        current_conditions = strtrim(regexp(strtrim(split{j}), '([^+]+(?=[^+]*\+)|(?<=\+\s*)[^+]+)', 'match'));
        if(isempty(current_conditions))
            current_conditions = {strtrim(split{j})};
        end
        if(j == 1)
            weight = 1.0/length(current_conditions);
        else
            weight = -1.0/length(current_conditions);
        end

        for k=1:length(current_conditions)
            condition_weights(strmatch(strtrim(current_conditions{k}), regressor_names, 'exact')) = weight;
        end
    end
end
