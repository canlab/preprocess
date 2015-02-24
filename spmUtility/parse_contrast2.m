function condition_weights = parse_contrast2(contrasts, conditions)
% condition_weights = parse_contrast(contrasts, conditions)
%
% Returns the weights of the contrasts for each condition, in the order the
% conditions are passed in. For use with SPM, if all conditions are present
% in each run in the same order, you can just tile the results across the runs. If
% not, you will have to figure out where to place the weights manually. Or,
% if you just duplicate the conditions across the whole experiment
% properly, the weights will be calculated correctly.
%
%
% contrasts: cellstr array along the lines of 'A > B'
%   or 'A + B > C + D + E'
% conditions: cellstr array of condition names
%
% condition_weights: matrix of weights for the conditions that represents the
% contrast, e.g. [ 0  1  0 -1  0]

warning('WARNING: parse_contrast2 is deprecated! All future development has been folded back into parse_contrast. Use that instead.');

condition_weights = parse_contrast(contrasts, conditions);
end