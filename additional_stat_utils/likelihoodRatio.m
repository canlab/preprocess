function [X2] = likelihoodRatio(table)

% [X2] = likelihoodRatio(table)
% Give the chi2 statistic from a contingency table

e = zeros(size(table));
for i = 1:size(table,1)
    for j = 1:size(table,2)
        e(i,j) = sum(table(i,:)) * sum(table(:,j)) / sum(sum(table));
    end
end

X2 = 2 * sum(sum( table .* log( table ./ e ) ));

