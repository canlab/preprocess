function [X2, e] = chi2(table)

% [X2, expected] = chi2(table)
% calculates the chisquare statistic on any 2d table
% expected is an optional output, which is the table of expected values

e = zeros(size(table));
X2 = 0;
for i = 1:size(table,1)
    for j = 1:size(table,2)
        e(i,j) = sum(table(i,:)) * sum(table(:,j)) / sum(sum(table));
        X2 = X2 + ((table(i,j) - e(i,j)) ^ 2) / e(i,j);
    end
end

