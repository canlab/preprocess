function F = compute_edf(y, t)

% computes empirical distribution function

% from Trisha van Zandt's "Analysis of Response Time Distributions" chapter
%
% input arguments: ordered column vector of RTs 'y' and a column vector of
% points 't' for which the EDF is desired
% output arguments: column vector equal in length to the vector 't'
% containing the points of the EDF

F=ones(length(t),1);
for i=1: length(t),
    F(i) = sum(y<=t  (i))/length(y);
end