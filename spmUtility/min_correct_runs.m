% Y = min_correct_runs(Y, num_runs, func)
%
% Given a vector of data and a number of runs/sections to divide it into, 
% min_correct_runs removes the min value from each run
%
% Mostly just a wrapper for func_correct_runs
%
% E.g.:
%   Y = min_correct_runs(amygdala_data, 6)

function Y = min_correct_runs(Y, num_runs)
    Y = func_correct_runs(Y, num_runs, @min);
end