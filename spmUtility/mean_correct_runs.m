% Y = mean_correct_runs(Y, num_runs)
%
% Given a vector of data and a number of runs/sections to divide it into, 
% mean_correct_runs removes the mean from each run
%
% Mostly just a wrapper for func_correct_runs
%
% E.g.:
%   Y = mean_correct_runs(amygdala_data, 6)

function Y = mean_correct_runs(Y, num_runs)
    Y = func_correct_runs(Y, num_runs, @mean);
end