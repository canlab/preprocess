% Y = func_correct_runs(Y, num_runs, func)
%
% Given a vector of data and a number of runs/sections to divide it into, 
% func_correct_runs corrects each run based on supplied function handle
%
% E.g.:
%   Y = func_correct_runs(amygdala_data, 6, @mean) % to remove the mean from each run
%   Y = func_correct_runs(amygdala_data, 6, @min) % to remove the min from each run

function Y = func_correct_runs(Y, num_runs, func)
    if(~isa(func, 'function_handle'))
        error('Input "func" is not a function handle');
    end

    num_vols_per_run = length(Y) / num_runs;
    for i=1:num_vols_per_run:length(Y)
        wh_run = i:i+num_vols_per_run-1;
        Y(wh_run) = Y(wh_run) - func(Y(wh_run));
    end
end