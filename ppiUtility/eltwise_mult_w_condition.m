% function crossed_data = eltwise_mult_w_condition(Y, design_matrix, wh_condition, num_runs, varargin)
%
% E.g.:
%   Y = mean(iimg_sphere_timeseries(images, XYZmm, radius)');   % to set Y to BOLD data around a certain voxel
%   design_matrix = SPM.xX.X;
%   wh_condition = 23;   % To multiply with condition 23 in each session in the design matrix
%   num_runs = length(SPM.Sess);
%   [crossed_data, regressor] = eltwise_mult_w_condition(Y, design_matrix, wh_condition, num_runs);

function [crossed_data, regressor] = eltwise_mult_w_condition(Y, design_matrix, wh_cond, num_runs, varargin)
    mean_correcting = 0;
    min_correcting = 0;

    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'mean_correct'
                    mean_correcting = 1;
                case 'nomean_correct'
                    mean_correcting = 0;
                case 'min_correct'
                    min_correcting = 1;
                case 'nomin_correct'
                    min_correcting = 0;
            end
        end
    end

    if(mean_correcting && min_correcting)
        error('Correcting for the mean and correcting for the min are mutually exclusive.');
    end

    regressor = combined_regressor(design_matrix, wh_cond, num_runs);

    if(length(regressor) ~= length(Y))
        error('Length of supplied data (%d) does not match length of design matrix (%d)', length(Y), length(regressor));
    end


    if(mean_correcting)
        Y = mean_correct_runs(Y, num_runs)
    elseif(min_correcting)
        Y = min_correct_runs(Y, num_runs)
    end

    crossed_data = Y .* regressor;
end
