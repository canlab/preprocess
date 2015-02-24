function actual_temp = tsa_temp_correction(nominal_temp, which_calib_set, varargin)
% actual_temp  = tsa_temp_correction(nominal_temp, which_calib_set, [enter any argument to plot])
% 
% This function adjusts nominal temperature entered into COVAS to reflect
% actual temp. delivered to skin. You need to use this, with appropriate
% calibration data, if you really care about exactly which temperature was
% delivered (e.g., when comparing subjects across different TSA setups).
% 
% A new calibration table should be created each time the TSA is
% calibrated.  The table will also be different for the 16 mm vs. 30 mm
% thermodes.
%
% Valid entries for which_calib_set:
% '16mm_TSA_no_ECU_April07'
% '16mm_TSA_with_ECU_April07'
% '16mm_TSA_with_ECU_Sep08'
% '16mm_TSA_no_ECU_Oct08'
% '30mm_TSA_no_ECU_Oct08'
% '16mm_TSA_with_ECU_Oct08_long'
% 
% NOTE: Prior to 12/4/07, the TSA used for scans was the ECU. The TSA used
% for behaviorals was non-ECU (this applies to Atlas and Leotti's studies.  
% DPSP may have overlapped as well). On 12/4/07, we switched these, so
% scans after that point are done on non-ECU, and behavioral sessions are
% done with ECU.
%
% 
% NOTE: could update to bilinear (broken stick) regression for small
% improvement in accuracy, esp with ECU
%
% Examples:
% actual_temp = tsa_temp_correction(nsfpX1data.temps{1}, '16mm_TSA_no_ECU_April07', 1)
% actual_temp = tsa_temp_correction(mytemps, '16mm_TSA_with_ECU_April07') % do not plot

doplot = 0;
if length(varargin) > 0, doplot = 1; end

switch which_calib_set

    case '16mm_TSA_no_ECU_April07'

        % calibration data
        nom = [32 40 45 46 46.5 47:.5:49 50];
        act = [31.65 39.05 43.7 44.95 45.4 45.9 46.4 46.94 47.4 47.94 48.75];

    case '16mm_TSA_with_ECU_April07'

        nom = [32 40 45 46 47 47.5 48 49 50];
        act = [31 39.1 44.9 45.4 46.8 47.25 47.7 48.6 49.8];

    case '16mm_TSA_with_ECU_Sep08'
        nom = [32    46    47    48    49    50];
        act =  [32.35   46.50   47.76   48.77   49.77   50.74];

    case '16mm_TSA_no_ECU_Oct08'
        nom = [32   38   40   42   44   46   46.50   47 ...
            47.50   48   48.50   49   49.50   50];
        act = [32.23   37.96   39.81   41.68   43.65   45.52   45.98   46.51 ...
            47   47.65   48.10   48.68   49.24   49.85];

    case '30mm_TSA_no_ECU_Oct08'
        nom = [32   38   40   42   44   46   46.50   47 ...
            47.50   48   48.50   49   49.50   50];
        act = [32.18   38.01   39.95   41.90   43.85   45.76   46.30   46.81 ...
            47.33   47.89   48.44   48.99   49.60   50.18];
        
    case '16mm_TSA_with_ECU_Oct08_long'
        nom = [32   38   40   42   44   45   46   46.50 ...
            47   47.50   48   48.50   49   49.50   50]
        act = [32.19   38.28   40.40   43.45   44.75   45.65   46.68   47.45 ...
            47.85   48.34   48.70   49.38   49.90   50.35   50.90]

    otherwise error('Invalid entry for which_calib_set')

end




% Get the linear transformation (scaling factor) for predicting actual
% from nominal, given the calibration data above.
% -----------------------------------------------
b = pinv(nom') * act';  % no-intercept model


% apply this linear scaling factor

actual_temp = nominal_temp' * b;

% Plot
if doplot
    create_figure('Calibration', 1, 2); title('Calibration values')
    plot(nom, act, 'ko')
    f = nom * b;
    hold on; plot(nom, f, 'r')

    subplot(1, 2, 2)
    plot(nominal_temp, 'ko-', 'LineWidth', 2);
    plot(actual_temp, 'rs-', 'LineWidth', 2);
end


end
