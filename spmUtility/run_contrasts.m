% run_contrasts(regressor_names, contrasts, SPMs, ['allow_unbalanced', 0], ['replace', replacement_range], ['remove'])
%   Runs contrasts for the list of SPM files given
%   By default, it appends contrasts to the existing list, but if the 'replace' option is used,
%       it'll replace existing contrasts
%
% Required Inputs
%   regressor_names - names of regressors in GLM
%   contrasts - contrasts between conditions in regressor_names
%   SPMs - array of SPM structures or cellstr of SPM.mat files to add contrasts to
%
% Optional Parameters
%   allow_unbalanced - boolean - allows contrasts to be unbalanced on one side or the other, e.g., 'condition1 > ' - default 0
%   replace - replace an existing range of contrasts in SPM.xCon - default off
%   remove - wipe out all existing contrasts - default off
%
% E.g.:
%   regressor_names = {'accept_mem_cue', 'accept_q1', 'accept_q2', 'accept_arrows', 'accept_strat_cue1', 'accept_strat_cue2', 'accept_strat_cue3', ...
%      'feel_mem_cue', 'feel_q1', 'feel_q2', 'feel_arrows', 'feel_strat_cue1', 'feel_strat_cue2', 'feel_strat_cue3', ...
%      'analyze_mem_cue', 'analyze_q1', 'analyze_q2', 'analyze_arrows', 'analyze_strat_cue1', 'analyze_strat_cue2', 'analyze_strat_cue3'};
%   contrasts = {};
%   contrasts{end+1}        = 'accept_strat_cue1 > feel_strat_cue1';
%   contrasts{end+1}        = 'accept_strat_cue1 > analyze_strat_cue1';
%   contrasts{end+1}        = 'feel_strat_cue1 > accept_strat_cue1';
%   contrasts{end+1}        = 'feel_strat_cue1 > analyze_strat_cue1';
%   contrasts{end+1}        = 'analyze_strat_cue1 > accept_strat_cue1';
%   contrasts{end+1}        = 'analyze_strat_cue1 > feel_strat_cue1';
%   contrasts{end+1}        = 'analyze_strat_cue1 + analyze_strat_cue1 + analyze_strat_cue1 > ';
%   SPMs = filenames('results_model5/afa*/SPM.mat');
%   run_contrasts(regressor_names, contrasts, SPMs, 'allow unbalanced', 1);

function run_contrasts(regressor_names, contrasts, SPMs, varargin)
    allow_unbalanced = 0;
    adding_contrasts = 1;
    removing_contrasts = 0;

    for i = 1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case {'allow_unbalanced' 'allow unbalanced'}
                    allow_unbalanced = varargin{i+1};
                case 'replace'
                    adding_contrasts = 0;
                    wh_cons = varargin{i+1};
                    if(length(contrasts) ~= length(wh_cons))
                        error('Number of contrasts (%d) does not match number to replace (%d) in SPM structure\n', length(contrasts), length(wh_cons));
                    end
                case 'remove'
                    removing_contrasts = 1;
            end
        end
    end

    regressor_names = cellstr(regressor_names);
    contrasts = cellstr(contrasts);
    if(ischar(SPMs))
        SPMs = cellstr(SPMs);
    end
    current_dir = pwd();


    try
        for i = 1:length(SPMs)
            if(isstruct(SPMs))
                SPM = SPMs(i);
                results_dir = SPM.swd;
                if(~isempty(results_dir))
                    cd(results_dir);
                end
            else
                results_dir = fileparts(SPMs{i});

                if(~isempty(results_dir))
                    cd(results_dir);
                end

                load('SPM');
            end

            check_spm_ver(SPM);

            expt_weights = compute_expt_weights(contrasts, regressor_names, SPM, allow_unbalanced);

            existing_num_contrasts = 0;
            if(removing_contrasts)
                fprintf('Removing previous contrasts\n');
                SPM.xCon = [];
            elseif(isfield(SPM, 'xCon'))
                existing_num_contrasts = length(SPM.xCon);
            end


            fprintf('Computing %d contrasts for %s\n', size(expt_weights, 1), SPM.swd);
            if(adding_contrasts)
                fprintf('Appending %d contrasts to end\n', size(expt_weights, 1));
                for j = 1:size(expt_weights, 1)
                    if(isempty(SPM.xCon))
                        SPM.xCon = spm_FcUtil('Set', contrasts{j}, 'T', 'c', expt_weights(j,:)', SPM.xX.xKXs);
                    else
                        SPM.xCon(end+1) = spm_FcUtil('Set', contrasts{j}, 'T', 'c', expt_weights(j,:)', SPM.xX.xKXs);
                    end
                end
                spm_contrasts(SPM, existing_num_contrasts+1:length(SPM.xCon));
            else
                fprintf('Replacing contrasts %s\n', num2str(wh_cons));
                for j = 1:size(expt_weights, 1)
                    SPM.xCon(wh_cons(j)) = spm_FcUtil('Set', contrasts{j}, 'T', 'c', expt_weights(j,:)', SPM.xX.xKXs);
                end

                spm_contrasts(SPM, wh_cons);
            end
            cd(current_dir);
        end
    catch
        cd(current_dir);
        error(lasterror());
    end
end

function check_spm_ver(SPM)
    spm_ver = spm('Ver');
    if(isempty(strfind(SPM.SPMid, spm_ver)))
        error('SPM.SPMid: "%s" does not match currently used SPM version: "%s"\n', SPM.SPMid, spm_ver);
    end
end

function expt_weights = compute_expt_weights(contrasts, regressor_names, SPM, allow_unbalanced)
    condition_weights = parse_contrast(contrasts, regressor_names);
    expt_weights = design_matrix_weights(regressor_names, condition_weights, SPM.xX.name);

    if(any(sum(abs(expt_weights), 2) == 0) || (~allow_unbalanced && (any(sum(expt_weights' ~= 0) == 0) || ~all(abs(sum(expt_weights')) < eps))))
        wh_bad = find((sum(expt_weights' ~= 0) == 0) | (abs(sum(expt_weights')) >= eps));
        for j = wh_bad
            fprintf('Contrast "%s" is empty or unbalanced\n', contrasts{j});
        end
        error('Some contrasts were empty or unbalanced');
    end
end