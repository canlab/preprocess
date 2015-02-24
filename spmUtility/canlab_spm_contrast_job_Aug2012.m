function [matlabbatch, connames, contrast_vectors] = canlab_spm_contrast_job_Aug2012(modeldir, input_contrasts, varargin)
%
% [matlabbatch, connames, contrast_vectors] = canlab_spm_contrast_job(modeldir, input_contrasts, [optional inputs])
%
% Specify and run an SPM8 contrast manager job
%
% This function takes in the name of a directory [modeldir] containing an SPM analysis
%(i.e., SPM.mat file after a model has been estimated) and a set of
%contrasts you would like to create specified using condition names
%[input_contrasts] (see below).
%
% It does 6 things by default:
% 1. Construct contrast vectors and names across all sessions and add to
% matlabbatch SPM job structure
%
% 2. Scale contrast vectors according to the number of regressors/sessions
% included.
%
% 3. Check the contrast vectors to make sure they sum to zero
%
% 4. Save the matlabbatch job file in the SPM directory
%
% 5. Runs the job, which adds contrasts to the SPM.xCon field and creates
% con*imgs, t*imgs
%
% 6. Old (pre-existing) contrasts are removed by default!
%
% You can turn off some of these options using optional inputs:
% 'nodelete' will ADD contrasts to SPM.mat instead of removing old ones
% 'nosave' will avoid saving the job file
% 'noscale' will avoid scaling contrast weights
% 'norun' will avoid actually running the job (maybe you want to add it to
%  an existing workflow before running.)
% 
% 'exact' will match the names you provide EXACTLY with the contrast names.
%  This is necessary when working parametric modulators, to isolate the
%  condition from the modulator.  Be advised that the contrast names are
%  not exactly what you might expect - for example, a condition you named
%  'rate' might be named 'rate*bf(1)' for the purposes of this function
% 'anywhere(default)' will match the keyword (or some letters) you provide with the 
%  contrast names that have the keyword *anywhere' in it. This is useful if 
%  you want a quick way to include some common words or letters in the 
%  contrast names all together. 
% The default is to use 'anywhere' options.
% 
% If you want to use your customized contrast (e.g., -5 -3 -1 1 3 5), you can
% use the 'custom' option. see below for the example.
%
% NOTES on the structure of the input_contrasts
% -------------------------------------------------------------------------
% - The idea is that you know what your regressors are named in the SPM.mat
% file you already have.  You enter input_contrasts in a cell array, one
% cell per T-contrast, based on those names.
% - Each input contrast is itself a cell array with one or two cells.  The
% first cell specifies the names of contrasts with positive weights, and
% the second cell (if it exists) specifies the names of contrasts with
% negative weights.
% - Within each of the positive and negative cells, a cell vector contains
% the names of conditions you want to get positive or negative weights. It
% uses string matching, so any condition name containing the string you
% enter will be included.
% - Sorry this is a little complex, but it will work once you get the hang
% of it.
%
% Here is a way to specify an average across one regressor or a set of them
% that share a common name:
% input_contrasts{1} = {  {positive effect name} }; 
% 
% e.g., if you have two regressors called 'soc_friend_pre' and 'soc_friend_post',
% This contrast will specify the average across both regressors:
% input_contrasts{1} = {  {'soc_friend'} };
%
% This will do the same thing:
% input_contrasts{1} = {  {'soc_friend_pre' 'soc_friend_post'}   }
%
% You can use only some strings in the middle of regressor names. 
% e.g., input_contrasts{1} = {  {'friend'} };
%
% This will specify a DIFFERENCE pre - post, because the first cell has one positive condition (pre)
% and the second cell has one negative condition (post):
% input_contrasts{2} = {  {'soc_friend_pre'} {'soc_friend_post'}   };
%
% ...Since we have entered it as input_contrasts{2}, it will be added as a
% second contrast after the first one.
%
% The contrast below specifies an interaction, % (Rej - Friend) x (Pre - Post)
% input_contrasts{3} = [{{'soc_rejector_pre' 'soc_friend_post'}} {{'soc_rejector_pre' 'soc_friend_post'}}]
%
% ...Since we have entered it as input_contrasts{3}, it will be added as a third contrast.
%
% 'Custom' option example: 
% This is an example to use your customized contrast (e.g., 1st contrast: soc_friend_pre = -3, soc_friend_post = -1,
% soc_rejector_pre = 1, soc_rejector_post = 3, 2nd contrast: soc_friend_pre = 1, soc_friend_post = 2, soc_rejector_pre =3). 
%
%   input_contrasts{1} = { {'soc_friend_pre' 'soc_friend_post' 'soc_rejector_pre' 'soc_friend_post'} };
%   input_contrasts{2} = { {'soc_friend_pre' 'soc_friend_post' 'soc_rejector_pre'} };
%   canlab_spm_contrast_job(modeldir, input_contrasts, 'custom', {[-3 -1 1 3] [1 2 3]});
%
% -------------------------------------------------------------------------
%
% NOTES on scaling
% The default is to scale the contrast weights
% so that all contrast weights sum to one for pos and neg effects
% This means that contrasts can be interpreted as condition means
% and/or differences between condition means, across runs/sessions
%
% This is important for stats in the case where not all subjects
% have the same number of sessions, in which case the unscaled
% contrast values will be on different scales depending on the
% number of sessions, artifactually adding variance in the group analysis.
%
% Tor Wager and Wani Woo, Feb 2012
% - 'custom' option and string match function are upgraded (Wani Woo and Marieke Jepma), Aug 2012
%
% For a sample batch script that specifies and runs this for a number of
% subjects, search for canlab_spm_contrast_job on the lab WIKI.

deleteold = 1;
savejob = 1;
runjob = 1;
doscale = 1;
exact = 0;
anywhere = 1;
custom = 0;

[matlabbatch, connames, contrast_vectors] = deal({});

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'nodelete', deleteold = 0;
            case 'nosave', savejob = 0;
            case 'norun', runjob = 0;
            case 'noscale', doscale = 0;
            case 'exact', exact = 1; anywhere = 0;
            case 'anywhere', anywhere = 1; exact = 0; % default
            case 'custom', custom = 1; custom_contrasts = varargin{i+1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end



% SETUP dir and SPM.mat file
% -----------------------------------------------


spmmatname = fullfile(modeldir,'SPM.mat');

if ~exist(modeldir, 'dir') || ~exist(spmmatname, 'file')
    fprintf('Skipping: No directory or no SPM.mat: %s\n', spmmatname);
    return
end

fprintf('\nContrast Specification for \n%s\n---------------------------------------\n', spmmatname)

% SETUP names
% -----------------------------------------------
load(spmmatname)
names = SPM.xX.name;
for i = 1:length(names)
    names{i} = deblank(names{i}(6:end));
    names{i}(isspace(names{i})) = [];
end


% Setup matlabbatch and attach the SPM mat file
% ------------------------------------------------------------------
matlabbatch = {};

matlabbatch{1}.spm.stats.con.spmmat(1) = {spmmatname};

% Delete previous contrasts?
if deleteold
    
    matlabbatch{1}.spm.stats.con.delete = 1; % Yes = 1, No = 0;
    
else
    
    matlabbatch{1}.spm.stats.con.delete = 0; % Yes = 1, No = 0;
    
end


% for each input contrast...
% -----------------------------------------------
for i = 1:length(input_contrasts)
    
    input_con = input_contrasts{i};
    if custom, custom_con = custom_contrasts{i}; end
    
    if custom
        [mycon, conname] = get_customized_contrast_weights(names, input_con, custom_con, exact, anywhere);
    else
        [mycon, conname] = get_contrast_weights(names, input_con, exact, anywhere);
    end
    
    
    if doscale
        % Scale so that all contrast weights sum to one for pos and neg effects
        % This means that contrasts can be interpreted as condition means
        % and/or differences between condition means, across runs/sessions
        %
        % This is important for stats in the case where not all subjects
        % have the same number of sessions, in which case the unscaled
        % contrast values will be on different scales depending on the
        % number of sessions, artifactually adding variance in the group analysis.
        
        mycon(mycon > 0) = mycon(mycon > 0) ./ sum(mycon(mycon > 0));
        
        mycon(mycon < 0) = mycon(mycon < 0) ./ abs(sum(mycon(mycon < 0)));
    end

    
    if deleteold
        contrastnum = i;
    else
        contrastnum = i + length(SPM.xCon);
    end
    
    fprintf('Contrast %04d: %s %3.0f Pos weights, %3.0f Neg weights\n', contrastnum, conname, sum(mycon>0), sum(mycon<0));
    
    % Check
    if (any(mycon < 0) && any(mycon > 0)) && (sum(mycon) > 0.00000001)
        disp('WARNING!!!! CONTRAST WEIGHTS DO NOT SUM TO ZERO.  CHECK SPECIFICATION.');
    end
    
    
    % Save for output
    
    connames{i} = conname;
    contrast_vectors{i} = mycon;
    
    
    % Add to matlabbatch
    
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = conname; % contrast name
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.convec = mycon;
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
    
    
end


% save
if savejob
    
    savename = fullfile(modeldir,'spm_specify_contrasts_job.mat');
    save(savename, 'matlabbatch');
    fprintf('\nSaved batch job (matlabbatch) in spm_specify_contrasts_job.mat\n');
    
end

% run
if runjob
    
    spm('defaults', 'fmri');
    spm_jobman('initcfg'); % initialize
    spm_jobman('run', matlabbatch); % run
    
end


end % function




% NOTES

% % % attach name and contrast vector
% % % ------------------------------------------------------------------
% %
% % matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = conname; % contrast name
% % %matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Contrast Name2'; % contrast name
% % %matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Contrast Name3'; % contrast name
% %
% % matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = mycon;
% %
% % % fcon/tcon.sessrep : Replicate over sessions
% % % Don't Replicate: 'none'; Replicate: 'repl'; Replicate&Scale: 'replsc'; Creat per session: 'sess';
% % % Replicate + Create per session: 'both'; Replicate&Scale + Create per session: 'bothsc'
% %
% % %     sessrep.help    = {
% % %                    'If there are multiple sessions with identical conditions, one might want to specify
% % %                    contrasts which are identical over sessions. This can be done automatically based on
% % %                    the contrast spec for one session.'
% % %                    'Contrasts can be either replicated (thus testing average effects over sessions)
% % %                    or created per session. In both cases, zero padding up to the length of each session
% % %                    and the block effects is done automatically. In addition, weights of replicated
% % %                    contrasts can be scaled by the number of sessions. This allows to use the same contrast
% % %                    manager batch for fMRI analyses with a variable number of sessions. The scaled contrasts
% % %                    can then be compared in a 2nd level model without a need for further adjustment of effect
% % %                    sizes.'}';
% %
% % matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% %
% % %matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'repl';



function [mycon, conname] = get_contrast_weights(names, input_con, exact, anywhere)

if length(input_con) == 1, input_con{2} = []; end % Pos then neg conditions -- always 2

c = zeros(1, length(names));

mycon = c;

npos = length(input_con{1});
nneg = length(input_con{2});

for i = 1:npos
    
    if exact
        wh = find(strcmp(input_con{1}{i}, names) == 1);
        
    elseif anywhere
        
        wh = find(strncmp(input_con{1}{i}, names, length(input_con{1}{i})) == 1);
        
        if isempty(wh)
            wh = find(cellfun(@isempty, strfind(names, input_con{1}{i})) == 0);
        end
        
    end
    
    if isempty(wh)
        fprintf('WARNING: No matching regressors found for %s\n', input_con{1}{i});
    else
        mycon(wh) = 1;
    end
end

for i = 1:nneg
    
    if exact
        wh = find(strcmp(input_con{2}{i}, names) == 1);
        
    elseif anywhere
        
        wh = find(strncmp(input_con{2}{i}, names, length(input_con{2}{i})) == 1);
        
        if isempty(wh)
            wh = find(cellfun(@isempty, strfind(names, input_con{2}{i})) == 0);
        end
        
    end
    
    if isempty(wh)
        fprintf('WARNING: No matching regressors found for %s\n', input_con{2}{i});
    else
        mycon(wh) = -1;
    end
end

conname = cell2conname(input_con);

end % function


function [mycon, conname] = get_customized_contrast_weights(names, input_con, custom_con, exact, anywhere)
    
mycon = zeros(1, length(names));

for i = 1:length(custom_con)
    
    if exact
        wh = find(strcmp(input_con{1}{i}, names) == 1);
    elseif anywhere
        
        wh = find(strncmp(input_con{1}{i}, names, length(input_con{1}{i})) == 1);
        
        if isempty(wh)
            wh = find(cellfun(@isempty, strfind(names, input_con{1}{i})) == 0);
        end
        
    end
    
    if isempty(wh)
        fprintf('WARNING: No matching regressors found for %s\n', input_con{1}{i});
    else
        mycon(wh) = custom_con(i);
    end
end

conname = cell2conname_custom(input_con, custom_con);

end % function



function conname = cell2conname(input_con)

for i = 1:length(input_con{1}) - 1
    input_con{1}{i} = [input_con{1}{i} '+'];
end

for i = 1:length(input_con{2}) - 1
    input_con{2}{i} = [input_con{2}{i} '+'];
end

conname = ['(' cell2mat(input_con{1}) ')'];

if ~isempty(input_con{2})
    
    conname = [conname '-(' cell2mat(input_con{2}) ')'];
    
end


end % function



function conname = cell2conname_custom(input_con, custom_con)

for i = 1:length(input_con{1})
    input_con{1}{i} = [input_con{1}{i} ':' num2str(custom_con(i)) ', '];
end

conname = ['(' cell2mat(input_con{1}) ')'];

end % function


