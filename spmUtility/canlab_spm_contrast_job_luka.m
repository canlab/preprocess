function [matlabbatch, connames, contrast_vectors] = canlab_spm_contrast_job_luka(modeldir, input_contrasts, varargin)
%
% [matlabbatch, connames, contrast_vectors] = canlab_spm_contrast_job_luka(modeldir, input_contrasts, [options])
%
% Specify and run an SPM8 contrast manager job in modeldir with input_contrasts.
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
% con_*.imgs, spmT_*.imgs
%
% 6. Old (pre-existing) contrasts are removed by default!
%
% You can turn off some of these options using optional inputs:
% 'nodelete' will ADD contrasts instead of removing old ones
% 'addcons' will ONLY ADD contrasts that don't already exist
% 'nosave' will avoid saving the job file
% 'noscale' will avoid scaling contrast weights
% 'norun' will avoid actually running the job (maybe you want to add it to
%  an existing workflow before running.)
% 
% Matching regressor names:
% 'Sn(?) ' will be stripped from the beginning of regressor names (as seen
%   in SPM.xX.name) unless regexp matching is selected.
% DEFAULT: will match the keyword (or some letters) you provide with the 
%	contrast names that have the keyword *anywhere* in it. This is useful if 
%	you want a quick way to include some common words or letters in the 
%	contrast names all together.
% 'exact' will match the names you provide EXACTLY with the contrast names.
%	This is necessary when working parametric modulators, to isolate the
%	condition from the modulator.  Be advised that the contrast names are
%	not exactly what you might expect - for example, a condition you named
%   'rate' might be named 'rate*bf(1)' for the purposes of this function
% 'regexp' match to regressor names using regular expressions (see help regexp)
%   '^' will be prepended to expressions
%     'heat' will become '^heat' and only match regressor names starting with 'heat'
%     '.*heat' will become '^.*heat' and match regressors names containing 'heat'
%	CAREFUL! regular expressions easily go awry
% 'defaultsuffix', sfx
%   will add sfx to the end of each provided keyword
%   only used in regexp matching mode
%   ex: 'defaultsuffix', '\*bf\(1\)'
%       keyword 'rate' will become regexp '^rate.* \*bf\(1\)$'
%
% 'names', customnames
%	use custom names (strings in the cell array customnames)
%	(empty or undefined cells get DEFAULT: built-in naming scheme)
%	EX: names{3} = 'task_quad'; % input_contrast{3} gets name 'quad'
% 'weights', customweights
%	use custom weights (arrays in cell array customweights)
%	(empty or undefined fields get DEFAULT weights [1 -1])
%	EX: weights{3} = [3 1 1 3];
%	  input_contrast{3} will test quadratic effect across 4 (sets of) regressors
%	  (note: input_contrasts{3} should have 4 cells)
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
% Custom weights example: 
% This is an example to use your customized contrast (e.g., 1st contrast: soc_friend_pre = -3, soc_friend_post = -1,
% soc_rejector_pre = 1, soc_rejector_post = 3, 2nd contrast: soc_friend_pre = 1, soc_friend_post = 2, soc_rejector_pre =3).  
%   input_contrasts{3} = { {'soc_friend_pre'} {'soc_friend_post'} {'soc_rejector_pre'} {'soc_friend_post'} };
%   customweights{3} = [-3 -1 1 3];
%   input_contrasts{4} = { {'soc_friend_pre'} {'soc_friend_post'} {'soc_rejector_pre'} };
%   customweights{4} = [1 2 3];
%   canlab_spm_contrast_job_luka(modeldir, input_contrasts, 'weights', customweights);
% 
% Custom name example:
%   input_contrasts{3} = { {'soc_friend_pre'} {'soc_friend_post'} {'soc_rejector_pre'} {'soc_friend_post'} };
%   customweights{3} = [-3 -1 1 3];
%   customnames{3} = 'soc_linear';
%   canlab_spm_contrast_job_luka(modeldir, input_contrasts, 'names', customnames, 'weights', customweights);
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
% - 'regexp' option, 'addcons' options, minor edits (Luka Ruzic, Sept 2012)
%
% For a sample batch script that specifies and runs this for a number of
% subjects, search for canlab_spm_contrast_job on the lab WIKI.

deleteold = 1;
savejob = 1;
runjob = 1;
doscale = 1;
addcons = 0;
matching = 'anywhere';
default_suffix = '';
for i=1:numel(input_contrasts), weights{i} = [1 -1]; end  %#ok
custom_connames = cell(numel(input_contrasts));

[matlabbatch, connames, contrast_vectors] = deal({});

i=1;
while i<=length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'nodelete', deleteold = 0;
            case 'nosave', savejob = 0;
            case 'norun', runjob = 0;
            case 'noscale', doscale = 0;
            case 'addcons', addcons = 1; deleteold = 0;
            case 'exact', matching = 'exact';
            case 'regexp', matching = 'regexp';
            case 'suffix', i=i+1; default_suffix = varargin{i};
            case 'weights'
                i=i+1;
                for j=1:numel(varargin{i})
                    if ~isempty(varargin{i}{j})
                        weights{j} = varargin{i}{j}; %#ok
                    end
                end
            case 'names'
                i=i+1;
                for j=1:numel(varargin{i})
                    custom_connames{j} = varargin{i}{j};
                end
            otherwise, error(['Unknown input string option: ' varargin{i}]);
        end
    else
        disp(varargin{i})
        error('Above input argument is unrecognized')
    end
    i=i+1;
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
if ~strcmp(matching,'regexp')        
    names = regexprep(names,'^\S* ','');
end

if addcons
    existing_connames = {};
    for i = 1:numel(SPM.xCon)
        existing_connames{i} = SPM.xCon(i).name; %#ok
    end
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
c=0;
for i = 1:length(input_contrasts)
    if (numel(input_contrasts{i}) > numel(weights{i}))
        error('Problem with contrast %d: number of elements (%d) is greater than the number of weights (%d).',i,numel(input_contrasts{i}),numel(weights{i}))
    end
        
    [mycon, conname] = get_contrast_weights(names, input_contrasts{i}, custom_connames{i}, weights{i}, matching, default_suffix);
    
    if addcons
        if numel(existing_connames)>=i & isempty(strmatch(conname,existing_connames{i},'exact'))
            warning('input contrast %d (%s) does not match existing contrast %d (%s).',i,conname,i,existing_connames{i}) %#ok
        end
        
        if ~isempty(strmatch(conname,existing_connames,'exact'))
            fprintf('Contrast %s exists: skipping.\n',conname)
            continue
        end
    end
    c=c+1;
    
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
        contrastnum = c;
    else
        contrastnum = c + length(SPM.xCon);
    end
    fprintf('Contrast %04d: %s %3.0f Pos weights, %3.0f Neg weights\n', contrastnum, conname, sum(mycon>0), sum(mycon<0));
    
    % Check
    if (any(mycon < 0) && any(mycon > 0)) && (sum(mycon) > 0.00000001)
        disp('WARNING!!!! CONTRAST WEIGHTS DO NOT SUM TO ZERO.  CHECK SPECIFICATION.');
    end
        
    % Save for output
    contrast_vectors{c} = mycon;    
    
    % Add to matlabbatch    
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = conname;  %#ok  % contrast name
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = mycon; %#ok
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; %#ok
end


if c == 0
    fprintf('\nNo contrasts to run: exiting.\n')
    return
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



function [mycon, conname] = get_contrast_weights(names, input_con, custom_conname, weights, matching, default_suffix)

mycon = zeros(1, length(names));

for i = 1:numel(input_con)
    for j = 1:numel(input_con{i})
        wh = [];  %#ok
        switch matching
            case 'anywhere'
                wh = find(strncmp(input_con{i}{j}, names, length(input_con{i}{j})) == 1);
                if isempty(wh)
                    wh = find(cellfun(@isempty, strfind(names, input_con{i}{j})) == 0);
                end
            case 'exact'
                wh = strmatch(input_con{i}{j}, names, 'exact');
            case 'regexp'
                expr = input_con{i}{j};
                %                 expr = ['^' input_con{i}{j}];
                if ~isempty(default_suffix) && ~strcmp(expr(end),'$')
                    expr = [expr '.*' default_suffix '$']; %#ok
                end                
                wh = find(~cellfun('isempty',regexp(names, expr)));
            otherwise
                error(['PROGRAMMER''s ERROR: matching code is unrecognized: ' matching])
        end
        
        if isempty(wh)
            error(['No matching regressors found for ' input_con{i}{j}]);
        end
        
        mycon(wh) = weights(i);
    end
end       

if ~isempty(custom_conname)
    conname = custom_conname;
else
    if strcmp(mat2str(weights),mat2str([1 -1]))
        conname = cell2conname(input_con);
    else
        conname = cell2conname_custom(input_con, weights);
    end
end

end % function


function conname = cell2conname(input_con)

for i = 1:length(input_con{1}) - 1
    input_con{1}{i} = [input_con{1}{i} '+'];
end
conname = ['(' cell2mat(input_con{1}) ')'];

if length(input_con) == 2    
    for i = 1:length(input_con{2}) - 1
        input_con{2}{i} = [input_con{2}{i} '+'];
    end
    conname = [conname '-(' cell2mat(input_con{2}) ')'];    
end

end % function



function conname = cell2conname_custom(input_con, custom_con)

for i = 1:length(input_con)
    for j = 1:length(input_con{i})
        if ~exist('conname','var')
            % initialize
            conname = sprintf('(%s:%s', input_con{i}{j}, num2str(custom_con(i)));
        else
            conname = sprintf('%s, %s:%s', conname, input_con{i}{j}, num2str(custom_con(i)));
        end
    end
end

conname = [conname ')'];

end % function


