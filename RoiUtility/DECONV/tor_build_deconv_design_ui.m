function EXPT = tor_build_deconv_design_ui(EXPT)
% EXPT = tor_build_deconv_design_ui(EXPT);
%
% Tor Wager
% 4/16/02
%
% Required fields:
% EXPT.evts_sf{i}   stick function of event onsets for each subject i
%                   created in read_ind_results.m 
%
% see also get_expt_info.m for menu of all design loading options
% see also group_fit_model, group_brain_nlfit, nl_contrasts 

disp('* Shape-free deconvolution model specification *')
disp('________________________________________________')
[num2str((1:length(EXPT.regnames))') repmat(' ',length(EXPT.regnames),1) str2mat(EXPT.regnames')]

if ~isfield(EXPT,'nsess')
        EXPT.nsess = input('How many runs (Sessions) in experiment? ');
end

EXPT.DX.nsess = EXPT.nsess;
EXPT.DX.regsofinterest = input('Enter vector of regressors of interest, not including baseline: ');
EXPT.DX.baseline = input('Enter number of baseline regressor, if any: ');
EXPT.DX.numframes = input('Enter number of time-points for each trial type to model (can be vector): ');
EXPT.DX.numframes = round(EXPT.DX.numframes);
if ~isfield(EXPT,'eres')
    EXPT.eres = input('Enter number of elements per second in Sess stick function (16 is default) ');
end
EXPT.DX.sfsamp = EXPT.eres; 
EXPT.DX.numbefore = input('Enter number of timepoints to model before start of trial (0 is start at start of trial):');
EXPT.DX.allsubssame = input('Do all subjects have exactly the same design matrix? (1/0) ');

if EXPT.DX.allsubssame, mylist = 1;, else, mylist = EXPT.snums;,end
EXPT.DX.sepsessions = input('Model conditions separately in different sessions? (1/0) ');

EXPT.DX.dxnames = EXPT.regnames(EXPT.DX.regsofinterest);
if ~isempty(EXPT.DX.baseline), EXPT.DX.dxnames{end+1} = 'Baseline';,end
EXPT.DX.dxnames{end+1} = 'Intercept';


for i = mylist		% my subjects...or 1 for all Ss the same.
    
    mysf = EXPT.evts_sf{i}(:,EXPT.DX.regsofinterest);

    % ------------------------------------------------------------------------------
    % * separate sessions stuff
    % ------------------------------------------------------------------------------
    if EXPT.DX.sepsessions, 
        [mysf,mydxnames] = sep_sessions(mysf,EXPT);, 
        if length(EXPT.DX.numframes) > 1
            EXPT.DX.numframes = repmat(EXPT.DX.numframes,1,EXPT.nsess);
        end
    end

    [EXPT.DX.DX{i}] = tor_make_deconv_mtx3(mysf,EXPT.DX.numframes,EXPT.DX.sfsamp,EXPT.DX.numbefore,EXPT.nsess);

    if ~isempty(EXPT.DX.baseline)
        % ------------------------------------------------------------------------------
        % * baseline is modeled as a single regressor, so that params may be adjusted relative to the baseline param.
        % ------------------------------------------------------------------------------
        EXPT.DX.baseframes = input(['Enter duration (time-points) of baseline points to include in baseline (subject ' num2str(i) '): ']);
        mysf = EXPT.evts_sf{i}(:,EXPT.DX.baseline);
        mydx = tor_make_deconv_mtx3(mysf,EXPT.DX.baseframes,EXPT.DX.sfsamp);
        mydx = sum(mydx(:,1:end-1),2);
        EXPT.DX.DX{i} = [EXPT.DX.DX{i}(:,1:end-EXPT.nsess) mydx EXPT.DX.DX{i}(:,end-EXPT.nsess+1:end)];
    end
    
    if length(EXPT.DX.numframes) == 1
        EXPT.DX.dxtrialonsets = 1:EXPT.DX.numframes:size(EXPT.DX.DX{i},2)-(EXPT.nsess+length(EXPT.DX.baseline));  % should be same for all
    else
        a=[1 cumsum(EXPT.DX.numframes)];
        a(2:end) = a(2:end)+1;
        a(end) = [];
        EXPT.DX.dxtrialonsets = a;
    end
    
    EXPT.DX.cf{i} = dx2cf(EXPT.DX.DX{i},EXPT.DX.dxtrialonsets);
    if EXPT.DX.sepsessions,
        EXPT.DX.baseparams = size(EXPT.DX.DX{i},2) - EXPT.nsess;
    else
        EXPT.DX.baseparams = size(EXPT.DX.DX{i},2) - 1;
    end
end

if EXPT.DX.sepsessions, EXPT.DX.dxnames = mydxnames;, end

return


function [cf] = dx2cf(DX,onsets)
% ------------------------------------------------------------------------------
% * make condition function for regs of interest
% ------------------------------------------------------------------------------
cf = zeros(size(DX,1),1);
%mytp = diff(onsets);
%if any(mytp - mean(mytp)), error('Not all trial types have same number of time points.'),end
%mytp = mytp(1);

for i = 1:length(onsets)
    %thison = onsets(i); 
    %whichel = find(sum(DX(:,thison:thison+mytp-1),2));
    whichel = find(DX(:,onsets(i)));
    cf(whichel) = i;
end
return



function [sfout,nmout] = sep_sessions(mysf,EXPT)
% ------------------------------------------------------------------------------
% * make separate columns for regs of interest in each session
% ------------------------------------------------------------------------------

nsess = EXPT.nsess;
spss = EXPT.nscans * EXPT.DX.sfsamp / nsess;
if round(spss) ~= spss, error('Num scans not equally divisible by num sessions.'), end
myz = zeros(spss, length(EXPT.DX.regsofinterest));
nmout = [];	

for i = 1:nsess
	ses_sf{i} = mysf((i-1)*spss+1:i*spss,:);
	sfout{i} = [repmat(myz,i-1,1); ses_sf{i}; repmat(myz,nsess-i)];
	
	mynm = EXPT.DX.dxnames(1:length(EXPT.DX.regsofinterest));	% already regs of int + intercept
	for j = 1:length(mynm), mynm{j} = [mynm{j} '_s' num2str(i)];, end
	nmout = [nmout mynm];
end

sfout = cell2mat(sfout);
nmout = [nmout EXPT.DX.dxnames(end)];


return

