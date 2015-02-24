

DATA.TRIAL.onsets = [];
DATA.TRIAL.onsets{1} = o;
DATA.TRIAL.firpoints = k;
DATA.TRIAL.spersess = DATA.SPEC.spersess;


numsub=DATA.SPEC.numsub;
numreg = size(DATA.DATA.filtered_dat{1},2);

% ---------------------------------------------------------------------
%    * Design building
% ---------------------------------------------------------------------
    
% look for existing design matrices
if isfield(DATA.TRIAL,'DX'),
    fprintf(1,'Found existing design matrices in TRIAL.DX . ');
    dobuild = input('Re-create from onsets? (1/0) : ');
else
    dobuild = 1;
end

if ~dobuild
    % already got em, don't rebuild
    if length(DATA.DX) < numsub
        fprintf(1,'*** DX seems to short-replicating first matrix for each subject . \n');
        for s = 1:numsub, DATA.DX{s} = DATA.DX{1};, end
    end
 
else
 % build design matrices  
 
    if ~isfield(DATA.TRIAL,'onsets')
        error('Enter onset vectors in DATA.SPEC.onsets or use already-created DATA.DX cells');
    end
     
    % check for length of trigs (onsets)
    if length(DATA.TRIAL.onsets) < numsub
        fprintf(1,'*** Fewer onset vector cells than subjects - replicating onsets for each subject . \n');
        tmp = DATA.TRIAL.onsets;
        for s = 1:numsub, DATA.TRIAL.onsets{s} = tmp;, end
        clear tmp
    end
        
    fprintf(1,'generating design matrices . ');
    for s=1:numsub                                  % for all subjects
        DATA.TRIAL.DX{s} = build_fir_model(DATA.SPEC.onsets{s},DATA.TRIAL.firpoints,DATA.TRIAL.spersess);     
    end
    
end

% ---------------------------------------------------------------------
%    * Model Fitting
% ---------------------------------------------------------------------

k = DATA.TRIAL.firpoints;
fprintf(1,'Trial-level betas for subject...');

% for each subject   
for i = 1:numsub                     % size(c.data,3)  in matrix version  
    fprintf(1,'%3.0f ',i);
    
    DX = DATA.TRIAL.DX{i};
    o = DATA.TRIAL.onsets{i};
    
    for r = 1:numreg
        
        y = DATA.DATA.filtered_dat{i}(:,r);
        
        [b,X] = trial_level_beta(y,o{1},k,DX);
        
        DATA.TRIAL.betas{i}(:,r) = b;
    end
end

fprintf(1,'\n');


