function EXPT = tor_setup_anova_ui(EXPT)

% get DX matrix
% fit DX matrix
% code DX betas by condition
% fit ANOVA
    % error must be based on subject * condition interaction
% save each m.e. and interaction p-value in map

% [P,T,STATS,TERMS]=anovan(Y,GROUP,MODEL,SSTYPE,GNAME,DISPLAYOPT)
% Subject is fixed effect
% MSB / MSE = F, so replace MSE with sub * condition interaction

factor_index = 2;
EXPT.ANOVA.numfactors = input('Enter number of factors in ANOVA, not including subject: ');
EXPT.ANOVA.factorcodes = [];
EXPT.ANOVA.factorcols = [];

% ----------------------------------------------------------------------------
% * set up subject factor
% ----------------------------------------------------------------------------
EXPT.ANOVA.factornames = {'subject'};
EXPT.ANOVA.numsubparams = EXPT.DX.numframes * length(EXPT.DX.regsofinterest);
myones = ones(EXPT.ANOVA.numsubparams,1);
EXPT.ANOVA.factorcodes{1} = myones;
for i = 2:length(EXPT.snums)
    EXPT.ANOVA.factorcodes{1} = [EXPT.ANOVA.factorcodes{1}; myones .* EXPT.snums(i)];
end

% ----------------------------------------------------------------------------
% * set up time factor
% ----------------------------------------------------------------------------

timefactor = input('Include time within trial as a factor? (1/0): ');
if timefactor
    tf = (1:EXPT.DX.numframes)';
    tf = repmat(tf,size(EXPT.ANOVA.factorcodes{1},1) ./ EXPT.DX.numframes,1);
    EXPT.ANOVA.factorcodes{factor_index} = tf;
    EXPT.ANOVA.factornames{factor_index} = 'Time';
    factor_index = factor_index + 1;
end

% ----------------------------------------------------------------------------
% * set up condition-related factors
% ----------------------------------------------------------------------------

% make cf for DX betas of interest - what will be cells in ANOVA
% for group ANOVA with subject as random effect
% --------------------------------------------------------------
factorcf = zeros(EXPT.ANOVA.numsubparams,1);
cindex = 1;
for j = 1:EXPT.DX.numframes:size(EXPT.DX.DX{1},2) - 2
    factorcf(j:j + EXPT.DX.numframes - 1) = cindex;
    cindex = cindex + 1;
end
EXPT.ANOVA.factorcf = repmat(factorcf,length(EXPT.snums),1);

disp('Trial types are:')
[num2str((1:length(EXPT.DX.dxnames))') repmat(' ',length(EXPT.DX.dxnames),1) str2mat(EXPT.DX.dxnames')]
disp(['Enter trial-type factors, or return to finish with this step.'])
myname = 'Factor1';

while ~isempty(myname)
    
    % get group membership of each trial type
    % --------------------------------------------------------------
    myname = input('Enter name for this factor: ','s');
    if isempty(myname),break,end
    EXPT.ANOVA.factornames{factor_index} = myname;
    grp = [];
    while length(grp) < length(EXPT.DX.regsofinterest)
        grp = input(['Enter vector of group membership for each trial type for factor ']);
        if ~(length(grp) == length(EXPT.DX.regsofinterest))
            disp(['Length of grp vector should be the same as index values in cond. function, which is ' num2str(length(EXPT.DX.regsofinterest))])
        else
            EXPT.ANOVA.factorcodes{factor_index} = zeros(size(EXPT.ANOVA.factorcodes{1}));
        end
    end
  
    % make revised condition function for the factor based on group
    % --------------------------------------------------------------
    for j = 1:max(EXPT.ANOVA.factorcf)      % j is the original trial type
        whichel = find(EXPT.ANOVA.factorcf == j);
        EXPT.ANOVA.factorcodes{factor_index}(whichel) = grp(j);
        
        % for within-subject ANOVA
        %whichel = find(EXPT.DX.cf{1} == grp(j));
        %EXPT.ANOVA.factorcodes{factor_index}(whichel) = j;
    end
    EXPT.ANOVA.factorcols{factor_index} = grp;
    factor_index = factor_index + 1;
    
end


%[P,T,STATS,TERMS]=anovan(Y,GROUP,'full',3,gname,'off')