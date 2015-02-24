function [DATA,SPEC] = cluster_mvroi(EXPT,cl)
% [DATA,SPEC] = function cluster_mvroi(EXPT,cl)
%
% Purpose: Take an EXPT and clusters structure and translate data 
% and parameters into the correct format for mvroi, the scnlab
% multivariate timeseries region-of-interest analysis tool.
%
% Load the EXPT and clusters files
% Make a subdirectory for the analysis and go there
% specs are saved in mvroi_specs

% remove just from here to save memory
DATA=struct([]);
if isfield(cl,'raw_data'), cl = rmfield(cl,'raw_data');, end
if isfield(cl,'INDIV'), cl = rmfield(cl,'INDIV');, end
if isfield(cl,'imP'), cl = rmfield(cl,'imP');, end
if isfield(cl,'prin_comps'), cl = rmfield(cl,'prin_comps');, end


% Options

SPEC = [];

% load existing spec file if possible
if exist('mvroi_specs.mat') == 2, disp('Loading mvroi_specs.mat'), load mvroi_specs, end

% names of regions
if isfield(cl,'shorttitle'),
    disp('Loading names of clusters from cl.shorttitle')
    for i = 1:length(cl), 
        SPEC.names{i} = cl(i).shorttitle;, 
    
        if isempty(SPEC.names{i})
            try
                cluster_orthviews(cl(i),{[1 0 0]});
                SPEC.names{i} = input('Enter string name for this cluster: ','s');
            catch
                disp('Error viewing cluster for naming');
            end
        end
    end
else
    disp('Looked for names of clusters in cl.shorttitle,but didn''t find any.')    
end


if isfield(EXPT,'TR'), SPEC.TR = EXPT.TR;, end
if isfield(EXPT,'HP'), SPEC.HP = EXPT.HP;, end 
    

if isfield(EXPT,'cov'), 
    if ~isempty(EXPT.cov)
        disp('Behavior: saving EXPT.cov column 1 in SPEC.beh')
        SPEC.beh = EXPT.cov(:,1);, 
    end
elseif isfield(EXPT,'behavior'), 
    if ~isempty(EXPT.behavior)
        disp('Behavior: aving EXPT.behavior column 1 in SPEC.beh')
        SPEC.beh = EXPT.behavior(:,1);,
    end
end


if isfield(EXPT,'DX')
    if isfield(EXPT.DX,'nsess'), SPEC.nruns = EXPT.DX.nsess;, 
        fprintf(1,'SPEC.nruns should contain number of runs.  Found in EXPT.DX.nsess: %s\n',num2str(SPEC.nruns)),
    end  
    
    if isfield(EXPT.DX,'nruns'), SPEC.spersess = EXPT.DX.nruns;, 
            fprintf(1,'SPEC.spersess should contain # of imgs in each run.  Found in EXPT.DX.nruns: %s\n',num2str(SPEC.spersess)),
    end
    if isfield(EXPT.DX,'numframes'), SPEC.firpoints = EXPT.DX.numframes;, end
    
    if isfield(EXPT.DX,'TR'), SPEC.TR = EXPT.DX.TR;, end
    if isfield(EXPT.DX,'HP'), SPEC.HP = EXPT.DX.HP;, end    
    
    if isfield(EXPT.DX,'DX'), SPEC.DX = EXPT.DX.DX;,DATA.DX = SPEC.DX;, end
    
    if isfield(EXPT.DX,'model'), SPEC.DX = EXPT.DX.model;,DATA.DX = SPEC.DX;, end
end


if isfield(EXPT,'FIR')
    if isfield(EXPT.FIR,'nsess'), SPEC.nruns = EXPT.FIR.nsess;, end  
    
    if isfield(EXPT.FIR,'nruns'), SPEC.spersess = EXPT.FIR.nruns;, end
    if isfield(EXPT.FIR,'numframes'), SPEC.firpoints = EXPT.FIR.numframes;, end
    
    if isfield(EXPT.FIR,'TR'), SPEC.TR = EXPT.FIR.TR;, end
    if isfield(EXPT.FIR,'HP'), SPEC.HP = EXPT.FIR.HP;, end    
    
    if isfield(EXPT.FIR,'DX'), SPEC.DX = EXPT.FIR.DX;,DATA.DX = SPEC.DX;, end

    if isfield(EXPT.FIR,'model'), SPEC.DX = EXPT.FIR.model;,DATA.DX = SPEC.DX;, end

end

if ~isfield(SPEC,'nruns') & isfield(SPEC,'spersess'), SPEC.nruns = length(SPEC.spersess);, end

    
% save  spec file 
disp('Saving mvroi_specs.mat'), save mvroi_specs SPEC


% Data

% define dat

doindiv = 0;
if isfield(cl(1),'indiv_timeseries'), 
    doindiv = input('indiv_timeseries found.  Run on individual peaks? (1/0) ');
elseif ~isfield(cl(1),'all_data'), 
        disp('Must have indiv_timeseries or all_data field in clusters.');
        return
end

if doindiv,
    for i = 1:length(EXPT.subjects),
        for j = 1:length(cl)
            dat{i}(:,j) = cl(j).indiv_timeseries(:,i);
        end
    end
    
else
    if size(cl(1).all_data,1) == length(EXPT.subjects)
        fprintf(1,'all_data field should be time x subjects for each region.\n')
        fprintf(1,'you have what seems to be individual contrast values .\n')
        fprintf(1,'Consider running cluster_nmdsfig instead.\n')
        return
    end
    
    for i = 1:length(EXPT.subjects),
        for j = 1:length(cl)
            dat{i}(:,j) = cl(j).all_data(:,i);
        end
    end
end



% check for not-OK subjects
if isfield(cl,'sbjctOK'),
    if any(~cl(1).sbjctOK), 
        wh = find(~cl(1).sbjctOK);
        disp(['Found bad subjects: ' num2str(wh) ' . Removing. ']);
        if isfield(SPEC,'beh'),
            SPEC.beh(wh) = [];
        end
        
        if isfield(SPEC,'DX'),
            SPEC.DX(wh) = [];
        end
        if isfield(DATA,'DX'),
            DATA.DX(wh) = [];
        end
        
        dat(wh) = [];
    end
end



% check for size match between models and data
% may not match if some subs are missing some sessions

if isfield (DATA,'DX'),if isempty(DATA.DX), DATA = rmfield(DATA,'DX');,end,end

if isfield(DATA,'DX'),
    for i = 1:length(dat),
        if size(DATA.DX{i},1) ~= size(dat{i},1)
              disp(['Model size and data size do not match for subj: ' num2str(i) ' . Truncating dat and states (if entered). ']);
              dat{i} = dat{i}(1:size(DATA.DX{i},1),:);
              
              % check for states
              if isfield(SPEC,'states'),
                  SPEC.states{i} = SPEC.states{i}(1:size(DATA.DX{i},1),:);
              end
        end
    end
end
              
              
% check data for missing values (NaNs) and all-0 voxels
allwh = [];
for i = 1:length(dat),
    wh = find(isnan(sum(dat{i})) | std(dat{i}) < eps*10);
    if ~isempty(wh)
        disp(['S ' num2str(i) ' has empty regions: ' num2str(wh)]);
        allwh = [allwh wh];
    end
end
if ~isempty(allwh)
    disp('Regions with missing data for any subject will be eliminated.');
    for i = 1:length(dat)
        dat{i}(:,allwh) = [];
    end
    
    if isfield(SPEC,'names'), SPEC.names(allwh) = [];,end
end

if isempty(DATA)
    DATA=struct('DATA',struct('dat',dat));
else
    DATA.DATA.dat = dat;
end

fprintf(1,'Attached data in DATA.DATA.dat.\n')
dosave = input('Save DATA.mat? (1/0) ');
if dosave, DATA.SPEC = SPEC;, save DATA DATA, end


fprintf(1,'You may still have to build a SPEC.states field, if you have multiple conditions.\n')

go = input('Try DATA = mvroi(DATA,''mvroi_specs'');   to run.  Run now? (1/0) ');

if go, 
    DATA = mvroi(DATA,'mvroi_specs');
end


return



    