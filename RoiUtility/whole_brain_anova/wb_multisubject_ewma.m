function EXPT = wb_multisubject_ewma(EXPT,varargin)
% EXPT = wb_multisubject_ewma(EXPT,[type],[startat subject #],[start slice])
%
% Multi-subject shell for performing EWMA analysis on a group of subjects
% Main functions called:
% whole_brain_ewma
% (this is the single-subject EWMA analysis)
%
% start slice reverts to 1 after the first subject
%
% Examples:
% hewma   % start GUI, then use buttons
%
% Command line:
% EXPT = CreateExpt('hewma');   % setup
% % run EWMA, starting at subject 29, slice 1: ('full' is deprecated and
%   isn't used)
% EXPT = wb_multisubject_ewma(EXPT,'full',29);
%
% Start at subject 12, slice 10, then continue to Sub. 13, slice 1...n
% EXPT = wb_multisubject_ewma(EXPT,'full',12,10);


% start at subject and slice (2nd subject starts at slice 1 auto)
if length(varargin) > 0, type = varargin{1};  else  type = 'full';  end
if length(varargin) > 1, startat = varargin{2};  else  startat = 1;  end
if length(varargin) > 2, firstslice = varargin{3};  else  firstslice = 1;  end

runwh = startat:length(EXPT.subjects);
if length(varargin) > 3, runwh = varargin{4};   end

imgnames = EXPT.FILES.im_files;       % for full model

if isfield(EXPT,'mask')
    mask = spm_read_vols(spm_vol(EXPT.mask));
elseif isfield(EXPT.FIR,'mask')
    mask = spm_read_vols(spm_vol(EXPT.FIR.mask));
else
    mask = [];
end


for i = startat:length(EXPT.subjects)  % [6 7 8 9 27] %length(EXPT.subjects)
    
    %spmname = fullfile(EXPT.studydir{i},EXPT.subjects{i},'SPM.mat');
    %warning off; DX = spm2dx(spmname); warning on

    
    mkdir(EXPT.subjects{i})
    
    cd(EXPT.subjects{i})
    
    if exist(EXPT.FILES.im_files{1}(1,:)) == 2
        % only if we can find images
        
        %Pw = whole_brain_ewma(P,DX,TR,HP,nruns,numframes,[doplot],[mask],[smoothlen])
        Pw = whole_brain_ewma( ...
        imgnames{i}, ...
        [], ...      %(:,1:end-1), ...
        EXPT.TR, ...
        Inf, ...
        EXPT.FIR.nruns, ...
        EXPT.FIR.numframes, ...
        0, ...                      % do graphics
        mask, ...
        EXPT.FIR.smoothlen, ...
        firstslice, ...
        type ...
        );
    
        
    else
        disp(['Missing images for ' EXPT.subjects{i}])
    end
            
    firstslice = 1;     % re-set first slice to 1
    
    cd ..
    
end

return
