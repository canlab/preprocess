function [cl_anova, EXPT] = clusters2anova(cwd,EXPT)
%
% Tor Wager
% 
% cl_summary is output of compile_clusters.m
%   use batch_compile_clusters
%
% EXPT is structure with DX and ANOVA info
%   use tor_setup_anova.m

if isempty(cwd), cwd = spm_get(-1,'*','Select dir with sub*_clusters.mat');,end

% set up filtering stuff
if ~isfield(EXPT,'FILT')
    EXPT.FILT.HP = EXPT.HP;
    EXPT.FILT.TR = EXPT.TR;
    EXPT.FILT.scanadjust = 1;
    EXPT.FILT.percent = 1;
    if EXPT.FILT.HP > 0, EXPT.FILT.doHP = 1;, else EXPT.FILT.doHP = 0;,end
    EXPT.FILT.filtertype = 'spm';
    EXPT.FILT.nruns = EXPT.nsess;
end

filtopt = EXPT.FILT;

% set up which model to use - same for all subjects, or different
if length(EXPT.DX.DX) == 1,
    subind = ones(size(EXPT.snums));
else
    subind = EXPT.snums;
end

% ---------------------------------------------------------------------
% * load all individual subjects and get HRF estimates for each
% ---------------------------------------------------------------------

fprintf(1,'\nLoading and fitting model - subject ')
D = dir;
subindex = 1;

for i = 1:length(D)
    if length(D(i).name) > 12
        if strcmp(D(i).name(1:3),'sub') & strcmp(D(i).name(end-11:end),'clusters.mat')
        
            fprintf(1,'%3.0f ',subindex)
            eval(['load ' D(i).name])
            
            for j = 1:length(clusters)
                
                % filter and adjust data
                % ---------------------------------------------------------------------
                filtopt.y = clusters(j).timeseries;
                y = filterAdjust(filtopt);
                
                % fit DX model to data, get timepoint estimate params
                % params: rows are deconv. parameter estimates, cols are subjects
                % ---------------------------------------------------------------------
                [params{j}(:,subindex), height, delay] = ts_stats(EXPT.DX.DX{subind(j)},y,EXPT.DX.paramsofinterest,EXPT.DX.baseparams);
                
            
            end
            
            subindex = subindex + 1;
            
                
        end % if sub*clusters.mat file   
    end % if namelen > 12
end % loop through dirs
    
    
% ---------------------------------------------------------------------
% * load all individual subjects and get HRF estimates for each
% ---------------------------------------------------------------------
fprintf(1,'\nRunning group anova - cluster ')

for i = 1:length(params)
    
    fprintf(1,'%3.0f ',i)
    [sig,T] = tor_mixed_anova(params{i},EXPT.ANOVA.factorcodes,EXPT.ANOVA.factornames);
    
    cl_anova(i).params = params{i};
    cl_anova(i).sig = sig;
    cl_anova(i).anovatable = T;
end


return
