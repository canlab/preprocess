function [b,n,bavg,navg,f,bste,OUT] = nl_hrf_plot(coords,EXPT,plottype,varargin)
% [b,n,bavg,navg,f,bste] = nl_hrf_plot(coords,EXPT,plottype,[opt] whichcon)
% Tor Wager 7/25/02
%
% coords are x,y,z list of mm (world space) coordinates to extract and average over
% in 3-row or 3-vector point list
%
% Plottype:
%   1   grand averages only
%   2   grand avgs plus individual subjects
%   3   grand averages plus contrasts
%
% whichcon: vector of length = # of contrasts, specifying 1 or 0 for each
%   1 plots contrast, 0 omits
%
% contrast HRFs are based on the deconvolution estimates, and are
% adjusted by removing the intercept (from nonlinear fitting) before plotting.
% contrast SE bars are based on the std. error of the difference
% where individual differences are computed as avg(+) - avg(-).
%

% Variables:
% Pb    cell array of dx_beta imgs - all HRF estimates concatenated 
% Pf    cell array of cond*imgs - height, delay, intercept estimates for each 
%       experimental condition
% b     matrix of betas (from dx_beta*img) for each subject.  
%       Subjects are rows, betas are columns
% n     cell array of length 3 (height, delay, intercept)
%       each cell contains parameter estimates averaged over voxels listed in coords
% bc    average of betas
% bcon  contrast HRF estimates for each contrast(cells).  {contrast}{1} is +, {2} is -
% n     nonlinear param estimates for each parameter/condition/subject, 
%       params are cells, subjects are rows, conditions are columns
% navg  average n across subjects
% indncon
%       contrast estimates for each subject(rows)/contrast(cells)
% nconse    
%       within Ss contrast error for each contrast (cells)

if length(varargin) > 0, 
    wcon = varargin{1};, 
else 
    wcon = ones(1,size(EXPT.DX.contrasts,1));
end

fprintf(1,'\nnl_hrf_plot.m -> ')

% ----------------------------------------------------------------------------
% define name of nonlinear fitting function to use
% ----------------------------------------------------------------------------

funcname = 'nlhrf3'; %['nlhrf2_TR' num2str(EXPT.TR)];  
global TR
TR = EXPT.TR;

% ----------------------------------------------------------------------------
% get names of beta img files and nlfit imgs for each condition for each subject
% ----------------------------------------------------------------------------

Pb = tor_list_files(EXPT.subjects,'dx_beta*img');

for p = 1:length(EXPT.DX.params)
    
    % save img names for each subject for this parameter
    
    for i = 1:length(EXPT.subjects)
    
        ncond = size(EXPT.DX.condP{1},1);
        Pf{i}{p} = [repmat(EXPT.subjects{i},ncond,1) repmat(filesep,ncond,1) EXPT.DX.condP{p}];

    end
    
end

% ----------------------------------------------------------------------------
% loop through subjects and get betas and nonlinear fit params
% ----------------------------------------------------------------------------

for i = 1:length(EXPT.subjects)
    
    fprintf(1,'%s ',EXPT.subjects{i})
    
    % extract betas from this subject
    b(i,:) = indiv_get_nlparams(coords,Pb{i});
    
    % extract nlfit params
    % n: rows are ss, cols are conditions, cells are params
    for p = 1:length(EXPT.DX.params)
        n{p}(i,:) = indiv_get_nlparams(coords,Pf{i}{p});
    end
    
end

fprintf(1,'\n ')

% ----------------------------------------------------------------------------
% average betas across ss and break into conditions
% ----------------------------------------------------------------------------

bavg = nanmean(b);
bste = nanstd(b) ./ sqrt(size(b,1));

bc = beta2conditions(bavg,EXPT.DX); % cell array
bcste = beta2conditions(bste,EXPT.DX); % cell array
bavg = bc;
clear bste

for i = 1:size(EXPT.DX.contrasts,1)
    bcon{i}{1} = mean(cell2mat(bc(EXPT.DX.contrasts(i,:) > 0)'));
    bcon{i}{2} = mean(cell2mat(bc(EXPT.DX.contrasts(i,:) < 0)'));
    
    %bste{i}{1} = mean(cell2mat(bcste(EXPT.DX.contrasts(i,:) > 0)'));
    %bste{i}{2} = mean(cell2mat(bcste(EXPT.DX.contrasts(i,:) < 0)'));
    %bste{1}(i) = mean(nanstd(n{1}) * EXPT.DX.contrasts(i,:)');
    %bste{1}(i) = 
end

% ----------------------------------------------------------------------------
% average nlfit params across ss
% ----------------------------------------------------------------------------

for p = 1:length(EXPT.DX.params),navg{p} = nanmean(n{p});, end


% ----------------------------------------------------------------------------
% calculate contrasts across conditions and standard errors of differences
% ----------------------------------------------------------------------------
for i = 1:size(EXPT.DX.contrasts,1)
    if wcon(i)        
    % individual contrast values and s.e. of difference
        indncon{i} = mean(n{1} .* repmat(EXPT.DX.contrasts(i,:),size(n{1},1),1),2);
        nconse{i} = nanstd(indncon{i}) ./ sqrt(size(indncon{i},1));
        
        try
            ncon{i}(1) = mean(navg{1}(EXPT.DX.contrasts(i,:) > 0) .* EXPT.DX.contrasts(i,EXPT.DX.contrasts(i,:) > 0))';
            ncintercept{i}(1) = mean(navg{3}(EXPT.DX.contrasts(i,:) > 0) .* EXPT.DX.contrasts(i,EXPT.DX.contrasts(i,:) > 0))';
            
        catch
        end
        try
            ncon{i}(2) = mean(navg{1}(EXPT.DX.contrasts(i,:) < 0) .* EXPT.DX.contrasts(i,EXPT.DX.contrasts(i,:) > 0))';
            ncintercept{i}(2) = mean(navg{3}(EXPT.DX.contrasts(i,:) < 0) .* EXPT.DX.contrasts(i,EXPT.DX.contrasts(i,:) > 0))';
        catch
        end
    end
 end
    

% ----------------------------------------------------------------------------
% plot betas
% ----------------------------------------------------------------------------

[f,h] = plot_dx_and_nlfit(bavg,navg,EXPT.DX.numframes,funcname);
subplot 131; title(['Est. HRF Grand Avgs'],'FontSize',14)
subplot 132; title(['Fitted Gamma Functions'],'FontSize',14)
legend(h,EXPT.DX.dxnames)

subplot 133; 
nste{1} = nanstd(n{1}) ./ sqrt(size(n{1},1));
bar(navg{1})
tor_bar_steplot(navg{1},nste{1},{'b'})
set(gca,'XTickLabel',EXPT.DX.dxnames)
title(['Parameter Estimates'],'FontSize',14)


if plottype == 2
    for i = 1:length(EXPT.subjects)
        for j = 1:3, n2{j} = n{j}(i,:);,end
        plot_dx_and_nlfit(beta2conditions(b(i,:),EXPT.DX),n2,EXPT.DX.numframes,funcname);
        subplot 121; title(['S' num2str(i)],'FontSize',14)
    end
end

if plottype == 3
    for i = 1:size(EXPT.DX.contrasts,1)
    if wcon(i)
        figure;hold on; grid on; set(gcf,'Color','w')
        subplot 121, hold on
        
        %plot((1:length(bcon{i}{1})).*TR,bcon{i}{1}-mean(bcon{i}{1}),'r','LineWidth',2)
        %plot((1:length(bcon{i}{2})).*TR,bcon{i}{2}-mean(bcon{i}{2}),'b','LineWidth',2)
        
        % Adjust for intercept, and plot
        if ~isnan(bcon{i}{1})
            plot((1:length(bcon{i}{1})).*TR,bcon{i}{1}-ncintercept{i}(1),'r','LineWidth',2)
        end
        if ~isnan(bcon{i}{2})
            plot((1:length(bcon{i}{2})).*TR,bcon{i}{2}-ncintercept{i}(2),'b','LineWidth',2)
        end
        title([EXPT.DX.connames{i}])
        legend({'Positive' 'Negative'})
        
        % se of fitted height difference (contrast)
        %tor_bar_steplot(bcon{i}{1},bste{i}{1},{'r'})
        %tor_bar_steplot(bcon{i}{2},bste{i}{2},{'b'})
        
        subplot 122, hold on
        
        bar(ncon{i})
        tor_bar_steplot(ncon{i},repmat(nconse{i},1,length(ncon{i})),{'b'});
        set(gca,'XTick',1:2)
        set(gca,'XTickLabel',{'+' '-'})
        axis auto %set(gca,'YLim',get(gca,'YLim') .* 1.5)
        title(['Contrast Param Estimates'],'FontSize',14)

        OUT.indncon = indncon;
        OUT.nconse = nconse;
    end
    end   
end

OUT.Pb = Pb;
OUT.Pf = Pf;
OUT.b = b;
OUT.bavg = bavg;
%OUT.bste = bste;
OUT.n = n;
OUT.navg = navg;
OUT.nste = nste;


return





function [b] = indiv_get_nlparams(coords,P)

% transform coords to VOXEL space from MM (world) space
V = spm_vol(P(1,:)); 
V.M = V.mat;

O.coords = mm2voxel(coords,V);

warning off
% extract betas from dx_beta_imgs or cond*imgs and average across voxels
%b = timeseries2('multi',P,O);
b = timeseries3(O.coords,P);
b = b.avg';
warning on
    
if ~any(b), warning('Empty timeseries!'), keyboard, end

return




function [f,h] = plot_dx_and_nlfit(bavg,navg,numframes,funcname)

global TR

figure;set(gcf,'Color','w'),colordef('white');
mycol = {'r-' 'g-' 'b-' 'm-' 'c-' 'k-' 'y-' 'b-'};

h1 = subplot(1,3,1); hold on; grid on;
for i = 1:length(bavg)
    plot((1:length(bavg{i})).*TR,bavg{i},mycol{i},'LineWidth',2)
end

xlabel(str2mat(mycol)')

h2 = subplot(1,3,2);
hold on; grid on;
for i = 1:length(navg{1})   % for each condition
    
    % define input delta function and params for this condition/subject
    xdata = zeros(1,numframes); xdata(1) = 1;
    for j = 1:length(navg)
        params(j) = navg{j}(i);
    end
    
    % construct fit
    f{i} = eval([funcname '(params,xdata);']);
    
    % plot fit
    h(i) = plot((1:length(f{i})).*TR,f{i},mycol{i}(1),'LineWidth',2);
    
    % plot individual points
    %for j = 1:size(n{1},1) % j is subject, i = condition, {} are params h,del,intercept
    %        plot(n{2}(j,i)+6,n{1}(j,i)+n{3}(j,i),[mycol{i}(1) 'o'],'LineWidth',2)
    %end
    
end

equalize_axes([h1 h2])


return
