function [DATA,betas,bmean] = mvroi_fit_model(DATA,varargin);
% [DATA,betas,bmean] = mvroi_fit_model(DATA,[plot]);
% Extracts betas given DATA.DX and DATA.dat for each condition in a
% {time x region} cell matrix, 1 cell per subject array in DATA.dat
%
% DATA.SPEC.numframes contains length of beta series (in TRs) 
%
% [opt last input:  plot -- default is 1, yes]
%
% Inputs:
% Requires DATA.DATA.filtered_dat field
%
% Outputs:
% DATA.DATA.resids, residuals of filtered data
% DATA.DATA.fits, timeseries of model fits for each region
% DATA.DATA.b, vector of FIR estimates concatenated
% DATA.DATA.fir is FIR HRF estimates for each event type; cells are {regions x event types},
% elements are (subjects x time)
% 

doplot = 1; f1 = [];
if length(varargin) > 0, doplot = varargin{1};,end

numsub=length(DATA.DATA.filtered_dat);    %number of subjects
numreg=size(DATA.DATA.filtered_dat{1},2);  %number of ROIs


fprintf(1,'FIR for subject...');

cumframes=[1 cumsum(DATA.SPEC.firpoints)+1]; 

% for each subject   
for i = 1:numsub                     % size(c.data,3)  in matrix version  
    
    fprintf(1,'%3.0f ',i);
     m = DATA.DATA.filtered_dat{i};
     b = pinv(DATA.DX{i}) * m;      % betas
     f = DATA.DX{i} * b;            % fit
     r = m - f;                     % residuals 
     
     % loop through event types, and get FIR betas for each event
     for n = 1:length(cumframes)-1      
         whichbetas=cumframes(n):cumframes(n+1)-1;
         
         % loop through regions
         for region = 1:numreg             
            betas{region,n}(i,:) = b(whichbetas,region)';
         end
         
     end   
     
     
    DATA.DATA.resids{i} = r;
    DATA.DATA.fits{i} = f;
    DATA.DATA.b{i} = b;
     
end

DATA.DATA.fir = betas;
    
fprintf(1,'\n');

% get means across subjects
for i = 1:size(betas,1),
    for j = 1:size(betas,2)
        bmean{j}(i,:) = mean(betas{i,j});
    end
end


if doplot
    [f1,colors,axh] = tor_fig(1,length(DATA.SPEC.firpoints));
    for n=1:length(DATA.SPEC.firpoints)
        subplot(1,length(DATA.SPEC.firpoints),n);
        plotdata=bmean{n};
        plot(plotdata','linewidth',2);
        axis auto
    end    
end


return
