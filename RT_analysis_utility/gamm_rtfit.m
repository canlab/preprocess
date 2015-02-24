function [phat,yfit,xx,yhist,rttmp,hand] = gamm_rtfit(rt,varargin)
% [phat,yfit,xx,yhist,rttmp,hand] = gamm_rtfit(rt,[fitopt],[phat],[plotopt],[color])
%
% phat: shape and scale parameters of gamma distribution
% yfit: fitted gamma values over xx intervals
% xx:   RT intervals for fitted values
% yhist: smoothed empirical PDF
% rttme: rt distribution with NaNs removed
% hand: handles to line objects for later use
%
% optional inputs:
% [fitopt]: fit option, 'full', '1st', or '2nd'
%   'full' fits both gamma parameters, 1st or 2nd fits only one
% [phat]: starting param estimates for gamma; required for '1st' or '2nd'
% [plotopt]: plotting options.
%       0 = no plot, 1 = plot only gamma on existing fig, >=2 = plot all
%
% tor wager, 10/10/05
%
% Examples:
%
% Basic fit with full plot
% [phat,yfit,xx,yhist,rttmp] = gamm_rtfit(rt);
%
% Full fit, plot only gamma fit in red
% [phat,yfit,xx,yhist,rttmp] = gamm_rtfit(rt,'full',[],1,'r');
%
% Fit only 2nd parameter, plot in green
% [phat,yfit,xx,yhist,rttmp] = gamm_rtfit(rt,'2nd',phat,1,'g');
%
% Fit only 1st parameer, Do not plot
% [phat,yfit,xx,yhist,rttmp] = gamm_rtfit(rt,'1st',phat,0);

fitopt = 'full'; phat = []; plotopt = 3; color = 'b';
hand = [];

if length(varargin) > 0, fitopt = varargin{1};, end
if length(varargin) > 1, phat = varargin{2};, end
if length(varargin) > 2, plotopt = varargin{3};, end
if length(varargin) > 3, color = varargin{4};, end

if plotopt > 1
    figure; set(gcf,'Color','w'); set(gca,'FontSize',16)
end

    % for histogram fitting
    % remove null values 
    rttmp = rt;
    rttmp(isnan(rttmp)) = []; 
    [hh,xx]=hist(rttmp,40); 
    
    % pad histogram with zeros
    diffxx = diff(xx);
    for i = 1:12, 
        xx = [xx(1)-diffxx(1) xx xx(end)+diffxx(1)];,
        hh = [0 hh 0];
    end
    
    % rescale to make a PDF
    hh = hh ./ sum(hh);  %.* length(hh) ./ sum(hh);
    
    if plotopt > 1
        h=bar(xx,hh);set(h,'FaceColor',[.8 .8 .8]);
    end
    
    
    
    % get smoothed empirical PDF
    
    yhist = [smooth_timeseries(hh,10)'];
    %xhist = [xx(1)-1 xx xx(end)+1];
    
    if plotopt > 1
        hold on; hand(1) = plot(xx,yhist,'k','LineWidth',3);
    end
    
    % fit gamma function
    switch fitopt
        case 'full'
            [phat] = gamfit(rttmp);
        case '1st'
            % we need an input phat w/one fixed param.  this is done above
            
            % construct a function handle with fixed parameter
            % this will evaluate the pdf with one fixed argument (phat(2)
            % is fixed)
            f = @(x,b) gampdf(x,b,phat(2));
            
            % max. likelihood est. of unknown parameter
            phat(1) = mle(rttmp,'pdf',f,'start',phat(1));
            
        case '2nd'
            % construct a function handle with fixed parameter
            % this will evaluate the pdf with one fixed argument
            f = @(x,b) gampdf(x,phat(1),b);
            
            % max. likelihood est. of unknown parameter
            phat(2) = mle(rttmp,'pdf',f,'start',phat(2));
            
        otherwise, 
            error('fitopt (1st var arg) must be ''full'' ''1st'' or ''2nd''');
    end
    
    
    %yfit = gampdf(1:max(rttmp),phat(1),phat(2)); 
    yfit = gampdf(xx,phat(1),phat(2));
    
    % scale gamma fit to be comparable to histogram
    %scaley = length(1:max(rttmp)) ./ length(xx);
    scaley = 1./sum(yfit);
    
    yfit = yfit.*scaley;
    
    if plotopt
        hold on; hand(2) = plot(xx,yfit,color,'LineWidth',3);
    
        if plotopt > 1
            set(gca,'FontSize',16)
            legend(hand,{'Empirical PDF' 'Gamma fit'},7);
            xlabel('Response Time (RT)'); ylabel('Probability');
        end
    end
    
    
    return
    
    