function [r,str,sig,myhandle] = plot_correlation(xvec,yvec,varargin)
% [r,infostring,sig,h] = plot_correlation(xvec,yvec,[textlabs],[color],[doquad],[dorobust])
% varargin is string of text labels
% This function makes a nice scatterplot of two variables with a regression line.  
% Options include colors, quadratic trend line, and robust estimation of regression
% line and correlation.  
%
% for text labels only, try:
% plot_correlation(beh1,mri1,highlow,'w.');
%
% color: 'bx', if empty or missing, creates new figure; otherwise uses existing current fig
% doquad: flag for quadratic correlations as well!
% dorobust: remove n outliers from data, using Min Cov Determinant (MCD) 
%          Rousseeuw, P.J. (1984), "Least Median of Squares Regression," 
%          Journal of the American Statistical Association, Vol. 79, pp. 871-88
%            outliers calculated using IRLS (robustfit.m) do not work well.
%           you enter n
%
% empty variable arguments are OK, defaults will be used
% color 'wo' plots open black circles

disp('OLD VERSION!!! SEE PLOT_CORRELATION_SAMEFIG FOR LATEST ROBUST CORREL!!!')

newfig = 0;
if length(varargin) > 0, mylabels = varargin{1};, else, mylabels = [];, end
if length(varargin) > 1, mycol = varargin{2};, else, mycol = 'ko';,newfig=1;,end
if length(varargin) > 2, doquad = varargin{3};, else, doquad = 0;, end
if length(varargin) > 3, dorobust = varargin{4};, else, dorobust = 0;, end
if isempty(mycol), mycol = 'ko';,newfig = 1;, end
if strcmp(mycol,'wo'), mycol = 'ko'; myfacec = 'w';,else,myfacec = mycol(1);,end

    if size(xvec,2) > size(xvec,1) & size(xvec,1)==1, xvec = xvec';, end
    if size(yvec,2) > size(yvec,1) & size(yvec,1)==1, yvec = yvec';, end
    
    X = [xvec ones(length(xvec),1)];
	X(isnan(xvec),:) = [];
      
	y = (yvec);
	y(isnan(xvec)) = [];
    
    b = X \ y;
    
    % robust regression - remove outliers
    if dorobust
        
        %if doquad
            %tmp = xvec .^2; tmp = tmp-mean(tmp);
            %[res]=fastmcd_noplot([tmp X(:,1) y]);
            %else
            %[res]=fastmcd_noplot([X(:,1) y]);
            [b,statsr]=robustfit(scale(X(:,1)),scale(y)); 
            nout = Inf; r= b(2); t = statsr.t(2);
            %end
        
        % remove n most extreme outliers and recompute correlation
        % OLD way - false positives!
        
        %wh = res.flag==0; nout = sum(res.flag==0);
        %y(wh) = []; X(wh,:) = []; xvec(wh) = []; yvec(wh) = [];
        
        %tmp = [(1:length(y))' stats.w]; tmp=sortrows(tmp,2); wh=tmp(1:dorobust,1);
        
        %keyboard
        %for i = 1:length(y),lab{i}=num2str(stats.w(i));,end
        %    makefigure(xvec,yvec,'wo',lab,doquad)
        %    hold on; plot(xvec(wh),yvec(wh),'ro')
        %end
        
        
        
    end
    
    
    
	%b = X \ y;
        
    myhandle = makefigure(X(:,1),y,mycol,mylabels,doquad,b,newfig,myfacec);

    
    % in regression, test significance of intercept parameter
    y2 = (y - X*b);   % subtract b1, the regression fit
    y2 = y2 + mean(y);  % add intercept back in
    try
        [h,b0p,ci,b0stats] = ttest(y2);
    catch   
        disp('No ttest.m: No stats toolbox?')
    end

    
    % separate correlation and t-test
    if ~dorobust
        try
            [h,p,ci,stats] = ttest(y);
        catch   
            disp('No ttest.m: No stats toolbox?')
        end
        
        r = corrcoef(y,X(:,1)); r= r(1,2);
        
        try
            [rci,sig,rZ,rp] = r2z(r,length(y),.05);
        catch   
            disp('Error in or missing r2z.m: No stats toolbox?')
        end
        
    else
        % robust irls
        %keyboard
            
            
    end
        
        text(min(X(:,1)),max(y),sprintf('r = %3.2f',r),'FontSize',16)
    
        try
            str = sprintf('Sig. of B0: u=%3.2f, t=%3.2f, p=%3.4f\n C: u=%3.2f, t=%3.2f, p=%3.4f   R: r=%3.2f, Z=%3.2f, p=%3.4f', ...
            mean(y2),b0stats.tstat,b0p,mean(y), stats.tstat, p, r, rZ, rp);
        catch
            str = ['Missing stats'];
            sig = NaN;
        end
    
        if dorobust, str=[str sprintf(' N_o_u_t=%3.0f',nout)];,end
        
    title(str,'FontSize',12)
    
    
    if doquad
        tmp = xvec .^ 2; tmp = tmp - mean(tmp);
        X = [tmp xvec ones(length(xvec),1)];
        X(isnan(xvec),:) = [];
	    y = (yvec);
	    y(isnan(xvec)) = [];
        
        [b,bint,dummy,rint,stats] = regress(y,X,.05);

        r0 = r;
        r = sqrt(stats(1));
        
        mysig = sum(sign(bint),2);
        for i = 1:size(b,1),if mysig(i),sigstr{i}='*';,else, sigstr{i}='';,end, end

        str = sprintf('y = %3.2fX^2%s + %3.2fX%s + %3.2f%s, r_0=%3.2f, Mult. r=%3.2f, Omni F=%3.2f, p=%3.2f', ...
            b(1),sigstr{1},b(2),sigstr{2},b(3),sigstr{3},r0, r, stats(2),stats(3));
        
        sig = stats(3) < .05;
        if dorobust, str=[str sprintf(' N_o_u_t=%3.0f',nout)];,end
       
        title(str,'FontSize',12)
        refcurve(b)

    end
    
    return
    
    
    
    
    
    function h = makefigure(xvec,yvec,mycol,mylabels,doquad,b,newfig,myfacec)

    if newfig,figure('Color','w'); 
        hold on; grid on; set(gca,'FontSize',18)
    end
   
        h = plot(xvec,yvec,mycol,'LineWidth',3,'MarkerSize',6,'MarkerFaceColor',myfacec);
        drawnow
        %x = min(xvec):max(xvec);
        %plot(x,b(1) * x + b(2),'k-','LineWidth',2)
       
    if doquad
        %refcurve(doquad)
    else
        try
            refline
        catch
            disp('No refline.m: no stats toolbox?')
            x = min(xvec):max(xvec);
            plot(x,b(1) * x + b(2),[mycol(1) '-'],'LineWidth',2)
        end
    end
        
	if ~isempty(mylabels)
		for j = 1:length(xvec)
			text(xvec(j),yvec(j),mylabels{j},'FontWeight','b')
		end
	end
    
    ylabel('Contrast beta value')
    xlabel('Behavioral score')
    
    return