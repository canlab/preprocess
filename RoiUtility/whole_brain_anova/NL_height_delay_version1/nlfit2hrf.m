function [beta,betas] = nlfit2hrf(Y,TRin,varargin)
%
% [beta,betas] = nlfit2hrf(Y,TR,[plotfit])
% Y can have multiple column vectors for different HRFs
% - each column is fit separately
% 
global TR
TR = TRin;

%bc = beta2conditions(b,DX);

% -------------------------------------------------------------------
% * nonlinear fitting options - hardcoded
% -------------------------------------------------------------------
myoptions = optimset('LevenbergMarquardt','on');
startestimates = [1 6 0];	% height, time to peak, intercept (baseline h)
xdata = zeros(size(Y,1),1); xdata(1) = 1;
%xdata(end+1) = DX.TR;

%lb = [-Inf 3 -Inf];
%ub = [Inf 8 Inf];

if length(varargin) > 0
    figure;
end

% custom starting estimates
for i = 1:size(Y,2)
    y = Y(:,i);
    
                t = abs(y - mean(y)); s1 = y(t == max(t)) - mean(y); s1 = s1(end);
                s2 = find(t == max(t)); s2 = s2(end);
                s3 = mean(y);
            
        	    %[beta,resnorm,residual,exitflag,output,lambda,J]= ...
        	    %lsqcurvefit('nlhrf3',[s1 s2 s3],xdata,y,lb,ub,myoptions);
    
                [beta,r] = nlinfit(xdata,y,'nlhrf3',startestimates);
                
    		    %if exitflag < 0, beta = beta * Inf / Inf;, end
        	    betas(:,i) = beta;   %{m}(i(k),j(k),:) = beta;

     if length(varargin) > 0
         subplot(1,size(Y,2),i); hold on;
         plot(y,'k'); plot(nlhrf3(beta,xdata),'r');
     end
end
            
            
return
            