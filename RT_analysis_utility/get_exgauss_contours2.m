% relevant output:
% pgo       predicted correct go trials based on gamma dist
% pstop     predicted stop accuracy based on gamma dist
% timeouts  predicted "too slow" trials 

%%% this version different from get_exgauss_contours.m due to changes in
%%% the fitted contour lines; attempting to plot the fitted EI and SSRT
%%% functions


clear pgo
clear pstop



    ssoffset = RTCUTOFF - 200;

    % ENTER varystop as the EI/SSRT measure
    % ENTER varyslow = phatmix(mu1,mu2,mu3)
    % ENTER thetavals in go-no-stop order (rows) by mu-sigma-tau (columns)
    % ENTER SSRTy = predicted stops for go-no-stop blocks
    % ENTER SSRTx = timeouts
    
    %varystop = 50:50:300;         %  rows, stop slowness param
    %varyslow = phatgo(1)./2:25:phat(1)+300;     %  cols, go slowdown param
    
    %thetavals = [phatmix(1,1),phatmix(1,2),phatmix(1,3); phatmix(2,1),phatmix(2,2),phatmix(2,3); phatmix(3,1), phatmix(3,2),phatmix(3,3)];
    %thetavals = sortrows(thetavals,1);
    
  
    
    %muvals = [linspace(mu1-(mu2-mu1),mu2,50), linspace(mu1,mu2,50), linspace(tau2,tau3,50), linspace(tau3,tau3+(tau3-tau2),50)];
    %sigvals = [linspace(sig1-(sig2-sig1),sig2,50),linspace(sig1,sig2,50),linspace(sig2,sig3,50),linspace(sig3,sig3+(sig3-sig2),50)];
    %tauvals = [linspace(tau1-(tau2-tau1),tau2,50),linspace(tau1,tau2,50),linspace(tau2,tau3,50),linspace(tau3,tau3+(tau3-tau2),50)];
    % these replace varyslow.
    
    
    %if isempty(varyslow), varyslow = phat(2):10:phat(2)+100; ,end
    varyx = 1:max(rt)+100;      %  variable to integrate over
    
    % initialize output
    pgo = zeros(length(varystop),length(varyslow));
    pstop = zeros(size(pgo));
    
    
    ssr = 1;

        for r = 1:length(varyslow)  % muvals
            
            %mu = phat(1);
            
            mu = thetavals(r,1);
            sig = thetavals(r,2);
            tau = thetavals(r,3);
           
            
            theta = [mu sig tau];
            % GO (FAST ENOUGH) probability
            pgo(ssr,r) = Iexgauss(RTCUTOFF,theta);
            %yy = gampdf(0:max(rt),phat(1),r);
            %hold on; plot(yy,'k');
        
        
            % stop probability      % integrate numerically
            % varslow = slowing parameter
            % varystop = stop speed parameter
            %pstop(ssr,r) = estimate_pstop(rt,theta,varystop(ssr),ssoffset);  
            
            
            pstop(ssr,r) = Calc_diff_exgauss2(varystop(ssr),phatgo(2),phatgo(3),mu,sig,tau,ssoffset);
            
        end                 % loop r (GO slowdown)
        
  
    
    b = glmfit(SSRTx,SSRTy);
    ssrtfit = b(1) + b(2)*SSRTx;
    
    stopnames = {'EI','SSRT'};
    %fprintf(1,'\n')
    
    %figure;plot(pgo',pstop');
    
    
     tor_fig;
     plot(SSRTx,pstop); hold on;
     set(gca,'YLim',[0 1],'XLim',[0 1]);
     % plot([timeouts(1) timeouts(1)],get(gca,'YLim'),'k');  % black line
     
     legstr = {};
     for i = 1:length(stopnames), legstr{i} = stopnames{i};,end
     
     legend(legstr)
     xlabel('Proportion Time-outs'); ylabel('Proportion correct stops');
     
     h3(1) = plot(truetimeout,truestop,'ko','MarkerSize',8,'MarkerFaceColor',[.9 .9 .9]);
     plot(SSRTx,SSRTy,'g');
     
     
     

     
     