% relevant output:
% pgo       predicted correct go trials based on gamma dist
% pstop     predicted stop accuracy based on gamma dist
% timeouts  predicted "too slow" trials 


clear pgo
clear pstop



    ssoffset = RTCUTOFF - 200;

    varystop = 1:10:60;         %  rows, stop slowness param
    varyslow = phatgo(2):10:phat(2)+100;     %  cols, go slowdown param
    if isempty(varyslow), varyslow = phat(2):10:phat(2)+100; ,end
    varyx = 1:max(rt)+100;      %  variable to integrate over
    
    % initialize output
    pgo = zeros(length(varystop),length(varyslow));
    pstop = zeros(size(pgo));
    
    
    for ssr = 1:length(varystop)

        for r = 1:length(varyslow)

            % GO (FAST ENOUGH) probability
            pgo(ssr,r) = gamcdf(RTCUTOFF,phat(1),varyslow(r));
            %yy = gampdf(0:max(rt),phat(1),r);
            %hold on; plot(yy,'k');
        
        
            % stop probability      % integrate numerically
            pstop(ssr,r) = estimate_pstop(rt,[phat(1) varyslow(r)],varystop(ssr),ssoffset);     
            
        end                 % loop r (GO slowdown)
        
    end
    %fprintf(1,'\n')
    
    %figure;plot(pgo',pstop');
    timeouts = 1 - pgo;
    
     tor_fig;
     plot(timeouts',pstop'); hold on;
     set(gca,'YLim',[0 1],'XLim',[0 1]);
     plot([timeouts(1) timeouts(1)],get(gca,'YLim'),'k');
     
     legstr = {};
     for i = 1:length(varystop), legstr{i} = num2str(varystop(i));,end
     
     legend(legstr)
     xlabel('Proportion Time-outs'); ylabel('Proportion correct stops');
     
     

     
     