function [stop_phat,c,X,Y] = create_SSRT_surface(SSmu,SSsig);
% x = 1:100; y = 0; p = .1:.1:.9; a = .5 ; b = .5; m = .7;
%
% Example: From Gonzales and Wu
% s1 = create_utility_surface(x,20,p,a,b,m);
% x = [25 50 75 100 150 200 400 800 50 75 100 150 150 200 200];
% y = [0  0  0  0   0   0   0   0   25 50 50  50  100 100 150];
% p = [.01 .05 .1 .25 .4 .5 .6 .75 .9 .95 .99];
% a = .6 ; b = .6; m = .8;
% figure;
% [u,s1] = create_utility_surface(x,y,p,a,b,m);
%
% % random a
%for i = 1:100
%a = 2.*rand;
%[u(:,:,i)] = create_utility_surface(x,y,p,a,b,m,0);
%end
 

 
    s1 = [];

nsamples = 50000;
A = normrnd(450,75,nsamples,1);
B = exprnd(10,nsamples,1);
responses = A+B; rt{1} = responses;

SOA = 200;

for i=1:length(SSmu),
    for k=1:length(SSsig),
        stop_phat(i,k) = Calc_diff_exgauss(SSmu(i),SSsig(k),0,SOA,rt{1});
    end
end


    [X,Y] = meshgrid(SSmu,SSsig);
    tmpX = X'; tmpY = Y';
    %surf(X,Y,stop_phat);
    %set(s1,'EdgeColor','none','FaceAlpha',.5)
    
    figure;
    v = [.2 .3 .4 .5 .6 .7 .8];
    [c,h] = contourf(tmpX,tmpY,stop_phat,v);
    %[C,H] = contour('v6',X,Y,stop_phat);
    
    colorbar
    colormap(gray)
    set(gcf,'Color',[1 1 1]);
    %hold on; wh = find(stop_phat > .29 & stop_phat < .31);
    %yfit = [smooth_timeseries(tmpY(wh),5)'];
    %hold on; plot(tmpX(wh),yfit,'w-','LineWidth',6)
    %plot(X(wh),Y(wh),'w-','LineWidth',3)
 
    title('Stop Bias')
    xlabel('Stop Latency Mu')
    ylabel('Stop Latency Variance')
    zlabel('False Alarm Rate')
 

 


