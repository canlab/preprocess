function [u,s1] = create_utility_surface(x,y,p,a,b,m,lo,doplot)
% x = 1:100; y = 0; p = .1:.1:.9; a = .5 ; b = .5; m = .7;
%
% Example: From Gonzales and Wu
% s1 = create_utility_surface(x,20,p,a,b,m);
% x = [25 50 75 100 150 200 400 800 50 75 100 150 150 200 200];
% y = [0  0  0  0   0   0   0   0   25 50 50  50  100 100 150];
% p = [.01 .05 .1 .25 .4 .5 .6 .75 .9 .95 .99];
% a = .6 ; b = .6; m = .8; lo = 1.2;
% figure;
% [u,s1] = create_utility_surface(x,y,p,a,b,m,lo,1);
%lo
% % random a
%for i = 1:100
%a = 2.*rand;
%[u(:,:,i)] = create_utility_surface(x,y,p,a,b,m,0);
%end

if nargin < 7, doplot = 1; end

    s1 = [];
    
ind = 1;
for pp = p
    %u(:,ind) = prospect_forward(x',y',pp,a,b,m);
    u(:,ind) = prospect_utility([a b m lo],x',y',pp);
    
    ind = ind + 1;
end

if doplot
    [X,Y] = meshgrid(p,x);

    s1 = surf(X,Y,u);
    set(s1,'EdgeColor','none','FaceAlpha',.5)


    xlabel('probability')
    ylabel('value')
    zlabel('utility')

end

end
