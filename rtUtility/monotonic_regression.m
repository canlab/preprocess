function [npos,df,rx,ry,i] = monotonic_regression(x,y)
%
% tor wager
% June 14, 2006
%
% Is the null hypothesis 1.5?
% Simulation: test proportion positive against 0.5, which it should be
% x = randn(20,5000);
% y = randn(20,5000);
% [npos,df,rx,ry,i] = monotonic_regression(x,y);
% dev = npos ./ (df) - .5;
% mean(dev)./ste(dev)
%
% Is the distribution binomial? (no...)
% y = binopdf(unique(npos),19,.5);
% figure;plot(y,'k');
% clear y2
% u = unique(npos); for i=1:length(u),y2(i) =sum(npos==u(i));, end
% y2 = y2./sum(y2);
% plot(u,y2,'ro-');

i = zeros(size(x)); rx = i; ry = i;
npos = zeros(1,size(x,2));
df = size(x,1) - 1;

for n = 1:size(x,2)
    
    rx(:,n) = rankdata(x(:,n)); 
    ry(:,n) = rankdata(y(:,n));
    [sortedx,i(:,n)] = sort(rx(:,n));
    npos(n) = sum(diff(  ry(i(:,n),n)  ) > 0);
    
    %plot(rx(i),ry(i),'ko-');
end



return
