function block_phats = interpolate_exgauss(rt,phatgo)
%
% block_phats = interpolate_exgauss(rt,phatgo)
%
% function used in concertion with expected earnings function; 
% interpolates between two points, by sampling from the two distributions
% at percentiles (.10)
%
% for example, at .10 between distribution 1 and 2, we will contatenate
% distribution 1 (9 times) and distribution 2(1 time) and then fit the rt
% data to exgaussian parameters mu,sigma, and tau.
% 
% input arguments: rt(cell array of distributions) and phatgo (starting estimates
% of exgaussian fit)

npoints = 21;    % number of interpolation points between distributions
block_phats = zeros(length(npoints),3);


for i=1:9,
    dist = repmat(rt{:,1},10-i,1);
    dist = [dist;repmat(rt{:,2},i,1)];
    block_phats(i+1,:) = exgauss_rtfit(dist,'full',phatgo,0);   
end

for i=1:9,
    dist = repmat(rt{:,2},10-i,1);
    dist = [dist;repmat(rt{:,3},i,1)];
    block_phats(i+11,:) = exgauss_rtfit(dist,'full',phatgo,0);
end


return