function [m,fit,x,b,bf] = fit_basis_set(y,TR,robust,varargin)
%[m,fit,x,b,bf] = fit_basis_set(y,TR)
% tor wager
%
% inputs
% y = data, 
% robust = 1 or 0 for robust fit
% optional delta function (matrix form)
% optional vector of number of basis groups (shifted) to model each column
% of delta function
%
% outputs
% -------
% returns fitted value, and 
% max amplitude (positive), time to peak, and half-height to half-height
% time
% model matrix x
%
% [m,fit,x] = fit_basis_set(roi.adjustedy,TR,delta);
%figure;plot(y);hold on; plot(fit,'r')
%[m,fit,x,b] = fit_basis_set(y,TR);
%figure;plot(y);hold on; plot(fit,'r')
%
% roi.adjustedy = trimts(roi.adjustedy,3,[]);
%roi.adjustedy = trimts(roi.adjustedy,3,[]);
% b1 =V2(1:size(DX,2),1:size(DX,2)) * pinv(DX) * roi.adjustedy;
%tmp3 = reshape(b1(1:end-1),60,6);
% figure('Color','w');subplot(1,3,1); plot(tmp3);
%legend({'1' '2' '5' '6' '10' '11'})
%[m,fit,x] = fit_basis_set(roi.adjustedy,TR,0,delta);
%subplot(1,3,2);hold on; plot(fit)
%legend({'1' '2' '5' '6' '10' '11'})
%[m,fit,x] = fit_basis_set(roi.adjustedy,TR,1,delta);
%subplot(1,3,3);hold on; plot(fit)


if length(varargin) == 0    
    % we have a delta matrix
    len = length(y);
else
    len = round(30 ./ TR);
end

if length(varargin) > 1
    xnum = varargin{1};
else
    xnum = ones(size(varargin{1},2));
end

% build basis functions

type = 'td8';

switch type
    case 'td8'
        
% HRF+time+dispersion, shifted by 8 s
xBF.dt = TR;
xBF.length = 30;
xBF.name = 'hrf (with time and dispersion derivatives)'; xBF = spm_get_bf(xBF);
ba = xBF.bf(1:len,:);
bb = [zeros(round(8./TR),size(ba,2)); ba];
bb = bb(1:len,:);
x = [ba bb];
bf = x;

case 'gamma4'
    dt = TR;
   u     = [0:(32/dt)] - 0/dt;
   hrf   = spm_Gpdf(u,6,dt)';
   hrf = [hrf; zeros(round(10./TR),1)];
   ind = 1;
   for i = round((0:4:12) ./ TR)
       h2 = [zeros(i,1); hrf]; 
       x(:,ind) = h2(1:len);
       ind = ind + 1;
   end
   x = [gradient(x(:,1)),x];
   bf = x;
   
otherwise, error('unknown option.')
end
    
    
% convolve, if necessary

if length(varargin) > 0    
    % we have a delta matrix
    ind = 1;
    d = varargin{1};
    for i = 1:size(d,2)
        for j = 1:size(x,2)
            tmp = conv(x(:,j),d(:,i));
            tmp = tmp(1:length(y));
            xx(:,ind) = tmp;
            ind = ind+1;
        end
    end
    x = xx;
else
    d = 1;
end

x = scale(x,1); x(:,end) = 1;

% model fitting and parameterization
if robust
    b = robustfit(x,y);
else
    x(:,end+1) = 1;
    b = pinv(x) * y;
end


for j = 1:size(d,2)
    
    wh = [(j-1)*size(bf,2)+1:j*size(bf,2)]; % which basis f/betas to use to get fit for this condition
    fit(:,j) = bf * b(wh);
    
    m(1,j) = max(fit(:,j));
    tmp = find(fit(:,j) > .5*m(1,j)); tmp2 = tmp(end); tmp = tmp(1);
    m(2,j) = find(fit(:,j)==m(1,j)) .* TR;
    m(3,j) = (tmp2 - tmp) .* TR;
end


return



