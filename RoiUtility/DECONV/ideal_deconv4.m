function [id1,b1,id2,b2] = ideal_deconv4(DX,TR,dosel,varargin)
% function [id1,b1,id2,b2] = ideal_deconv4(DX,TR,dosel,[opt] conditions2convolve)
% tests deconvolution matrix directly against idealized data
% you put in the exact temporal sequence to be deconvolved,
% in the form of the DX matrix.
%
% Tor Wager, 4/19/02
%
% inputs:
%	DX	deconvolution matrix
%	tp	time points estimated for each condition in DX
%	TR	repetition time of scan
%	dosel	1 or 0, test selective average method also
%   contrast    row vector of contrast weights for conditions
%       e.g., [1 3 5]

if length(varargin) > 0
    cond2c = varargin{1};
else
    cond2c = zeros(size(DX,2)); % too long, but we don't know how many conditions...
    cond2c(1) = 1;
end

mycolors = {'ro-' 'b^--' 'go:' 'y^--' 'mo-' 'r^-' 'g^-' 'b^-' 'k^-' 'm^-' 'y^-'};

[delta, wd, wb] = DX_find_delta(DX);
delta = delta(:,wb > 0);
disp(['Conditions detected: ' num2str(size(delta,2))])

hrf = spm_hrf(TR);
hrf = hrf ./ max(hrf);


% plot dimensions - number of plots
ncols = 4;
if dosel, nrows = 2;, else, nrows = 1;, end


% -------------------------------------------------------------------
% * set up response to conditions specified in cond2c
% -------------------------------------------------------------------

%ideal_data = conv(hrf,delta(:,1));
%ideal_data = ideal_data(1:size(DX,1));

ideal_data = [];
for i = 1:size(delta,2)
    if cond2c(i) ~= 0
        mydata = conv(hrf .* cond2c(i),delta(:,i));
        mydata = mydata(1:size(DX,1));
        ideal_data(:,i) = mydata;
    end
end
ideal_data = sum(ideal_data,2);

% autocorrelated noise, SNR = 1
myvar = var(ideal_data);
load('myscannerxc.mat')
noise = noisevector(length(ideal_data),myscannerxc,myvar)';
%noise = rand(size(ideal_data));
%noise = noise - mean(noise);
ideal_data = ideal_data + noise;

b = pinv(DX) * ideal_data;

id1 = ideal_data;
b1 = b;

% -------------------------------------------------------------------
% * plot single response
% -------------------------------------------------------------------
figure;
set(gcf,'Color','w')

subplot(nrows,ncols,1)
plot(hrf,'LineWidth',2)
hold on; 
%plot(b(1:end-1),'rs-')

for i = 1:length(wd)
    if i > length(mycolors),mycolors{i} = mycolors{i-7};,end
    plot(b(wd(i):wd(i)+wb(i)),mycolors{i})
    myleg{i} = ['Condition' num2str(i)];
end
grid on
title(['Contrast ' num2str(cond2c) ' SNR = 1'])

% -------------------------------------------------------------------
% * set up multiple response
% -------------------------------------------------------------------

ideal_data = [];
for i = 1:size(delta,2)
        mydata = conv(hrf,delta(:,i));
        mydata = mydata(1:size(DX,1));
        ideal_data(:,i) = mydata;
end
ideal_data = sum(ideal_data,2);

b = pinv(DX) * ideal_data;

id2 = ideal_data;
b2 = b;

% -------------------------------------------------------------------
% * plot multiple response
% -------------------------------------------------------------------
plot_multiple(2,hrf,b,wd,wb,DX,mycolors,nrows,ncols,'Ideal response, all conditions')

% -------------------------------------------------------------------
% * set up noisy response
% -------------------------------------------------------------------

noise = 1 * rand(size(ideal_data));
myvar = var(ideal_data);
noise = noisevector(length(ideal_data),myscannerxc,myvar)';
%noise = noise - mean(noise);
%ideal_data = ideal_data + noise;
b = pinv(DX) * ideal_data;

id3 = ideal_data;
b3 = b;

% -------------------------------------------------------------------
% * plot noisy response
% -------------------------------------------------------------------
plot_multiple(3,hrf,b,wd,wb,DX,mycolors,nrows,ncols,'Noisy signal, all conditions')


% -------------------------------------------------------------------
% * set up random response
% -------------------------------------------------------------------

ideal_data = rand(size(ideal_data));
b = pinv(DX) * ideal_data;

id4 = ideal_data;
b4 = b;

% -------------------------------------------------------------------
% * plot random response
% -------------------------------------------------------------------
plot_multiple(4,hrf,b,wd,wb,DX,mycolors,nrows,ncols,'No signal',1)


if dosel
	window = [0 round(32./TR)];
	% delta = DX(:,1:tp:size(DX,2)-1);

	subplot(nrows,ncols,5)
	trialavg2(id1,delta,window,'plot',1,'title','Selective Average, single response');
    legend off
	subplot(nrows,ncols,6)
	trialavg2(id2,delta,window,'plot',1,'title','Multiple response');
    legend off;title('')
	subplot(nrows,ncols,7)
	trialavg2(id3,delta,window,'plot',1,'title','Noisy multiple response');
    legend off;title('')
	subplot(nrows,ncols,8)
	trialavg2(id4,delta,window,'plot',1,'title','No signal');
    legend off;title('')
end
	






function plot_multiple(plotindex,hrf,b,wd,wb,DX,mycolors,nrows,ncols,mytitle,varargin)

subplot(nrows,ncols,plotindex)

plot(hrf)
myleg{1} = 'canonical HRF';
hold on; 

for i = 1:length(wd)
    h(i) = plot(b(wd(i):wd(i)+wb(i)),mycolors{i});
    myleg{i} = ['Condition' num2str(i)];
end

grid on
title(mytitle)
if length(varargin) > 0, legend(h,myleg), end

return
