function [id1,b1,id2,b2] = ideal_deconv4(DX,TR,dosel,varargin)
% function [id1,b1,id2,b2] = ideal_deconv4(DX,TR,dosel,[opt] contrast, [opt] shiftby)
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

if length(varargin) > 1
    shiftby = varargin{2};
else
    shiftby = .5;
end

mycolors = {'ro-' 'b^--' 'go:' 'y^--' 'mo-' 'r^-' 'g^-' 'b^-' 'k^-' 'm^-' 'y^-'};

[delta, wd, wb] = DX_find_delta(DX);
delta = delta(:,wb > 0);
disp(['Conditions detected: ' num2str(size(delta,2))])

hrf = spm_hrf(.1);
hrf = hrf ./ max(hrf);


% plot dimensions - number of plots
ncols = 1;
if dosel, nrows = 2;, else, nrows = 1;, end




for mysh = 1:length(shiftby)
    
% -------------------------------------------------------------------
% * set up response to conditions specified in cond2c
% -------------------------------------------------------------------

%ideal_data = conv(hrf,delta(:,1));
%ideal_data = ideal_data(1:size(DX,1));

ideal_data = [];
myd = sum(delta * cond2c',2);
myd = sampleinseconds(myd,TR);  % upsample to .1 s

myd = conv(hrf,myd);

% shift here - by .5 TR
if shiftby(mysh) > 0
    myd = [zeros(TR*shiftby(mysh)/.1,1); myd];
elseif shiftby(mysh) < 0
    myd(1:(TR*abs(shiftby(mysh)./.1))) = [];
end

% resample at TR
mydata = resample(myd,1,TR*10);
ideal_data = mydata(1:size(DX,1));


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
plot(resample(hrf,1,TR*10),'LineWidth',2)
hold on; 
%plot(b(1:end-1),'rs-')

for i = 1:length(wd)
    if i > length(mycolors),mycolors{i} = mycolors{i-7};,end
    myh(i) = plot(b(wd(i):wd(i)+wb(i)),mycolors{i});
    myleg{i} = ['Condition' num2str(i)];
end
grid on
title(['Contrast ' num2str(cond2c) ' SNR = 1 shiftby: ' num2str(shiftby(mysh)) ' TRs'])
legend(myh,{'1' '2' '3' '4' '5' '6' '7' '8'})


if dosel
	window = [0 round(32./TR)];
	% delta = DX(:,1:tp:size(DX,2)-1);

	subplot(nrows,ncols,2)
	trialavg2(id1,delta,window,'plot',1,'title','Selective Average, single response');
	
end
	

end % shiftby

return



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
