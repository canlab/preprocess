function [id1,b1,id2,b2] = ideal_deconv3(DX,tp,TR,dosel)
% function [id1,b1,id2,b2] = ideal_deconv3(DX,tp,TR,dosel)
% tests deconvolution matrix directly against idealized data
% you put in the exact temporal sequence to be deconvolved,
% in the form of the DX matrix.
%
% Tor Wager, 10/24/01
%
% inputs:
%	DX	deconvolution matrix
%	tp	time points estimated for each condition in DX
%	TR	repetition time of scan
%	dosel	1 or 0, test selective average method also

mycolors = {'ro-' 'b^--' 'go:' 'y^--' 'mo-' 'r^-' 'g^-' 'b^-' 'k^-' 'm^-' 'y^-'};

hrf = spm_hrf(TR);
hrf = hrf ./ max(hrf);


% plot dimensions - number of plots
ncols = 4;
if dosel, nrows = 2;, else, nrows = 1;, end


% -------------------------------------------------------------------
% * set up single response
% -------------------------------------------------------------------

ideal_data = conv(hrf,DX(:,1));
ideal_data = ideal_data(1:size(DX,1));

b = pinv(DX) * ideal_data;

id1 = ideal_data;
b1 = b;

% -------------------------------------------------------------------
% * plot single response
% -------------------------------------------------------------------
figure;
set(gcf,'Color','w')

subplot(nrows,ncols,1)
plot(hrf)
hold on; 
%plot(b(1:end-1),'rs-')
index = 1;
for i = 1:tp:size(DX,2)-1
    if index > length(mycolors),mycolors{index} = mycolors{index-7};,end
    plot(b(i:i+tp-1),mycolors{index})
    myleg{index} = ['Condition' num2str(index)];
    index = index + 1;
end
grid on
title('Single condition response')

% -------------------------------------------------------------------
% * set up multiple response
% -------------------------------------------------------------------

ideal_data = [];
for i = 1:tp:size(DX,2)-1
    mydata = conv(hrf,DX(:,i));
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
plot_multiple(2,hrf,b,tp,DX,mycolors,nrows,ncols,'Ideal response, all conditions')

% -------------------------------------------------------------------
% * set up noisy response
% -------------------------------------------------------------------

noise = .2 * rand(size(ideal_data));
ideal_data = ideal_data + noise;
b = pinv(DX) * ideal_data;

id3 = ideal_data;
b3 = b;

% -------------------------------------------------------------------
% * plot noisy response
% -------------------------------------------------------------------
plot_multiple(3,hrf,b,tp,DX,mycolors,nrows,ncols,'Noisy signal, all conditions')


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
plot_multiple(4,hrf,b,tp,DX,mycolors,nrows,ncols,'No signal')


if dosel
	window = [0 round(16./TR)];
	delta = DX(:,1:tp:size(DX,2)-1);

	subplot(nrows,ncols,5)
	trialavg2(id1,delta,window,'plot',1,'title','Selective Average, single response');

	subplot(nrows,ncols,6)
	trialavg2(id2,delta,window,'plot',1,'title','Multiple response');

	subplot(nrows,ncols,7)
	trialavg2(id3,delta,window,'plot',1,'title','Noisy multiple response');

	subplot(nrows,ncols,8)
	trialavg2(id4,delta,window,'plot',1,'title','No signal');
end
	






function plot_multiple(plotindex,hrf,b,tp,DX,mycolors,nrows,ncols,mytitle,varargin)

subplot(nrows,ncols,plotindex)

plot(hrf)
myleg{1} = 'canonical HRF';
hold on; 
index = 1;
for i = 1:tp:size(DX,2)-1
    plot(b(i:i+tp-1),mycolors{index})
    myleg{index+1} = ['Condition' num2str(index)];
    index = index + 1;
end

grid on
title(mytitle)
if length(varargin) > 0, legend(myleg), end

return
