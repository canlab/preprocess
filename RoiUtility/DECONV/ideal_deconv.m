function [ideal_hrf, ideal_data, DX, null_hrf, rand_data] = ideal_deconv(DX,sf,tp,TR,eres,doplot)
% [ideal_hrf, ideal_data,DX] = ideal_deconv(DX,sf,tp,TR,eres,doplot)
%
% DX is deconv matrix; if empty input, makes DX. 
% sf is hi-res stick function cell array, 1 cell per condition
% tp is time points in ideal HRF and estimate
% TR is repetition time of scans 
% eres is number of elements per TR in sf
% plot is 1 or 0, plot or don't
%
% Tor Wager, 10/24/01

if isempty(DX)    
    % eres = input(['Length of sf is ' size(sf{1},1) '.  Enter samples per second in sf: ']);
    [DX,short_sf] = tor_make_deconv_mtx2(sf,tp,eres);
end

hrf = spm_hrf(TR ./ eres);
hrf = hrf ./ max(hrf);


% -------------------------------------------------------------------
% * make ideal data
% -------------------------------------------------------------------

ideal_data = make_ideal_data(hrf,sf,eres,size(DX,1));

% -------------------------------------------------------------------
% * fit deconv
% -------------------------------------------------------------------
try    
    b = pinv(DX) * ideal_data;
catch
    whos ideal_data
    whos DX
    b = pinv(DX) * ideal_data;
end

% -------------------------------------------------------------------
% * break betas into hrf estimates
% -------------------------------------------------------------------

index = 1;
for i = 1:tp:length(b)-1
	ideal_hrf{index} = b(i:i+tp-1);
	index = index + 1;
end

% -------------------------------------------------------------------
% * fit deconv to RANDOM data
% -------------------------------------------------------------------
rand_data = rand(size(DX,1),1);
b2 = pinv(DX) * rand_data;
index = 1;
for i = 1:tp:length(b2)-1
	null_hrf{index} = b2(i:i+tp-1);
	index = index + 1;
end


% -------------------------------------------------------------------
% * plot
% -------------------------------------------------------------------
if doplot

    mycolors = {'ro-' 'b^--' 'go:' 'y^--' 'mo-' 'r^-' 'g^-' 'b^-' 'k^-' 'm^-' 'y^-'};

    figure
    subplot(1,2,1)
    plot(hrf,'LineWidth',2)
    hold on
    resamp_hrf = resample(hrf,1,eres);
    resamp_x = resample(1:length(hrf),1,eres);
    plot(resamp_x(1:tp),resamp_hrf(1:tp),'ks','MarkerFaceColor','k')
    grid on
    title('Ideal response')
    myleg{1} = 'ideal hrf';
    myleg{2} = 'best estimates';
    
    hold on 
    grid on
    title('Estimated hemodynamic responses')
    for i = 1:length(ideal_hrf)
    	try
            plot(resamp_x(1:tp),ideal_hrf{i},mycolors{i},'LineWidth',2),hold on
    	    myleg{i+2} = ['Condition' num2str(i)];
        catch
            warning('Can''t plot: trying to plot beyond index values?')
            tp
            whos resamp_x
            ideal_hrf
            mycolors
        end
    end
    legend(myleg,0)
    
    subplot(1,2,2)
    hold on 
    grid on
    title('Estimated responses to null data')
    plot(hrf,'LineWidth',2)
    plot(resamp_x(1:tp),resamp_hrf(1:tp),'ks','MarkerFaceColor','k')
    for i = 1:length(null_hrf)
    	plot(resamp_x(1:tp),null_hrf{i},mycolors{i},'LineWidth',2),hold on
    end
end

return