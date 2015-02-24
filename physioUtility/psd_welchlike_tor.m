function yy1 = psd_welchlike_tor(x,binlen,binspace,samprate,minfreq,doplot)
% yy1 = psd_welchlike_tor(x,binlen,binspace,samprate,minfreq,doplot)
%
% Peak spectral power in a series of overlapping windows in a timeseries
% Uses Hanning window
% Resolution increases as length of bins (binlen)
%
% x = data
% samprate = samples/sec
% binlen = samples in a bin
% binspace = how far apart bin centers are, in elements
% minfreq = do not consider power below this frequency
%
% for heart rate, plausible Hz range is .8 to 2.5, 48 - 150
%
% Example: Heart rate sampled at 100 Hz
% 60 s Hanning window, get estimates every 1 s (100 samples)
% yy1 = psd_welchlike_tor(physio_data.pulse,100*60,100,100,50/60,1);

% [SMap, Sig, CP, OOCcnt, TWidMap, W1Map, W2Map, W3Map, Tsq,stats] =
% ewma4(X, .1 , .05, 'WN',120,0);

h = hanning(binlen+1);
len = length(x);

% Not quite right for this for some reason
%TR = 1 / samprate;
%nyq = 1 / (2 * TR);
%secs = (1:len) ./ samprate;
%Hz = (1:len) / (len * TR);
%Hz2 = linspace(1,samprate,length(x));

[Pxx,W] = pwelch(x,[],[],[],samprate);  % get W, frequencies

x = x - nanmean(x);

pad1 = x(1:binlen);
pad2 = x(end-binlen:end);
x = [pad1;x;pad2];

xx = binlen+1:binspace:length(x)-binlen;        % elements at which to sample

if doplot, f1 = tor_fig; hold off; f2 = tor_fig; end;
    
ind = 1;
for i = xx

    sig = x(i-binlen./2:i+binlen./2);   % centered around i
    sig = sig - nanmean(sig);
    
    %y = psd(sig .* h);
    [y,WW] = pwelch(sig .* h,[],[],[],samprate);    % WW = frequencies
    
    y(WW <= minfreq) = 0;                         % do not consider freqs below or at minfreq
    
    yy(ind) = WW(find(y == max(y)));              % frequency at max power
    
    
    % output
    if ind == 1, fprintf(1,'Resolution: %3.2f Hz, %3.2f cycles/min\n',WW(2)-WW(1),(WW(2)-WW(1))*60), end
    
    if doplot, figure(f1); plot(sig); title(['Signal, Epoch ' num2str(ind) ' of ' num2str(length(xx))]); 
        figure(f2);   pwelch(sig .* h,[],[],[],samprate); title(['Freq: ' num2str(yy(ind))]);
        drawnow
        doplot = input('Press 1 to keep plotting, or 0 to skip:');
    end
    
    ind = ind+1;
    
end

% linear interpolation
yy1 = interp1(xx',yy,1:len,'linear','extrap');




return