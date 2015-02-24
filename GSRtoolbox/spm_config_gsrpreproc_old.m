function c = spm_config_gsrpreproc)old(varargin)
% Configuration file for concatenation jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_3Dto4D.m 245 2005-09-27 14:16:41Z guillaume $

addpath(fullfile(spm('dir'),'toolbox','GSRtoolbox'));


data.type = 'files';
data.name = 'Raw GSR signal';
data.tag  = 'gsrsig';
data.num  = [1 Inf];
data.filter = '.*';
data.help = {'Select the signals to preprocess'};


fs.type = 'entry';
fs.name = 'Frequency of acquisition';
fs.tag  = 'fs';
fs.strtype = 'n';
fs.num  = [1 1];
fs.help = {'Enter the acquisition frequency'};


swin.type = 'entry';
swin.name = 'Smoothing window';
swin.tag  = 'swin';
swin.strtype = 'n';
swin.num  = [1 1];
swin.help = {'Enter window size for Gaussian smoothing on signal and derivitives'};


sigcol.type = 'entry';
sigcol.name = 'Signal column';
sigcol.tag  = 'sigcol';
sigcol.strtype = 'n';
sigcol.num  = [1 1];
sigcol.help = {'Enter that column number in the text file that corresponds to the GSR signal'};



trigcol.type = 'entry';
trigcol.name = 'Trigger column';
trigcol.tag  = 'trigcol';
trigcol.strtype = 'n';
trigcol.num  = [1 1];
trigcol.help = {'Enter that column number in the text file that corresponds to the event triggers'};



c.type = 'branch';
c.name = 'GSR preprocessing';
c.tag  = 'gsrpreproc';
c.val  = {data, fs, swin, sigcol, trigcol};
c.prog = @spm_gsrpreproc;
c.help = {'Normalise and smooth GSR signal and construct smoothed derivitives.',...
'Requires ascii text input.'};
end
%_______________________________________________________________________

%_______________________________________________________________________
function spm_gsrpreproc(job)
    
    swin = job.swin;
    fs = job.fs;
    spm_progress_bar('init');
    for i=1:length(job.gsrsig)
        fname = job.gsrsig{i};
        dat = load(fname);
        x = dat(:,job.sigcol);
        n = length(x);
        x = normalise(x);
        smoothx = smoothy(x,swin);
        dx = gradient(smoothx,1/fs);
        smoothdx = smoothy(dx,swin);
        ddx = gradient(smoothdx,1/fs);
        smoothddx = smoothy(ddx,swin);
        outdat = [x smoothx smoothdx smoothddx dat(:,job.trigcol)];
        
        dirdelim = find(fname=='\');
        fdir = fname(1:dirdelim(length(dirdelim)));
        fname = fname(dirdelim(length(dirdelim))+1:length(fname));
        dlmwrite([fdir 'prep' fname],outdat, ' ')
        
        spm_progress_bar('set', i/length(job.gsrsig));
    end 
    spm_progress_bar('clear');
end





function outsignal = smoothy(insignal, smoothwin)
    if smoothwin/2~=floor(smoothwin/2), smoothwin = smoothwin + 1; end;
    w = gausswin(smoothwin);
    outsignal = conv(insignal,w);
    outsignal = outsignal((smoothwin/2):(length(outsignal)-smoothwin/2))./sum(w);
end





function y = normalise(x)
    y = (x-mean(x))/std(x);
end