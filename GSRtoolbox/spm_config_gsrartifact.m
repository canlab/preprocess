function c = spm_config_gsrartifact(varargin)
% Configuration file for concatenation jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_3Dto4D.m 245 2005-09-27 14:16:41Z guillaume $

addpath(fullfile(spm('dir'),'toolbox','GSRtoolbox'));


data.type = 'files';
data.name = 'GSR signal';
data.tag  = 'gsrsig';
data.num  = [1 Inf];
data.filter = '.*';
data.help = {'Select the signals to view'};


sigcol.type = 'entry';
sigcol.name = 'Signal column';
sigcol.tag  = 'sigcol';
sigcol.strtype = 'n';
sigcol.num  = [1 1];
sigcol.help = {'Enter that column number in the text file that corresponds to the GSR signal'};


dwin.type = 'entry';
dwin.name = 'Zoomed display size';
dwin.tag  = 'dispwin';
dwin.strtype = 'n';
dwin.num  = [1 1];
dwin.help = {'Enter that size of the artifact selection window.',...
    'This number should be greater than the width of the largest artifact, but small enough to allow for the selection of small artifacts.'};


c.type = 'branch';
c.name = 'GSR artifact definition';
c.tag  = 'gsrartifact';
c.val  = {data, sigcol, dwin};
c.prog = @spm_gsrartifact;
c.help = {'Tool for viewing a GSR signal and identifying artifacts.',...
'Requires ascii text input.'};
end
%_______________________________________________________________________

%_______________________________________________________________________
function spm_gsrartifact(job)
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)], 'UserData', num2str(job.dispwin));
    A = subplot(2,1,1);
    dat = load(job.gsrsig{1});
    P = plot(dat(:,job.sigcol));
    
    
    set(P, 'ButtonDownFcn', 'zoomIn;');
    set(A, 'ButtonDownFcn', 'zoomIn;');

    

    
end


