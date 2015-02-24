function V = spm_write_4dvol(V,ofname,prefix)
% Write an image volume to disk, setting scales and offsets as appropriate
% FORMAT V = spm_write_4dvol(V,Y)
% V (input)  - a structure containing image volume information (see spm_vol)
%              including the cdf structure for time info
% ofname     - new name for file if desired. V contains the file to be read 
%              and split into timepoints
% V (output) - data structure after modification for writing.
%_______________________________________________________________________
% (c) SpeechLab, Boston University, 2002 Author: Satrajit Ghosh

if nargin<2,
    ofname = V.fname;
end;
if nargin<3,
    prefix = '';
end;

% Create a suffix template for filename
exttemplate = sprintf('_%s%%0%dd',prefix,5); %floor(log10(NTime))+2);

% Call default function when cdf struct is absent or contains only 1 time point
if ~isfield(V,'cdf') | ~isfield(V.cdf,'time_pts') | V.cdf.time_pts == 1,
    Y = spm_read_vols(V);
    if nargin>1,
        [pth,nm,xt] = fileparts(ofname);
        if isempty(pth),
            pth = pwd;
        end;
        V.fname = [pth,filesep,nm,sprintf(exttemplate,1),'.img'];
    end;
    V = spm_write_vol(V,Y);
    return;
end;


NTime = V.cdf.time_pts;

% read the data
Y = spm_read_vols(V);

% Modify header and reshape data
V.dim(3) = V.dim(3)/NTime;
V.pinfo = V.pinfo(:,1:NTime:end);

Y = reshape(Y,[V.dim(1:3),NTime]);


[pth,nm,xt] = fileparts(ofname);
if isempty(pth),
    pth = pwd;
end;

% Write out each individual time points
for ntime=1:NTime,
    V.fname = [pth,filesep,nm,sprintf(exttemplate,ntime),'.img'];
    spm_write_vol(V,squeeze(Y(:,:,:,ntime)));
end;