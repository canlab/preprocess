% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
% Routine to create analysis mask.
%
% function results = SOM_CreateMask(P)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_CreateMask(P)

fprintf('Create mask...\n');

dVol = spm_read_vols(spm_vol(P(1,:)));

dVol = ones(size(dVol));

for iP = 1:size(P,1)
    fprintf('\b\b\b%03d',iP)
    t(iP) = spm_global(spm_vol(P(iP,:)));
    vol   = spm_read_vols(spm_vol(P(iP,:)));
    dVol  = dVol.*(vol>t(iP));
    clear vol;
end

fprintf('\b\b\bdone\n');

theMask = dVol;
iMask = find(theMask);

hdr = spm_vol(P(1,:));

[pn fn] = fileparts(hdr.fname);

hdr.fname = fullfile(pn,'mask.img');

hdr.pinfo = [1;0;0];
hdr.dim(4) = 4;

spm_write_vol(hdr,theMask);

results = theMask;

clear theMask;

%
% All done.
%
