% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
% Modified by Tor Wager to create a meta-analysis appropriate mask.
% Use voxels in which AT LEAST 3 studies report some result, rather than ALL.
%
% Routine to create analysis mask.
%
% function results = SOM_CreateMask(P)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_CreateMask(P)

fprintf('Create mask...\n');

dVol = spm_read_vols(spm_vol(P(1,:)));

dVol = zeros(size(dVol));

for iP = 1:size(P,1)
    fprintf('\b\b\b%03d',iP)
    %t(iP) = spm_global(spm_vol(P(iP,:)));
    vol   = spm_read_vols(spm_vol(P(iP,:)));
    dVol  = dVol + (vol>0);
    clear vol;
end

dVol = double(dVol >= 3);   % thresh. for how many nonzero values there are.

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
