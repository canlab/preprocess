% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
% Routine to create analysis mask.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_CreateMask(P)

dVol = spm_read_vols(spm_vol(P(1,:)));

dVol = ones(size(dVol));

for iP = 1:size(P,1)
    fprintf('\b\b\b%03d',iP)
    vol = spm_read_vols(spm_vol(P(iP,:)));
    t(iP) = mean(reshape(vol,[1 prod(size(vol))]));
    dvol = dVol.*(vol>t(iP)/8);
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

%
% All done.
%
