% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
% Routine to read time-series data.
%
% function [results, maskInfo] = SOM_PrepData(P,PMask)
%
%     P is array of file names (like that returned from spm_get)
%     PMask is a file name for a binary mask.
% 
% Only those voxels that are included in the mask are read.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [results, maskInfo] = SOM_PrepData(P,PMask)

fprintf('Reading Data\n\n');

maskInfo.header = spm_vol(PMask);

theMask = spm_read_vols(maskInfo.header);

maskInfo.iMask = find(theMask);

results = zeros(length(maskInfo.iMask),size(P,1));

for iP = 1:size(P,1)
    fprintf('\b\b\b%03d',iP);
    theVol = spm_read_vols(spm_vol(P(iP,:)));
    results(:,iP) = theVol(maskInfo.iMask);
    clear theVol;
end

fprintf('\nDone\n');
clear theVol;

return

%
% All done.
%