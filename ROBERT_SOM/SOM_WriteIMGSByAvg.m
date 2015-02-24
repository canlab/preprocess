% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
% Ann Arbor MI.
%
% function results = SOM_Write(SOMResults)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function results = SOM_WriteIMGSByAvg(SOMResults)

curDIR = pwd;

[wts awt] = SOM_Weighted_SUM(SOMResults);

[oWTS oI] = sort(awt);

iIMG = 0;

for ix = 100:-1:91
    iIMG = iIMG + 1;
    hdr = SOMResults.header;
    hdr.fname = fullfile(curDIR,sprintf('som_byAvg_%03d.img',iIMG));
    hdr.pinfo = [1;0;0];
    vol = zeros(hdr.dim(1:3));
    vol(:,:,:) = nan;
    ii = find(SOMResults.IDX == oI(ix));
    vol(SOMResults.iMask(ii)) = SOMResults.WTS(ii);
    spm_write_vol(hdr,vol);
end

clear SOMResults
clear vol

return

%
% All done.
%
    
