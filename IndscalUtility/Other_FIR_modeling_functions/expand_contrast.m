function con = expand_contrast(contrast,numframes,DX)
% function con = expand_contrast(contrast,numframes,DX)
%
% expand contrast to cover all basis functions
% in an FIR response
% tor wager
%
% e.g., 
%c.contrast = [.5 .5 0 0 0 0; ...    % antic
%              0 0 .5 .5 0 0];       % pain
%c.numframes = [16 16 16 16 16 16];
% DX is 600 x 97 deconvolution FIR matrix (ONLY uses size)
% see tor_make_deconv_mtx3
% con = expand_contrast(c.contrast,c.numframes,DX);

con = [];       % all contrasts, expanded
for cind = 1:size(contrast,1) % for each contrast

    tmp = [];
    for i = 1:size(contrast,2),
        tmp = [tmp repmat(contrast(cind,i),1,numframes(i))]; % expand to cover all basis functions
    end
    con1 = tmp;
    z = zeros(1,size(DX,2) - length(con1));
    con1 = [con1 z];        % add zeros for intercepts, etc.
    con = [con; con1];
end

return
