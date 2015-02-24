function nan_mask

maskP = spm_get(1,'*img','Select mask image');

P = spm_get(Inf,'*img','Select images to mask');

for i = 1:size(P,1)
    
    % f is input, Q is output file name
    f = deblank(P(i,:));

    [d,fn,e] = fileparts(f);
    Q = fullfile(d,['m' fn e]);
    
    f = str2mat(f,maskP);
        
    % apply brain mask
    spm_imcalc_ui(f,Q,'i1 .* i2',{0 1 16 0});

    % turn zeros into NaNs
    % though, instead, i just turned the brain mask into zeros!
    spm_imcalc_ui(Q,Q,'i1 + 0 ./ i1',{0 1 16 0});
    
end

return