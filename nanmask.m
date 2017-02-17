function out = nanmask(P,varargin)
% out = nanmask(image_to_mask,masking_image)
%
% tor wager
% turns zeros into NaN values in the first input image
% input = a string with the image file name
%
% a second input will apply a mask to the first image,
% including only voxels in the output image where the mask image values are
% greater than zero.
% an empty 2nd argument calls the spm gui to ask for a mask image
if length(varargin) > 0, 
    maskP = varargin{1}; ,
    if isempty(maskP)
        maskP = spm_get(1,'*img','Select mask image');
    end
end

if isempty(P)
    P = spm_get(Inf,'*img','Select images to mask');
end

for i = 1:size(P,1)
    
    % f is input, Q is output file name
    f = deblank(P(i,:));

    [d,fn,e] = fileparts(f);
    Q = fullfile(d,['m' fn e]);
    
    if exist('maskP') == 1
        
        f = str2mat(f,maskP);
        
        % apply brain mask
        spm_imcalc_ui(f,Q,'i1 .* (i2 > 0)',{0 1 spm_type('float') 0});

        % turn zeros into NaNs
        spm_imcalc_ui(Q,Q,'i1 + 0 ./ i1',{0 1 spm_type('float') 0});
    
    else
 
        % turn zeros into NaNs
        spm_imcalc_ui(f,Q,'i1 + 0 ./ i1',{0 1 spm_type('float') 0});
        
    end
        
    if i == 1, 
        out = Q;
    else
        out = str2mat(out,Q);
    end
    
end
return

