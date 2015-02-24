% flip_if_neg_dim(imgs, [verbose])
% 
% Flips an Analyze image's X dimension if it is negative
%
% Only influences the .mat file associated with an image, not the image
% data itself.
%
% tor modified 8/2012 to handle 4-D files and missing M variable
% but does not always seem to work right. use with caution.
% this is mainly a legacy script from including FSL preprocessing.

function flip_if_neg_dim(imgs, varargin)
    flip_X_pixdim(imgs, 'flip_negative_only', 1, varargin{:});
end