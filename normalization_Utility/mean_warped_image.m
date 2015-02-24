function [meandata,fitness,prm] = mean_warped_image(imagenames,prmnames,varargin)
% [meandata,fitness,prm] = mean_warped_image(imagenames,prmnames,varargin)
%
% Create weighted mean of best 50% of normalized images
% 1) Apply normalization parameters in prmnames sn3d .mat files to images
%   in imagenames
% 2) Write individual warped images (optional.  Enter 'write' keyword)
% 3) Create weighted mean image and write
%
% imagenames : string matrix of images to apply warps to
% prmnames : str matrix of sn3d .mat file names
% 
% Optional inputs:
% 'write'
%  volTemplate structure, volume info for template/mask.  See spm_vol
%       Defaults to this template image: prm{1}.VG.fname
%       Entering a mask image of the same dims will give masked fitness
%       weights
%
% Write only trimmed mean image:
% [meandata,fitness] = mean_warped_image(imagenames,prmnames);
%
% Write individual images:
% [meandata,fitness] = mean_warped_image(imagenames,prmnames,'write');
%
% Enter a tempate/mask volume to define dimensions
% 
% tor wager,  nov. 06, version 12/2


fprintf(1,'\n-----------------------------\n');
fprintf(1,'mean_warped_image.m\n');    
fprintf(1,'1) Applying norm params\n');
fprintf(1,'2) Weighted mean of best 50%%\n\n');


t1 = clock;

% Load volume info, warping parameters
% ------------------------------------
n = size(imagenames,1);
clear dat
clear volInfo
clear prm

V = spm_vol(imagenames);
for i = 1:n
    
    % read norm params
    prm{i} = load(deblank(prmnames(i,:)));
    
end

% Get mask and template info
% ------------------------------------
% this should be in the same space as the template
% % % mask = which('scalped_avg152T1_graymatter_smoothed.img');
% % % [volTemplate] = iimg_read_img(mask,2);

for i = 1:length(varargin)
    if strcmp(varargin{i},'template') || strcmp(varargin{i},'mask')
        template = varargin{i+1}; 
    end
end
    
if ~exist('template','var')
    [dd,ff,ee] = fileparts(prm{1}.VG.fname);
    template = which([ff ee]);
end

[volTemplate,maskdat] = iimg_read_img(template,2);
% % not needed because volTemplate contains list  : maskdat = maskdat ~= 0 && ~isnan(maskdat);  % logical


% Apply warps
% ------------------------------------
clear normdat
fprintf(1,'Applying warps: ');

for i = 1:n
    fprintf(1,'%03d',i);
    
    % apply warps
    voldata = spm_read_vols(V(i));
    
    ndat = norm_apply_params(prm{i},'data',voldata,'mat',V(i).mat);
    
    if any(strcmp(varargin,'write'))
        Vout = volTemplate;
        [dd,ff,ee] = fileparts(V(i).fname);
        Vout.fname = sprintf('w%03d_%s%s',i, ff, ee);
        fprintf(1,' Writing: %s\n   ',Vout.fname);
        spm_write_vol(Vout,ndat);
    end
        

    normdat(:,i) = ndat(volTemplate.wh_inmask);
    
    fprintf(1,'\b\b\b');
end
fprintf(1,'Done.\n');

fprintf(1,'Creating mean: ');

% Provisional mean
% ------------------------------------
fprintf(1,'Initial. ');

normdat(isnan(normdat)) = 0;

% convert in-mask voxels to mean = gm, so averaging is not distorted by
% overall intensity differences among images
% and fitness weights determine relative contribution of images
means = mean(normdat);
gm = mean(means);           % global mean
scalef = gm ./ means;       % scaling factor to norm each image to grand mean

for i = 1:n
    normdat(:,i) = normdat(:,i) .* scalef(i);
end

meandat = mean(normdat,2);
vox = size(meandat,1);


fprintf(1,'Weighted mean of best 50%%. ');
% Weights based on closeness to mean
% ------------------------------------
for i = 1:n
    fitness(i) = norm_eval_fitness(meandat,normdat(:,i),0);
end

mfit = nanmedian(fitness);
weights = fitness;
weights(weights < mfit) = 0;

weights = weights ./ sum(weights);

if all(weights == 0)
    warning('All weights are zero!')
    weights = ones(1,n);
end

% Weighted mean
% ------------------------------------
meandat = zeros(vox,1);
for i = 1:n
    
    if weights(i) > 0   % just to save time
        meandat = meandat + weights(i) .* normdat(:,i);
    end
    
end

fprintf(1,'Done: Writing\n');
% Write output image
% This will be the new normalization template
% ------------------------------------
meandata = iimg_reconstruct_vols(meandat,volTemplate,'outname','mean.img');

fprintf(1,'Mean warped image completed: %3.0f s\n',etime(clock,t1));

fprintf(1,'-----------------------------\n');

return