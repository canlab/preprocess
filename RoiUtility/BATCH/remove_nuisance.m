function EXPT = remove_nuisance(EXPT)
%
% This function removes nuisance covariates from contrast images before random
% effects analysis.
%
% It saves another set of contrast imgs, Rcon*.img, in individual subject results
% directories.  R is for "Residual"
%
% Nuisance covs are removed voxel by voxel, based on the set of contrast images
% from each subject of the same name.
%
% All covariates are centered in the function, and should be a matrix of column vectors.
%
% File names should be in EXPT.SNPM.P

covt = EXPT.cov;
if length(covt) > size(covt,1), covt = covt';, end
covt = covt - repmat(mean(covt),size(covt,1),1);

EXPT.cov = covt;
covt(:,end+1) = 1;

for i = 1:length(EXPT.SNPM.P)
    
    EXPT.SNPM.P{i} = rm_nuis(EXPT.SNPM.P{i},covt);
    
end



return



function newP = rm_nuis(P,covt)

V = spm_vol(P); v = spm_read_vols(V);
wh=sum(v,4);
wh = ~isnan(wh);
vo = zeros(size(v)) .* NaN;
fprintf(1,'\nworking on %s: %6.0f voxels, %3.0f planes\n\t',P(1,:),sum(sum(sum(wh))),size(v,1))

for i = 1:size(v,1)
    for j = 1:size(v,2)
        for k = 1:size(v,3)
            
            if wh(i,j,k)
                t = squeeze(v(i,j,k,:));
                b = pinv(covt) * t;
                r = t - (covt * b) + b(end);    % subtract fits but add intercept back in (preserve mean)
                vo(i,j,k,:) = r;
            end
            
        end
    end
    fprintf(1,'.')
end
fprintf(1,' done. ')
fprintf(1,'\twriting volumes.\n')

for i = 1:length(V)
    [d,f,e] = fileparts(V(i).fname);
    V(i).fname = fullfile(d,['R' f e]);
    V(i).descrip = 'Residuals from removed covariate';
    
    if i == 1, newP = V(i).fname;
    else
        newP = str2mat(newP,V(i).fname);
    end
    
    spm_write_vol(V(i),squeeze(vo(:,:,:,i)));
end


return