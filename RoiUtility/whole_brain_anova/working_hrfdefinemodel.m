% ----------------------------------------------------
% * Get the AVERAGE of all DX beta images for each condition
% ----------------------------------------------------

for i = 1:length(EXPT.subjects)
    
    d=dir('dx_beta*img');
    d = str2mat(d.name);
    v = spm_read_vols(spm_vol(d));

    % compute and save running sum of v
    % average v, and write avg dx beta images
    
    
end



% load the average dx beta images in v
    
    
    
olen = size(EXPT.DX.DX{1},1);

% ----------------------------------------------------
% * Loop through voxels in mask
% ----------------------------------------------------

[x,y,z] = ind2sub(size(EXPT.DX.mask),find(EXPT.DX.mask));

for i = 1:length(x)
    
    X{x,y,z} = get_X(EXPT,v,x(i),y(i),z(i),olen);
    

function get_X(EXPT,v,x(i),y(i),z(i),olen)

% ----------------------------------------------------
% * get HRF and build model
% ----------------------------------------------------

for i = 1:length(EXPT.DX.dxtrialonsets), 
    hrf{i} = squeeze(v(x,y,z,EXPT.DX.dxtrialonsets(i):EXPT.DX.dxtrialonsets(i)+EXPT.DX.numframes-1));
    
    m = conv(hrf{i},EXPT.DX.DX{1}(:,EXPT.DX.dxtrialonsets(i));
    X(:,i) = m(1:olen);
    X(:,i) = X(:,i) - mean(X(:,i));
    
end

X(:,end+1) = 1;

return