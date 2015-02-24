function cleanup_masks
% clean up masks

% threshold masks and apply gray matter mask to all
% start in analysis_masks directory
in1 = input('Enter selection string 1: ','s');
in2 = input('Enter selection string 2: ','s');

D = dir;
ind = zeros(1,length(D));
ind2 = ind;

for i = 1:length(D),
    a = findstr(D(i).name,in1); 
    if ~isempty(a),ind(i) = a;,end
    
        a = findstr(D(i).name,in2); 
    if ~isempty(a),ind2(i) = a;,end
end

D = str2mat(D.name);
P = D(ind>0 & ind2 > 0,:);
V = spm_vol(P);

A1 = spm_get(1,'*img','Choose gray matter mask');
Av = spm_vol(A1);
Avol = spm_read_vols(Av);

vols = spm_read_vols(V);


  
% reslice anat if necessary
siz = size(vols);
if any(size(Avol) - siz(1:3))
    Avol = reslice_anat(P,A1);
end
        
% multiply masks by anatomical mask   
for i = 1:size(P,1)
    
    vols(:,:,:,i) = vols(:,:,:,i) .* Avol;
    
end

% threshold masks
vols = vols > .001;

% write output images   
for i = 1:size(P,1)
    
    spm_write_vol(V(i),vols(:,:,:,i));
    
end

return
    




function vol = reslice_anat(P,A1)
% returns volume of resliced anatomical image A1
    nms = reslice_imgs(P(1,:),A1,1);
                
    [d f e] = fileparts(A1);
               
    outnm = fullfile(d,['r' f e]); 
                
    V = spm_vol(outnm);
    vol = spm_read_vols(V);

return
