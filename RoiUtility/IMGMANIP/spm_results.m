function [clusters,SPM,VY] = spm_results(u,varargin)
% function [clusters,SPM,VY] = spm_results(u,varargin)
%
% u is height threshold (p value)
% k is cluster size extent threshold

if length(varargin) > 0,k=varargin{1};,else,k=0;,end
TOR = struct('u',u,'k',k,'resdir',pwd)

TOR.multCompCorrect = 'uncorrected';
 TOR.maskWOtherContrasts = 0;
	% mask with other contrasts; 1 or 0

if isfield(TOR,'Ic')
    if ismepty(TOR.Ic)
        TOR.Ic = 2;
    end
else, TOR.Ic = 2;
end

[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(TOR);
spm_list('List',SPM,VOL,[],[],'',hReg);	
print -dpsc2 -painters -append -noui ../rfx_results.ps

% write results mask 
mask = voxel2mask(SPM.XYZ',VOL.DIM');
V = spm_vol('mask.img');   
V.fname = 'res_mask.img';
V = spm_write_vol(V,mask); 

% write filtered and clusters
V=spm_vol('spmT_0002.img');V2=spm_vol('res_mask.img');Vo=V;Vo.fname='spmT_filtered_0002.img';
Vo=spm_imcalc([V,V2],Vo,'i1 .* i2');

load SPM
p = str2mat(VY.fname);

% reslice into 2 x 2 x 2
%pt = which('scalped_single_subj_T1.img');
%[p1,p2]=reslice_imgs(pt,'spmT_filtered_0002.img');
%V=spm_vol(p2(2:end));vol=spm_read_vols(V); thr=abs(tinv(u,size(str2mat(VY.fname),1)-1));
%vol(vol < thr) = 0;
%spm_write_vol(V,vol);



clusters = tor_extract_rois(p,SPM,VOL);

for i = 1:length(clusters),clusters(i).Z = spm_t2z(clusters(i).Z,size(p,1)-1),end


spm_image('init',which('scalped_single_subj_T1.img'))
spm_orthviews('AddBlobs',1,SPM.XYZ,SPM.Z,VOL.M)

%h=findobj('Tag','Graphics'); figure(h); set(h,'Tag','old'); tmp=get(gcf,'Position');
%figure; set(gcf,'Tag','Graphics','Position',tmp); 
%spm_list('List',SPM,VOL,[],[],'',hReg);

%montage_clusters([],clusters,[2 2]); %set(gcf,'Position',tmp);

%figure('Color','w'); plot(tmp','bo','MarkerFaceColor','y');hold on; plot(mean(tmp),'rs','MarkerFaceColor','r','MarkerSize',8)
%xlabel('Cluster'),ylabel('Contrast values')
%clusters = cluster_table(clusters);

return
