%-----------------------------------------------------------------------
% Job saved on 24-May-2017 12:58:08 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%

				% Adding SPM12 to the path
addpath('/work/ics/data/projects/wagerlab/Resources/spm12');

		  %Defining inputs for a general matlabbatch function.


spm('defaults', 'FMRI');
spm_jobman('initcfg');
basedirname= '/work/ics/data/projects/wagerlab/labdata/collab/Mackey_Duloxetine/';
load preproc_filenames.mat;

for i=1:length(subnames)
  clear('matlabbatch','volnames','warpvols','nvols');
  subname= subnames{i};
  warp_fn= strcat(basedirname,towarp_fns{i});
  structural_fn= strcat(basedirname,structural_fns{i},',1');
  functional_fn= strcat(basedirname, functional_fns{i});

				%Deriving matlabbatch inputs
  nvols=length(spm_vol(functional_fn));
  for j=1:nvols
    volnames{j}=strcat(functional_fn, ',',num2str(j));
    warpvols{j}=strcat(warp_fn,',',num2str(j));
  end
  volnames=volnames';
  warpvols=warpvols';
  destdirname=strcat(basedirname,'/Imaging/',subname);
				%Actual matlabbatch script
  matlabbatch{1}.spm.temporal.st.scans = {volnames};
  %%
  matlabbatch{1}.spm.temporal.st.nslices = 30;
  matlabbatch{1}.spm.temporal.st.tr = 2;
  matlabbatch{1}.spm.temporal.st.ta = 1.93333333333333;
  matlabbatch{1}.spm.temporal.st.so = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];
  matlabbatch{1}.spm.temporal.st.refslice = 1;
  matlabbatch{1}.spm.temporal.st.prefix = 'a';
  matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
  matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
  matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
  matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
  matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
  matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
  matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
  matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
  matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
  matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
  matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
  matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
  matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
  matlabbatch{3}.spm.spatial.preproc.channel.vols = {structural_fn};
  matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001;
  matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60;
  matlabbatch{3}.spm.spatial.preproc.channel.write = [0 1];
  matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm = {'/work/ics/data/projects/wagerlab/Resources/spm12/tpm/TPM.nii,1'};
  matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1;
  matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm = {'/work/ics/data/projects/wagerlab/Resources/spm12/tpm/TPM.nii,2'};
  matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1;
  matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = {'/work/ics/data/projects/wagerlab/Resources/spm12/tpm/TPM.nii,3'};
  matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
  matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = {'/work/ics/data/projects/wagerlab/Resources/spm12/tpm/TPM.nii,4'};
  matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3;
  matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [1 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = {'/work/ics/data/projects/wagerlab/Resources/spm12/tpm/TPM.nii,5'};
  matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
  matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [1 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = {'/work/ics/data/projects/wagerlab/Resources/spm12/tpm/TPM.nii,6'};
  matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2;
  matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
  matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
  matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
  matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1;
  matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
  matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
  matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
  matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
  matlabbatch{3}.spm.spatial.preproc.warp.write = [1 1];
  matlabbatch{4}.spm.util.imcalc.input(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
  matlabbatch{4}.spm.util.imcalc.input(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
  matlabbatch{4}.spm.util.imcalc.input(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
  matlabbatch{4}.spm.util.imcalc.output = 'in_brain_mask';
  matlabbatch{4}.spm.util.imcalc.outdir = {strcat(basedirname,'/Imaging/',subname,'/Structural/SPGR')};
  matlabbatch{4}.spm.util.imcalc.expression = 'double(i1|i2|i3)';
  matlabbatch{4}.spm.util.imcalc.var = struct('name', {}, 'value', {});
  matlabbatch{4}.spm.util.imcalc.options.dmtx = 0;
  matlabbatch{4}.spm.util.imcalc.options.mask = 0;
  matlabbatch{4}.spm.util.imcalc.options.interp = 0;
  matlabbatch{4}.spm.util.imcalc.options.dtype = 4;
  matlabbatch{5}.spm.util.imcalc.input(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
  matlabbatch{5}.spm.util.imcalc.input(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
  matlabbatch{5}.spm.util.imcalc.input(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
  matlabbatch{5}.spm.util.imcalc.input(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
  matlabbatch{5}.spm.util.imcalc.output = 'masked_anatomical';
  matlabbatch{5}.spm.util.imcalc.outdir = {strcat(basedirname,'/Imaging/',subname,'/Structural/SPGR')};
  matlabbatch{5}.spm.util.imcalc.expression = 'i1.*double(i2|i3|i4)';
  matlabbatch{5}.spm.util.imcalc.var = struct('name', {}, 'value', {});
  matlabbatch{5}.spm.util.imcalc.options.dmtx = 0;
  matlabbatch{5}.spm.util.imcalc.options.mask = 0;
  matlabbatch{5}.spm.util.imcalc.options.interp = 1;
  matlabbatch{5}.spm.util.imcalc.options.dtype = 4;
  matlabbatch{6}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Image Calculator: ImCalc Computed Image: masked_anatomical', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
  matlabbatch{6}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
  matlabbatch{6}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
  matlabbatch{6}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Image Calculator: ImCalc Computed Image: in_brain_mask', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
  matlabbatch{6}.spm.spatial.coreg.estimate.other(3) = cfg_dep('Image Calculator: ImCalc Computed Image: masked_anatomical', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
  matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
  matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
  matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
  matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
  matlabbatch{7}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
  matlabbatch{7}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
  matlabbatch{7}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
							    78 76 85];
  matlabbatch{7}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
  matlabbatch{7}.spm.spatial.normalise.write.woptions.interp = 4;
  matlabbatch{7}.spm.spatial.normalise.write.woptions.prefix = 'w';
  matlabbatch{8}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
  %%
  matlabbatch{8}.spm.spatial.normalise.write.subj.resample = warpvols;
  %%
  matlabbatch{8}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
							    78 76 85];
  matlabbatch{8}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
  matlabbatch{8}.spm.spatial.normalise.write.woptions.interp = 4;
  matlabbatch{8}.spm.spatial.normalise.write.woptions.prefix = 'w';
  matlabbatch{9}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
  matlabbatch{9}.spm.spatial.smooth.fwhm = [6 6 6];
  matlabbatch{9}.spm.spatial.smooth.dtype = 0;
  matlabbatch{9}.spm.spatial.smooth.im = 0;
  matlabbatch{9}.spm.spatial.smooth.prefix = 's';
  %Save the matlabbatch file for this subject.
  save(strcat('preproc_batches/duloxetine_preproc_batch_',num2str(i),'.mat'),'matlabbatch','destdirname');
end
