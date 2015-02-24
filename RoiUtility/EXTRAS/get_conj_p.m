% start in rfx directory

load SPM
u = spm_uc_rf(.05,[1 11],'T',R,12)
p = 1- spm_tcdf(.5715,[11])