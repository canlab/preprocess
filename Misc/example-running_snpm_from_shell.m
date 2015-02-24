[torw@linux-jjonides01 snpm0002]$ matlab << EOF >& snpm_run.txt
cwd = pwd; cd ~; setup; cd(cwd) 
load SnPMcfg.mat; vFWHM = [10 10 10]; bVarSm = 1;
save SnPMcfg
spm_snpm('.')
exit
EOF

cd /data/placebo/RESULTS/model5_test/whole_brain_n23/var_smoothed/snpm0002

cd /data/placebo/RESULTS/model5_test/whole_brain_n23/var_smoothed/snpm0003

matlab << EOF >& snpm_run.txt
cwd = pwd; cd ~; setup; cd(cwd) 
load SnPMcfg.mat; vFWHM = [10 10 10]; bVarSm = 1;
save SnPMcfg
spm_snpm('.')
exit
EOF

tail -f snpm_run.txt

% robust 2
cd /data/sim
try,mkdir robust2,catch,end
cd robust2
matlab << EOF >& robust2.out
cwd = pwd; cd ~; setup; cd(cwd) 
disp('starting')
robust_vs_ols_corr2(100,[10 20 30 50 100]);
disp('done')
exit
EOF

