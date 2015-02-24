function spm_rfx(u,k,P)
%function spm_rfx(u,k,P)

%try,mkdir tmp_results,catch,end
%try,cd tmp_results
%catch,
%end

%BAT.covariate = EXPT.behavior;
BAT.multCompCorrect = 'uncorrected';
BAT.u = u;
BAT.k = k;
BAT.P = P;
BAT.npairs = size(BAT.P,1);
BAT.estnow = 1;
BAT.estcons = 1;    % estimate contrasts
BAT.dores = 1;      % get results

BAT.resdir = pwd
tor_spm_random_effects(BAT);