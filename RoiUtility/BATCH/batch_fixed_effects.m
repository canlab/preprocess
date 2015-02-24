function batch_fixed_effects(u,k)
% function batch_fixed_effects(u,k)
%
% Tor Wager
%
% run in directory above individual results directories

load EXPT

for i = EXPT.SNPM.connums
    
    clusters{i} = tor_fixed_effects2(i,u,k);
    
end



for i = 1:length(clusters)
    if ~isempty(clusters{i})
        conname = deblank(EXPT.SNPM.connames(EXPT.SNPM.connums == i,:));
        fprintf(1,'Contrast: %s',conname)
        if isfield(EXPT,'behavior')
            clusters{i} = simple_correl_clusters(clusters{i},EXPT.behavior,1);
        end
        cluster_table(clusters{i})
    end
end


if ~(exist('fixed_effect_results') == 2)
    mkdir fixed_effect_results
end

!move mean_spmT* fixed_effect_results\
!move fixed_results* fixed_effect_results\
cd fixed_effect_results
save fixed_clusters clusters

return

