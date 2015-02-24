function EXPT = get_nlcon_names_for_rfx(EXPT)
%
% tor wager
% start in directory above indivudal results
% extracts the following images:
%   dx_beta*img     hrf deconvolution estimates for each voxel
%   cond*img        nonlinear fits to height, delay, and intercept for each
%                   trial type (cond indexes trial type)
%   nlcon*img       contrast across parameters estimated by nonlinear fits
%
% see group_fit_model, group_brain_nlfit, nl_contrasts

for i = 1:length(EXPT.subjects)
    
    cd(EXPT.subjects{i})
    
    d = dir('dx_beta*.img');
    for j = 1:length(d)
        
        if j == 1
            EXPT.dx_beta_imgs{i} = [pwd filesep d(j).name];
        else
            EXPT.dx_beta_imgs{i} = str2mat(EXPT.dx_beta_imgs{i},[pwd filesep d(j).name]);
        end
        
    end
    
    for j = 1:length(d)
        
        if i == 1
            EXPT.NLCON.P{j} = [pwd filesep d(j).name];
        else
            EXPT.NLCON.P{j} = str2mat(EXPT.NLCON.P{j},[pwd filesep d(j).name]);
        end
        
    end
    
    d = dir('nlcon*.img');
    
    for j = 1:length(d)
        
        if i == 1
            EXPT.NLCON.P{j} = [pwd filesep d(j).name];
        else
            EXPT.NLCON.P{j} = str2mat(EXPT.NLCON.P{j},[pwd filesep d(j).name]);
        end
        
    end
    
    d = dir('cond*.img');
    
    for j = 1:length(d)
        
        if i == 1
            EXPT.NLCON.condP{j} = [pwd filesep d(j).name];
        else
            EXPT.NLCON.condP{j} = str2mat(EXPT.NLCON.condP{j},[pwd filesep d(j).name]);
        end
        
    end
    
    cd ..
    
end

return