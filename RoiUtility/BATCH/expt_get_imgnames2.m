function EXPT = expt_get_imgnames2(EXPT,str)
% function expt_get_imgnames2(EXPT,str)
% alternate way to load contrast or other images into EXPT.SNPM.P for data extraction 
% or random effects model fitting
%
% str is a wildcard string for the images - e.g., 'Rcon*img'
% all images found will be put in cells of EXPT.SNPM.P, in row order of EXPT.subjects

cwd = pwd;

EXPT.SNPM.P = [];
for subj = 1:length(EXPT.subjects)
    
    cd(EXPT.subjects{subj})
    
    dd = dir(str);
    
    if subj == 1
        for con = 1:length(dd), EXPT.SNPM.P{con} = which(dd(con).name);, end
    else
        for con = 1:length(dd), EXPT.SNPM.P{con} = str2mat(EXPT.SNPM.P{con},which(dd(con).name));, end
    end
    
    cd ..
    
end

cd(cwd)

return