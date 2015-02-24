function batch_montage2
% batch_montage2
%
% by Tor Wager
%
% start this in the directory above each random effects/spm contrasts directory
% containing spm rfx analyses
%
% this script works on pairs of contrasts assumed to be activations and deactivations
% respectively
%
%

ovl = 'c:\tor_scripts\spm\templates\scalped_avg152T1.img';

d = dir;

if(exist('EXPT.mat') == 2), load EXPT, else, error('No EXPT.mat file.  Create with get_expt_info.m'),end
gonow = 0;
index = 1;

for i = 3:length(d)

    if d(i).isdir & strcmp(d(i).name(1:3),'rfx'),eval(['cd ' d(i).name]),end
    
    if(exist('SPM.mat') == 2)
        
        a = pwd;
        a = a(end-6:end);
        a = [a '_clusters'];
        load(a)

        if mod(index,2) == 0
            montage_clusters(ovl,clusters1,clusters);
        else
            clusters1 = clusters;
            name1 = d(i).name;
        end
        index = index + 1;
    end
    
    if d(i).isdir & strcmp(d(i).name(1:3),'rfx'), 
        saveas(gcf,[name1 '_' d(i).name '_montage'],'jpg')
        cd ..
    end
    
    
end

return

        

