function batch_montage3
% batch_montage3
%
% by Tor Wager
%
% start this in the directory above each random effects/spm contrasts directory
% containing spm rfx analyses; dirs w/sig results must contain *clusters.mat file.
%
% the script gets structures and clusters, and correlates each cluster
% with behavioral interference.
% This version works with correlation directories.

ovl = 'c:\tor_scripts\spm\templates\scalped_avg152T1.img';
dotext = input('Plot cluster numbers on figures? (1/0) ');

d = dir;

if(exist('EXPT.mat') == 2), load EXPT, else, error('No EXPT.mat file.  Create with get_expt_info.m'),end
prefix = input('Enter prefix of directories: ','s');
mylen = length(prefix);

for i = 3:length(d)

    myclust = [];
    if d(i).isdir & strcmp(d(i).name(1:mylen),prefix),
        
        eval(['cd ' d(i).name]),

        cli = 1;
        d2 = dir;
        for i = 3:length(d2), 
            if strcmp(d2(i).name(end-min(12,length(d2(i).name)-1):end),'_clusters.mat')

                load(d2(i).name)
                if exist('clusters') == 1
                    if ~isempty(clusters)
                        if length(clusters) < 100
                            myclust{cli} = clusters;
                            cli = cli + 1;
                        else
                            warning([d2(i).name ': More than 100 clusters - skipping.'])
                        end
                    end
                end
            end
        end
        

           
        if ~isempty(myclust)
            
            str = 'myclust{1}';
            for i = 2:length(myclust)
                str = [str ',myclust{' num2str(i) '}'];
            end
            str = [str ');']; 
            
            if dotext
                eval(['montage_clusters_text(ovl,' str])
            else
                eval(['montage_clusters(ovl,' str])
            end
            [dummy,myf,dummy] = fileparts(pwd);
            saveas(gcf,['..' filesep myf '_montage'],'jpg')
        else
            disp(['Warning: clusters variable does not exist in ' pwd])
        end
                
        
        cd ..
    end
    
    
end

return

        

