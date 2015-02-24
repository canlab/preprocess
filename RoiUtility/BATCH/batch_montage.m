function batch_montage
% batch_montage
%
% by Tor Wager
%
% start this in the directory above each random effects/spm contrasts directory
% containing spm rfx analyses; dirs w/sig results must contain *clusters.mat file.
%
% the script gets structures and clusters, and correlates each cluster
% with behavioral interference.
%
% see threshold_spm_t for info on fast thresholding and saving of clusters structs.

ovl = 'c:\tor_scripts\spm\templates\scalped_avg152T1.img';
ovl = which('scalped_single_subj_T1.img');
dotext = input('Plot cluster numbers on figures? (1/0) ');
prefix = input('Enter prefix of directories, e.g., rfx*: ','s');

d = dir(prefix);

if(exist('EXPT.mat') == 2), load EXPT, else, error('No EXPT.mat file.  Create with get_expt_info.m'),end

%mylen = length(prefix);

for i = 1:length(d)

    if d(i).isdir
        
        %eval(['cd ' d(i).name]),

        %d2 = dir;
        %for i = 3:length(d2), 
            %if strcmp(d2(i).name(end-min(12,length(d2(i).name)-1):end),'_clusters.mat')

                %load(d2(i).name)
        
        try
                if findstr(prefix,'rfx'),
                    load([d(i).name filesep 't_rfx0002_clusters'])
                elseif findstr(prefix,'corr')
                    load([d(i).name filesep 't_corr0002_clusters'])
                else
                    disp('Not sure which filename to load.')
                    load([d(i).name filesep 't_rfx0002_clusters'])
                end
                
                if exist('clusters') == 1
                    if ~isempty(clusters)
                        
                        if dotext
                            montage_clusters_text(ovl,clusters);
                        else
                            montage_clusters([],clusters,[2 2]);
                        end
                        %[dummy,myf,dummy] = fileparts(d2(i).name);
                        %saveas(gcf,['..' filesep myf '_montage'],'tif')
                        %[dummy,d2] = fileparts(pwd);
                        saveas(gcf,[filesep d(i).name '_montage_colbar'],'tif'),close
                        saveas(gcf,[filesep d(i).name '_montage'],'tif')
                        saveas(gcf,[filesep d(i).name '_montage'],'fig')
                    end
                else
                    disp(['Warning: clusters variable does not exist in ' pwd])
                    disp(['Name of clusters file is ' d2(i).name])
                end
                
            catch
                disp(['not found: ' d(i).name filesep 't_rfx0002_clusters'])
            end
                
            %end
        %end
        
       
        %cd ..
    end
    
    
end

return

        

