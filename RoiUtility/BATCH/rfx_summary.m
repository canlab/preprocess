function rfx_summary
% function rfx_summary
%
% Start in dir above SPM random effects or SnPM analyses
% reports critical thresholds for bonferroni, random fields,
% FDR and/or SnPM t or pseudo-t height/cluster size
%
% Prints table of output for each subdirectory


mypwd = pwd;

look4 = input('Enter wildcard for analysis dirs (e.g., rfx*): ','s');
D = dir(look4);

diary([mypwd filesep 'rfx_summary.txt'])

if strcmp(look4,'rfx*') | strcmp(look4,'corr*')
            fprintf(1,'_______________________________________________________________________________________\n')
            fprintf(1,'directory\tdf\tvoxels\tPcorr\tbonfu\tGRFu\tFDRu\tHochu\tSidaku\tmax_obsT\tbonf_n\tGRF_n\tFDR_n\tHoch_n\tSidak_n\tSmoothness\t\t\n')
            fprintf(1,'_______________________________________________________________________________________\n')
		look4 = 1;
    
elseif strcmp(look4,'snpm*')
            fprintf(1,'_______________________________________________________________________________________\n')
            fprintf(1,'directory\tdf\tvoxels\tPmin\tu\tk\tvar_smooth\tmax_obsT\tcrit_n\tST_voxels\t\n')
            fprintf(1,'_______________________________________________________________________________________\n')
		look4 = 2;

else error('Enter rfx*, corr*, or snpm*')
end


for i = 1:length(D)
    
    if D(i).isdir, 
        
        eval(['cd(''' D(i).name ''')'])

        
        
        
        % -----------------------------------------------------------
        % * SPM
        % -----------------------------------------------------------
        
        if exist('SPM.mat') == 2 & look4 == 1
            load SPM
            if exist('spmT_0002.img') == 2
                Vs = spm_vol('spmT_0002.img');

            
            diary off
            TOR.maskWOtherContrasts = 0;
            TOR.multCompCorrect = 'FWE'; 
            TOR.u = .05;
            TOR.k = 0;
            TOR.resdir = pwd;
            TOR.Ic = 2;
            [SPM,VOL] = tor_spm_getSPM(TOR);
            diary on
            
            df = SPM.df(2); % length(VY) - 1;
            vFWHM = [];

            u = spm_uc_Bonf(.05,[1 df],'T',S,1);
            u2 = spm_uc_RF(.05,[1 df],'T',R,1);
            u3 = spm_uc_FDR(.05,[1 df],'T',1,Vs,0);
	        u4 = spm_uc_Hoch(.05,[1 df],'T',1,Vs,0);
            u5 = spm_uc_Sidak(.05,[1 df],'T',1,Vs,0);
         
            % power stuff
            vol = spm_read_vols(Vs);
            mT = max(max(max(vol)));
            n1 = power_calc(mT,u,df+1);
            n2 = power_calc(mT,u2,df+1);
            n3 = power_calc(mT,u3,df+1);
	        n4 = power_calc(mT,u4,df+1);
            n5 = power_calc(mT,u5,df+1);
            
	    % counting suprathreshold voxels          
	    vol = vol(:); vol(vol == 0 | isnan(vol)) = [];
	    count1 = sum(vol >= u);
	    count2 = sum(vol >= u2);
  	    count3 = sum(vol >= u3);
	    count4 = sum(vol >= u4);
        count5 = sum(vol >= u5);
        
            p = tdist(mT,round(df));
            
            k = []; % extent thresh
            S = VOL.S; % voxels
            
            [dummy mywd] = fileparts(pwd);
            
            %fprintf(1,'%s\n',pwd)

            % directory\tdf\tvoxels\tPcorr\tbonfu\tGRFu\tFDRu\tHochu\tSidaku\tmax_obsT\tbonf_n\tGRF_n\tFDR_n\tHoch_n\tSidak_n\tSmoothness\t\t\n
            fprintf(1,'%s\t%3.0f\t%6.0f\t%2.6f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%3.2f\t%3.2f\t%3.2f\t\n',mywd,df,S,p,u,u2,u3,u4,u5,mT,count1,count2,count3,count4,count5,VOL.FWHM(1),VOL.FWHM(2),VOL.FWHM(3))
            %fprintf(1,'\n')
            
  
            else
                diary off
                disp([pwd ' No spmT_0002.img? skipping...'])
            end
        end
        
        % -----------------------------------------------------------
        % * SnPM
        % -----------------------------------------------------------
        
        if exist('SnPM.mat') == 2 & look4 == 2
            load SnPM
            u = prctile(MaxT(:),95);     
            
            
            if exist('SnPM_pp.mat') == 2
                load SnPM_pp
                k = prctile(MaxSTCS,95);
            else
                k = 0;
            end
            
            load SnPMcfg
            if ~exist('vFWHM') == 1, vFWHM = 0;, end
            if isempty('vFWHM'),vFWHM = 0;, end
                
            % power stuff
            Vs = spm_vol('SPMt.img');
            vol = spm_read_vols(Vs);
            mT = max(max(max(vol)));
            n1 = power_calc(mT,u,df+1);
            
	    % counting suprathreshold voxels          
	    vol = vol(:); vol(vol == 0 | isnan(vol)) = [];
	    count1 = sum(vol >= u);

            p = sum(MaxT(:) > mT) ./ length(MaxT(:));
            
            %fprintf(1,'%s\n',pwd)
            [dummy mywd] = fileparts(pwd);
            fprintf(1,'%s\t%3.0f\t%6.0f\t%2.6f\t%3.2f\t%3.2f\t%3.0f\t%3.4f\t%3.0f\t%3.0f\t\n',mywd,df,S,p,u,k,vFWHM(1),mT,n1,count1)
            %fprintf(1,'\n')
            
        end
          
        
        
        cd ..
        
    end
    
end
 diary off
 
 return
 