function tor_spm_random_effects(BAT)

 % MODEL SPEC OPTIONS
 
% see tor_basic_models.m for explanations
BAT.multivariate = 'n';
BAT.destype = 1;
        % Basic stats module
        % 1 = 1 sample t-test
        % 2 = 2 sample t-test
        % 3 = paired t-test
        % etc.
%BAT.npairs = 11;
        % for one-sample t-test, # of images
        % for paired t-test, # of pairs             
BAT.gmscale = 9;
        % grand mean scaling--goes to sGMsca
        % 9 = none
        
 % specified in biman5_run_randfx
 % for paired t-test use both, for one-sample use A
 % A = getfiles2('*con_0003.img');
 % B = getfiles2('*con_0004.img');
      
 BAT.Tmask = -Inf;
    %threshold masking--goes to D.M_.T and M_T
    % -Inf = none
    % 0 = absolute
    % 0.8*sqrt(-1) = proportional
        
 if ~isfield(BAT,'expMask'),BAT.expMask = 0;,end
    %explicit mask?--goes to D.M_.X
    % 1 = yes
    % 0 = no
    
 % specified in biman5_run_randfx
 % BAT.MaskImg = {'011205tk_con_0015.img'};
    %explicit mask to use--goes to M_P
    %type in path
    % e.g., BAT.MaskImg = {'mask.img'};
    % use {} if no mask to be used
    
 BAT.gcalc = 1;
    %global calculation--goes to sGXcalc
    % 1 = none
    % 2 = user specified
    % 3 = mean voxel value (within per image fullmean/8 mask)
    
 %BAT.estnow = 1;
    % estimate now; 1 = yes, 0 = no
   
 % RESULTS OPTIONS
 
 BAT.maskWOtherContrasts = 0;
	% mask with other contrasts; 1 or 0

%BAT.multCompCorrect = 'FWE'; 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

%BAT.u = .05;
	% height threshold - T or p value

%BAT.k = 000;
	% extent threshold - number of voxels

if isfield(BAT,'Ic')
    if ismepty(BAT.Ic)
        BAT.Ic = 2;
    end
else, BAT.Ic = 2;
end
	% optional: which contrast to show in results (number)
	% or 'all'.  If 'all', uses all contrasts, in order.
    % Ic = 2 is the default for random effects 
    % (1 is "effects of interest", 2 is first contrast)

% BAT.resdir = pwd;   
    
    
    
    
% -----------------------------------------------------------------------------  
% * End User Input
% -----------------------------------------------------------------------------
        %clear C
        %if BAT.destype == 3
        %    if size(A,1) ~= BAT.npairs, error('Num of subjects mismatch!'), end
            % interleave
        %    index = 1;
        %    for i = 1:size(A,1)
        %        C{index} = A(i,:);
                
        %        C{index} = B(i,:);
        %        index = index + 1;
        %    end
            %D = str2mat(C)
            %else
            %if BAT.destype ~= 3
            
            index = 1;
            for i = 1:size(BAT.P,1)
                C{index} = deblank(BAT.P(i,:));
                index = index + 1;
            end
            %end
        C = C';     % P must have images in cell array rows
        BAT.P = C;
        % file names
       
if BAT.estnow    
    tor_spm_spm_ui(BAT)
end

if BAT.estcons
% - - - - - - - - - - - - - - - - - - - - - - - - -
% results
% - - - - - - - - - - - - - - - - - - - - - - - - -
clear c

cnum = 1;
c{cnum} = [1];
names{cnum} = pwd;
type{cnum} = 'T';

estimateContrasts(c,names,type,0);

end

if BAT.dores
    try
        [hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(BAT);
        spm_list('List',SPM,VOL,[],[],'',hReg);	
        print -dpsc2 -painters -append -noui ../rfx_results.ps

        % write results mask 
        mask = voxel2mask(SPM.XYZ',VOL.DIM');
        V = spm_vol('mask.img');   
        V.fname = 'res_mask.img';
        V = spm_write_vol(V,mask); 

        % write filtered and clusters
        V=spm_vol('spmT_0002.img');V2=spm_vol('res_mask.img');Vo=V;Vo.fname='spmT_filtered_0002.img';
        Vo=spm_imcalc([V,V2],Vo,'i1 .* i2');
        [clusters,clum]=mask2clusters('spmT_filtered_0002.img');
        save spmT0002_clusters clusters
        
    catch
         disp('Cannot make results -- update for SPM2?');
    end
end
            
            
return
