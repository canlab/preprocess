clear

TOR.multivariate = 'n';
        % Multivariate analysis?
        % choices: y/n
        
TOR.destype = 99;
        % Basic stats module
        % 1 = 1 sample t-test
        % 2 = 2 sample t-test
        % 3 = paired t-test
        % 99 = 2 sample t-test, no intercept - for rfx conjunction
        % etc.
        
TOR.npairs = 11;
        % for one-sample t-test, # of images
        % for 2-sample t-test, total number images
        % for paired t-test, # of pairs   
        % for correlation, length of covariate

% ----- for 2-sample t-test    -----
% TOR.groupvec defined automatically below based on image lists A and B
%TOR.groupvec = [1 1 1 1 1 1 2 2 2 2 2];
% ---------------------------------- 

        
% ----- for covariate analysis -----
TOR.covariate = [1 2 3 4 5 6 7]';
        % for option 5, regression, this is the covariate vector        
        % it must be a column vector
        
TOR.covname = 'int_pure_pos';
        % name of covariate analysis
% ----------------------------------  
                
TOR.gmscale = 9;
        % grand mean scaling--goes to sGMsca
        % 9 = none
        % 1 = 'scaling of overall grand mean'
		% 2 = 'scaling of sF1 grand means'
		% 3 = 'scaling of sF2 grand means'
        % 4 = 'scaling of sF3 grand means'			
		% 5 = 'scaling of sF4 grand means'
		% 6 = 'scaling of sF2 (within sF4) grand means'
		% 7 = 'scaling of sF3 (within sF4) grand means'
		% 8 = '(implicit in PropSca global normalisation)'

        
 % for paired/2-sample t-test use both, for one-sample use A
 A = getfiles2('../../../../imgs/*con_0003.img');
 B = getfiles2('../../../../imgs/*con_0004.img');
      
 TOR.Tmask = -Inf;
    %threshold masking--goes to D.M_.T and M_T
    % -Inf = none
    % 0 = absolute
    % 0.8*sqrt(-1) = proportional
        
 TOR.expMask = 0;
    %explicit mask?--goes to D.M_.X
    % 1 = yes
    % 0 = no
    

 TOR.MaskImg = {'011205tk_con_0015.img'};
    %explicit mask to use--goes to M_P
    %type in path
    % e.g., TOR.MaskImg = {'mask.img'};
    % use {} if no mask to be used
    
 TOR.gcalc = 1;
    %global calculation--goes to sGXcalc
    % 1 = none
    % 2 = user specified
    % 3 = mean voxel value (within per image fullmean/8 mask)
    
 TOR.estnow = 1;
    % estimate now; 1 = yes, 0 = no
   
% -----------------------------------------------------------------------------  
% * End User Input
% -----------------------------------------------------------------------------
clear C;   
        if TOR.destype == 3
            if size(A,1) ~= TOR.npairs, error('Num of subjects mismatch!'), end
            % interleave
            index = 1;
            for i = 1:size(A,1)
                C{index} = A(i,:);
                
                C{index} = B(i,:);
                index = index + 1;
            end
            %D = str2mat(C)
        elseif TOR.destype == 2 |  TOR.destype == 99
            % stack next to each other
            index = 1;
            for i = 1:size(A,1)
                C{index} = A(i,:);
                TOR.groupvec(index) = 1;
                index = index + 1;
            end
            for i = 1:size(B,1)
                C{index} = B(i,:);
                TOR.groupvec(index) = 2;
                index = index + 1;
            end
        else
            %if TOR.destype ~= 3 or 2
            index = 1;
            for i = 1:size(A,1)
                C{index} = A(i,:);
                index = index + 1;
            end
        end
        C = C';     % P must have images in cell array rows
        TOR.P = C;
        % file names
       
        % setup for 2-sample no intercept, from Andrew Holmes
if TOR.destype == 99
    TOR.destype = 2;
            D = struct(...
          'DesName','Two sample t-test (no constant term)',...
          'n',    [Inf 2 1 1],    'sF',{{'obs','group','',''}},...
          'Hform',                'I(:,2),''-'',''group''',...
          'Bform',                '[]',...
          'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
          'iGXcalc',[-1,2,3],'iGMsca',[1,-9],'GM',[],...
          'iGloNorm',9,'iGC',12,...
          'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
          'b',struct('aTime',1))

    tor_spm_spm_ui(TOR,D)
    
else
    tor_spm_spm_ui(TOR)
end