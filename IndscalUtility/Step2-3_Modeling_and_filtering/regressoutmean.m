function DATA=regressoutmean(DATA)

% old cell data format
%numsub=size(DATA.dat,1);
%numreg=size(DATA.dat,2);

% newer vector of cells for each subject format
numsub=length(DATA.dat);    %number of subjects
numreg=size(DATA.dat{1},2);  %number of ROIs


for s = 1:numsub
    
    clear mdat;clear krr;clear gmm;
    
    dat = DATA.resids{s};
    gm  = mean(dat')';     % mean across regions
    gm(:,end+1) = 1; 
    
    gb = pinv(gm) * dat;    % global betas 
    gf = gm * gb;           % global fit
    rr = dat - gf;          % remove global fit
    
    DATA.resids{s} = rr;
    DATA.resids_descrip = 'Global mean across regions removed by regression.';
    
end


return




    
    % old cell data format
    
    
   for sub = 1:numsub
    
    clear mdat;clear krr;clear gmm;
    
       for roi=1:numreg
       mdat(roi,:,:)=DATA.r{sub,roi}';      
       end
       gmm=squeeze(mean(mdat));     %gmm is blocks*scans; mean across regions for that subject
       numblock=size(gmm,1);        %number of blocks for that subject
       for k=1:numblock;
       clear gm;
       gm = squeeze(gmm(k,:))';
       gm(:,end+1) = 1; %DX with 2 regressors: mean (across regions) and intercept
       dat=squeeze(mdat(:,k,:))';
        gb = pinv(gm) * dat;    %global betas for each block
        gf = gm * gb;           %global fit
        rr = dat - gf;          %remove global fit
        krr(:,k,:)=rr';
       end      
       for roi=1:numreg
       DATA.r{sub,roi}=squeeze(krr(roi,:,:))';
       end
end
