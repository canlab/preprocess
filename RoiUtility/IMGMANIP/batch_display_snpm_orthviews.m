function batch_display_snpm_orthviews
%
% start in the directory above individual snpm* analysis directories
% Tor Wager

D = dir;

P      = spm_get(1,'.img','please select canonical image',[],0);

for i = 3:length(D)
    
    myname = [D(i).name filesep D(i).name '_clusters.mat'];
    if exist(myname) == 2, 
        eval(['load ' myname])
        spm_image('init',P);
        spm_orthviews('AddColouredBlobs',1,CLU.XYZ,CLU.Z,CLU.M,[0 0 1]);
        disp(['This image is from ' myname])
        input('Press <RETURN> to continue.')
    end
    
end