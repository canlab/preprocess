% update_spm_vol_structs(SPMs)
%
% Goes through the spm_vol-related fields of an SPM and reloads them
%
% Useful when your contrasts (or other images) appear to be gibberish because 
% the SPM.mat file was run on a machine with a different endianness, or when you
% set the origin *after* running your analysis.
%
% NB: not legitimate if you've chanegd the voxel size.

function update_spm_vol_structs(SPMs)
    SPMs = cellstr(SPMs);
    orig_wd = pwd();

    try
        SPM_fields = {'Vbeta' 'VHp' 'VResMS' 'VM' 'xVol.VRpv'};
        xCon_fields = {'Vcon' 'Vspm'};

        for i = 1:length(SPMs)
            load(SPMs{i});
            cd(SPM.swd);
            
            breakpoints = dbstatus();
            dbclear if caught error
            for j=1:length(SPM_fields)
                try
                    % NB: eval() is usually a bad idea, but Matlab's inability to dynamically string together field names
                    % precludes a less painful option
                    eval(sprintf('SPM.%s = spm_vol(strvcat(SPM.%s.fname));', SPM_fields{j}, SPM_fields{j}));
                catch
                    fprintf('Skipping missing field %s\n', SPM_fields{j});
                end
            end
            
            SPM.xVol.M = SPM.xVol.VRpv.mat;
            SPM.xVol.iM = pinv(SPM.xVol.M);

            for j=1:length(xCon_fields)
                for k=1:length(SPM.xCon)
                    try
                        SPM.xCon(k).(xCon_fields{j}) = spm_vol(SPM.xCon(k).(xCon_fields{j}).fname);
                    catch
                        fprintf('Skipping missing field SPM.xCon(%d).%s\n', k, xCon_fields{j});
                    end
                end
            end
            dbstop(breakpoints);

            save SPM SPM
            cd(orig_wd);
        end
    catch
        cd(orig_wd);
        error(lasterror());
    end
end