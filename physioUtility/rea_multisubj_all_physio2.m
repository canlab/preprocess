% load('../analysis/onsets');
% load('physio');

condition_names = {'AntOnly_Cue' 'AntStim_Cue' 'StimOnly_Cue'};
ni_codes = {'rea15' 'rea16' 'rea17' 'rea18' 'rea19' 'rea20' 'rea21'};
scanspersess = [192 196 196 184 190 192];     % images per session; length is num sessions
tp = 30;

DX = [];
bgsr = [];
bcengsr = [];
bhr = [];
bcenhr = [];

%setup deconv model
[physio, physio_onsets] = setup_deconv(physio, onsets, condition_names);

% convert pulse to hr
for i=1:length(ni_codes)
    [hr hrs] = physio_pulse2hr(physio.(ni_codes{i}).task_data.pulse);
    physio.(ni_codes{i}).task_data.hr = hr;
    physio.(ni_codes{i}).task_data.hrs = hrs;
end

% deconv
for i=1:length(ni_codes)
    ons = {};
    for j=1:length(condition_names)
        ons{end+1} = physio_onsets.(ni_codes{i}).(condition_names{j})';
    end
    [DX(:,:,i),delta] = physio_make_DX(ons,scanspersess,tp);
    
    [bgsr(:,:,i),bcengsr(:,:,i),gsr,gsrf] = gsr_deconv(physio.(ni_codes{i}).task_data.gsr,scanspersess,DX(:,:,i),tp);
    [bhr(:,:,i),bcenhr(:,:,i),hr,hrf] = hr_deconv(physio.(ni_codes{i}).task_data.hrs,scanspersess,DX(:,:,i),tp);
end

% save results
save('physio_onsets2', 'physio', 'physio_onsets');
save('physio_model2', 'bgsr', 'bcengsr', 'bhr', 'bcenhr');