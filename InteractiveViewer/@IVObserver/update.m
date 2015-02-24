function update(obs, mmPos, voxPos)
    if(~isempty(obs.updateFun))
        obs.updateFun(obs.iv, mmPos, voxPos);
    else
        fprintf('No update function\n');
    end
    %plot(sin(0:.01:3));
end