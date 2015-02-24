function V = spmVol(img)
    if(iscellstr(img))
        V = spm_vol(char(img)); % to avoid SPM's behavior of returning cells if you pass in a cellstr
    else
        V = spm_vol(img);
    end
end


