function outsignal = smoothy(insignal, smoothwin)
    if smoothwin/2~=floor(smoothwin/2), smoothwin = smoothwin + 1; end;
    w = gausswin(smoothwin);
    outsignal = conv(insignal,w);
    outsignal = outsignal((smoothwin/2):(length(outsignal)-smoothwin/2))./sum(w);
end
