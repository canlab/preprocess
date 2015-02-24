function humps = findhumps(signal, firstval, humpmin)
    dsig = gradient(signal);
    %ddsig = gradient(dsig);
    zero_inds = find_intercept(dsig);
    
    humpmin = humpmin * (max(signal)-min(signal));
    %
    %max_inds = zero_inds(find(ddsig(zero_inds)<0));
    %
    max_inds = [];
    
    if ~isempty(zero_inds)
        if zero_inds(1) ~= 1
            zero_inds = [1; zero_inds];
        end;
        if zero_inds(length(zero_inds))~=length(zero_inds)
            zero_inds = [zero_inds; length(zero_inds)];
        end;
    end;
    
    for i = 2:length(zero_inds)-1
        if (signal(zero_inds(i))-signal(zero_inds(i-1))>humpmin) || (signal(zero_inds(i))-signal(zero_inds(i+1))>humpmin)
            max_inds = [max_inds; zero_inds(i)];
        end;
    end;
    
    max_vals = signal(max_inds) - signal(1);
    
    humps = [max_inds+firstval-1 max_vals max_inds];
end
    




