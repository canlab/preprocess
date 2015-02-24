function intercepts = find_intercept(signal)
    n = length(signal);
    possints = find(((sign(signal(2:n))~=sign(signal(1:n-1))) .* (signal(2:n)~=0) .* (signal(1:n-1)~=0)) + (signal(1:n-1) == 0));
    
    
    
    for i = 1:length(possints)
        if possints(i) ~= n && (abs(signal(possints(i))) > abs(signal(possints(i)+1)));
            possints(i) = possints(i)+1;
        end;
    end;
    
    intercepts = possints;
    
end
