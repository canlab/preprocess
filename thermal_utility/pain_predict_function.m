function r = pain_predict_function(p, x)
    
    % constrain p(3), exponent, to be >= 1 ?
    
    if x - p(1) - 32 < 0
        r = 0;
        
    else
       
       r =  p(2).*(x - p(1) - 32) .^ p(3);
       
    end
    
end

