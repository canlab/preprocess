function int = interaction(a,b)
% function int = interaction(a,b)
% 
% calculates centered interaction term for predictors a and b


int = scale(scale(a,1) .* scale(b,1),1);

return