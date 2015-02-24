function [t,p] = tor_bonf(p,numsearch,df)

% bonferroni corrected critical t value

p = 1 - (p ./ numsearch);

t = tinv(p,df);

p = 1-p;

return

