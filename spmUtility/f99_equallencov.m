function A = f99_equallencov(type,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20)
%------------------------------------------------------------------------------
% We wncountered the error when covariates are not all of equal length
% e.g. [cov1, cov2, cov3] then returns an error.
% You can use f96_equallencov(cov1,cov2,cov3) 
% this function generates a matrix of the covariates with the end
% filled with zeros 
% these zeros should be ignored by fmri_stat.
% NOTE : only 20 covariates can be specified
%-----------------------------------------------------------------------------
ncov = nargin-1;

lenmax = 0;
for i = 1 : ncov
  str = ['lenmax = max(lenmax, length(c' mat2str(i) '));'];
  eval(str)
end

A = zeros(ncov,lenmax);
if type
    A = (A+1).*type;
end

for i = 1 : ncov
  str = ['A(i,1:size(c' mat2str(i) ',2)) = c' mat2str(i) ';'];
  eval(str)
end
