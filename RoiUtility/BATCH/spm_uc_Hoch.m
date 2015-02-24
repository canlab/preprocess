function [u,Ps,Ts] = spm_uc_Hoch(q,df,STAT,n,Vs,Vm)
% False Discovery critical height threshold
% FORMAT [u,Ps,Ts] = spm_uc_Hoch(q,df,STAT,n,Vs[,Vm])
%
% q     - critical expected False Discovery Rate
% df    - [df{interest} df{residuals}]
% STAT  - Statisical feild  (see comments below about FWER and EHoch)
%		'Z' - Gaussian feild
%		'T' - T - feild
%		'X' - Chi squared feild
%		'F' - F - feild
% n     - number of component SPMs in conjunction
% Vs    - Mapped statistic image(s)
%          -or-
%         Vector of sorted p-values (saves i/o w/ repeated calls)
% Vm    - Mask in 1 of 3 forms
%           o Scalar, indicating implicit mask value in statistic image(s)
%           o Vector of indicies of elements within mask
%           o Mapped mask image
%
% u     - critical height
% Ps    - Sorted p-values
% Ts    - Sorted statistic values
%
%___________________________________________________________________________
%
% The Benjamini & Hochberch (1995) False Discovery Rate (Hoch) procedure
% finds a threshold u such that the expected Hoch is at most q.
% spm_uc_Hoch returns this critical threshold u. 
%
% For repeated use of a given statistic image, return Ps in the place
% of Vs:
%      [P Ps] = spm_uc_Hoch(Z1,df,STAT,n,Vs);  %-Initialization, image read
%      P      = spm_uc_Hoch(Z2,df,STAT,n,Ps);  %-Image not read, Ps used
%
% Note that a threshold of Inf is possible if q is very small.  This
% means that there is no threshold such that Hoch is controlled at q.
%
%
% Background
%
% For a given threshold on a statistic image, the False Discovery Rate
% is the proportion of suprathreshold voxels which are false positives.
% Recall that the thresholding of each voxel consists of a hypothesis
% test, where the null hypothesis is rejected if the statistic is larger
% than threshold.  In this terminology, the Hoch is the proportion of
% rejected tests where the null hypothesis is actually true.
%
% A Hoch proceedure produces a threshold that controls the expected Hoch
% at or below q.  The Hoch adjusted p-value for a voxel is the smallest q
% such that the voxel would be suprathreshold.
%
% In comparison, a traditional multiple comparisons proceedure
% (e.g. Bonferroni or random field methods) controls Familywise Error
% rate (FWER) at or below alpha.  FWER is the chance of one or more
% false positives *anywhere* (not just among suprathreshold voxels).  A
% FWER adjusted p-value for a voxel is the smallest alpha such that the
% voxel would be suprathreshold.
%
% 
% If there is truely no signal in the image anywhere, then a Hoch
% proceedure controls FWER, just as Bonferroni and random field methods
% do. (Precisely, controlling E(Hoch) yeilds weak control of FWE).  If
% there *is* some signal in the image, a Hoch method will be more powerful
% than a traditional method.
%
%
% References
%
% Benjamini & Hochberg (1995), "Controlling the False Discovery Rate: A
% Practical and Powerful Approach to Multiple Testing". J Royal Stat Soc,
% Ser B.  57:289-300.
%
% Benjamini & Yekutieli (2001), "The Control of the false discovery rate
% in multiple testing under dependency". To appear, Annals of Statistics.
% Available at http://www.math.tau.ac.il/~benja 
%___________________________________________________________________________
% @(#)spm_uc_Hoch.m	2.2 Thomas Nichols 01/08/07

if (nargin<6), Vm = []; end

% Set Benjamini & Yeuketeli cV for independence/PosRegDep case
%-----------------------------------------------------------------------
cV = 1; 


% Load, mask & sort statistic image (if needed)
%-----------------------------------------------------------------------
if isstruct(Vs)
  if (n ~= length(Vs))
    error(sprintf('n & number of mapped images doesn''t match (%d,%d)',...
		  n,length(Vs)));
  end
  Ts = spm_read_vols(Vs(1));
  for i = 2:n
    Ts = min(Ts,spm_read_vols(Vs(i)));
  end
  if ~isempty(Vm)
    if isstruct(Vm)
      Ts(spm_read_vols(Vm)==0) = [];
    elseif (prod(size(Vm))==1)
      Ts(Ts==Vm) = [];
    else 
      Ts = Ts(Vm);
    end
  end
  Ts(isnan(Ts)) = [];
  Ts = flipud(sort(Ts(:)));
end


% Calculate p values of image (if needed)
%-----------------------------------------------------------------------
if isstruct(Vs)
  if      STAT == 'Z'
    Ps = (1-spm_Ncdf(Ts)).^n;
  elseif  STAT == 'T'
    Ps = (1-spm_Tcdf(Ts,df(2))).^n;
  elseif  STAT == 'X'
    Ps = (1-spm_Xcdf(Ts,df(2))).^n;
  elseif  STAT == 'F'
    Ps = (1-spm_Fcdf(Ts,df)).^n;
  end
else
  Ps = Vs;
end

S = length(Ps);


% Calculate Hoch inequality RHS
%-----------------------------------------------------------------------
Fi  = 1./(S:-1:1)'*q;


% Find threshold
%-----------------------------------------------------------------------
I   = max(find(Ps<=Fi));
if isempty(I)
  u = Inf;
else
  u   = Ts(I);
end



