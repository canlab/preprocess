function [u,Ps,Ts] = spm_uc_Sidak(q,df,STAT,n,Vs,Vm)
% False Discovery critical height threshold
% FORMAT [u,Ps,Ts] = spm_uc_Sidak(q,df,STAT,n,Vs[,Vm])
%
% q     - critical expected False Discovery Rate
% df    - [df{interest} df{residuals}]
% STAT  - Statisical feild  (see comments below about FWER and ESidak)
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
% spm_uc_Sidak returns a FWER-controlling critical threshold u which can be
% more sensitive than Bonferroni.  It makes an assumption of positive
% correlation, or positive regression dependency on subsets.
%
% For repeated use of a given statistic image, return Ps in the place
% of Vs:
%      [P Ps] = spm_uc_Sidak(Z1,df,STAT,n,Vs);  %-Initialization, image read
%      P      = spm_uc_Sidak(Z2,df,STAT,n,Ps);  %-Image not read, Ps used
%
% Note that a threshold of Inf is possible if q is very small.  This
% means that there is no threshold such that FWER is controlled at q.
%
%
% References
%
% Sidakberg (),  ****
%
% Benjamini & Yekutieli (2001), "The Control of the false discovery rate
% in multiple testing under dependency". Annals of Statistics. ***
%___________________________________________________________________________
% Based on spm_uc_Hoch.m, ver 1.2
% @(#)spm_uc_Sidak.m	1.1 03/04/17

if (nargin<6), Vm = []; end


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


% Calculate Sidak inequality RHS
%-----------------------------------------------------------------------
Fi  = 1-(1-q).^(1./(S:-1:1)');


% Find threshold
%-----------------------------------------------------------------------
I   = min(find(Ps>Fi));  % This is one above the threshold
if isempty(I)
  u = Ps(end);  % Very suspicious; everything significant!
elseif I==1
  u = Inf;
else
  I   = I-1;  
  u   = Ts(I);
end


