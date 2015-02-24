function [b,k] = ridge(y,X,k)
%RIDGE Ridge regression.
%   B = RIDGE(Y,X,K) returns the vector of regression coefficients, B.
%   Given the linear model: y = Xb, 
%   X is an n by p matrix, y is the n by 1 vector of observations and k
%   is a scalar constant.
%
%	B is scaled - so it may be different than the b you expected, unless
%		your model X is scaled to have std = 1.
%
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.5 $  $Date: 1997/11/29 01:46:39 $
%
% 	Modified by Tor Wager, 1/10/01, to show the df-trace criterion
%		for unknown k (Tripp, 1983; Myers, 1986/90, 2nd Ed.)
%
%	If k is not a scalar, runs df-trace to choose k.

if  nargin < 3,              
    error('Requires at least three input arguments.');      
end 

% Make sure k is a scalar if not, do the k-choosing algorithm.
if prod(size(k)) ~= 1
   	disp('Choosing k with df-trace method.');
	disp('	...setting up xtx, identity matrix')
    Xtmp = X; Xtmp(:,sum(Xtmp) == size(Xtmp,1)) = [];
 	for i = 1:size(Xtmp,2)
		   Xc(:,i) = (X(:,i) - mean(X(:,i))) / std(X(:,i));		% centered X
	end
	xtx = full(Xc' * Xc);
	ident = eye(size(xtx,1));
	
	disp('	...running for selected k')
	index = 1;
    go = 1;
    
	x = 0;
  %while x < 2 & go
  for k = .001:.001:1.5
  %k = x;
		xtxi = inv(xtx + k*ident);
	   traceh(index) = trace(Xc * xtxi * Xc');
	   vif(index) = trace(xtxi * xtx * xtxi)/size(xtx,1);
	   disp(['		...k = 	' num2str(k) ',	Hk = ' num2str(traceh(index)) ',	vif = ' num2str(vif(index))])
	   index = index + 1;
       
       %if index == 3, crit = (traceh(2) - traceh(1)) .* .02;,end % 2% change
       %if index > 8
       %    if traceh(index-1) - mean(traceh(3:7)) < crit, go = 0;,end
       %end
       
       %x = x + .001;
       
  	end
	b.traceh = traceh; b.vif = vif;
	figure;subplot 121; plot(x,traceh); title('traceh');subplot 122; plot(x,vif);title('vif')
    pause(4); close
    k = x(traceh == max(traceh))
    
	return
end

% Check that matrix (X) and left hand side (y) have compatible dimensions
[n,p] = size(X);

[n1,collhs] = size(y);
if n~=n1, 
    error('The number of rows in Y must equal the number of rows in X.'); 
end 

if collhs ~= 1, 
    error('Y must be a column vector.'); 
end

% Normalize the columns of X to mean, zero, and standard deviation, one.
mx = mean(X);
stdx = std(X);
idx = find(abs(stdx) < sqrt(eps)); 
if any(idx)
  stdx(idx) = 1;
end

MX = mx(ones(n,1),:);
STDX = stdx(ones(n,1),:);
Z = (X - MX) ./ STDX;
if any(idx)
  Z(:,idx) = ones(n,length(idx));
end

% Create matrix of pseudo-observations and add to the Z and y data.
pseudo = sqrt(k) * eye(p);
Zplus  = [Z;pseudo];
yplus  = [y;zeros(p,1)];

% Solve the linear system for regression coefficients.
b = Zplus\yplus;